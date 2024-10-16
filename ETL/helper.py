import click
import psycopg2
import traceback
import typing as t
from pathlib import Path
from tqdm import tqdm
import json
import os

FileDescriptor = t.Union[int, str]

def copy_from_tsv(conn: 'psycopg2.connection', table: str, columns: list[str], r: FileDescriptor):
  ''' Copy from a file descriptor into a postgres database table through as psycopg2 connection object
  :param con: The psycopg2.connect object
  :param r: An file descriptor to be opened in read mode
  :param table: The table top copy into
  :param columns: The columns being copied
  '''
  import os
  with conn.cursor() as cur:
    with os.fdopen(r, 'rb', buffering=0, closefd=True) as fr:
      columns = fr.readline().strip().split(b'\t')
      cur.copy_expert(
        sql=f'''
        COPY {table} ({",".join(f'"{c.decode()}"' for c in columns)})
        FROM STDIN WITH CSV DELIMITER E'\\t'
        ''',
        file=fr,
      )
    conn.commit()

def copy_from_records(conn: 'psycopg2.connection', table: str, columns: list[str], records: t.Iterable[dict]):
  ''' Copy from records into a postgres database table through as psycopg2 connection object.
  This is done by constructing a unix pipe, writing the records with csv writer
   into the pipe while loading from the pipe into postgres at the same time.
  :param con: The psycopg2.connect object
  :param table: The table to write the pandas dataframe into
  :param columns: The columns being written into the table
  :param records: An iterable of records to write
  '''
  import os, csv, threading
  r, w = os.pipe()
  # we copy_from_tsv with the read end of this pipe in
  #  another thread
  rt = threading.Thread(
    target=copy_from_tsv,
    args=(conn, table, columns, r,),
  )
  rt.start()
  try:
    # we write to the write end of this pipe in this thread
    with os.fdopen(w, 'w', closefd=True) as fw:
      writer = csv.DictWriter(fw, fieldnames=columns, delimiter='\t')
      writer.writeheader()
      writer.writerows(records)
  finally:
    # we wait for the copy_from_tsv thread to finish
    rt.join()

def import_gene_set_library(
  plpy,
  library: Path | str,
  species: str = 'human',
  prefix='',
  postfix=''
):
  import re
  import json
  import uuid

  # fetch the gene_set, gather the background genes
  new_gene_sets = []
  background_genes = set()
  n_geneset_genes = 0
  with Path(library).open('r') as fr:
    for line in tqdm(fr, desc='Loading gmt...'):
      line_split = line.strip().split('\t')
      if len(line_split) < 3: continue
      term, _description, *raw_genes = line_split
      genes = [
        cleaned_gene
        for raw_gene in map(str.strip, raw_genes)
        if raw_gene
        for cleaned_gene in (re.split(r'[;,:\s]', raw_gene)[0],)
        if cleaned_gene
      ]
      new_gene_sets.append(dict(
        term=prefix+term+postfix,
        genes=genes,
      ))
      background_genes.update(genes)
      n_geneset_genes += len(genes)

  # get a mapping from background_genes to background_gene_ids
  gene_map, = plpy.cursor(
    plpy.prepare(
      '''
        select coalesce(jsonb_object_agg(g.gene, g.gene_id), '{}') as gene_map
        from app_public_v2.gene_map($1) as g
      ''',
      ['varchar[]']
    ),
    [list(background_genes)]
  )
  gene_map = json.loads(gene_map['gene_map'])

  # upsert any new genes not in the mapping & add them to the mapping
  new_genes = {
    id: dict(id=id, symbol=gene)
    for gene in tqdm(background_genes - gene_map.keys(), desc='Preparing new genes...')
    for id in (str(uuid.uuid4()),)
  }
  if new_genes:
    copy_from_records(
      plpy.conn, 'app_public_v2.gene', ('id', 'symbol',),
      tqdm(new_genes.values(), desc='Inserting new genes...'))
    gene_map.update({
      new_gene['symbol']: new_gene['id']
      for new_gene in new_genes.values()
    })

  existing_terms = {
    row['term']
    for row in plpy.cursor('select term from app_public_v2.gene_set', tuple())
  }

  copy_from_records(
    plpy.conn, 'app_public_v2.gene_set', ('term', 'gene_ids', 'n_gene_ids', 'species',),
    tqdm((
      dict(
        term=gene_set['term'],
        gene_ids=json.dumps({gene_map[gene]: None for gene in gene_set['genes']}),
        n_gene_ids=len(gene_set['genes']),
        species=species,
      )
      for gene_set in new_gene_sets
      if gene_set['term'] not in existing_terms
    ),
    total=len(new_gene_sets) - len(existing_terms),
    desc='Inserting new genesets...'),
  )

def import_pb_info(plpy):
  import requests
  import re
  from datetime import datetime
  import ast

  to_ingest = [
    r['pmid']
    for r in plpy.cursor(
      f'''
        select pmid
        from app_public_v2.gse_info
        where gse not in (
          select pmid
          from app_public_v2.pmid_info
        )
      '''
    )
  ]

  to_ingest = list(set(to_ingest))

  to_ingest = [id for id in to_ingest if id is not None]

  # deal with multiple pmids in one entry
  formated = []
  for val in to_ingest:
    if '[' in val:
      to_ingest.remove(val)
      parsed_data = ast.literal_eval(val)
      formated += parsed_data
    else: formated.append(val)
  to_ingest = formated

  if os.path.exists('data/pb_info_to_ingest.json'):
    with open('data/pb_info_to_ingest.json') as f:
      pb_info = json.load(f)
    to_pull = list(set(to_ingest).difference(list(pb_info.keys())))
  else:
    pb_info = {}
    to_pull = to_ingest

  for i in tqdm(range(0, len(to_pull), 250), 'Pulling titles...'):
    for j in range(10):
      try:
        ids_string = ",".join(to_pull[i:i+250])
        res = requests.get(f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&retmode=json&id={ids_string}')
        ids_info = res.json()
      except KeyboardInterrupt:
        raise
      except Exception as e:
        traceback.print_exc()
        print('Error resolving info. Retrying...')
        continue
      for id in ids_info['result']['uids']:
        pb_info[id] = {}
        try:
          pb_info[id]['title'] = ids_info['result'][id]['title']
        except KeyError:
          pass
        try:
          possible_formats = ['%Y %b', '%Y %b %d']
          for format_str in possible_formats:
              try:
                  date_object = datetime.strptime(ids_info['result'][id]['pubdate'], format_str)
              except ValueError:
                  pass
          formatted_date = date_object.strftime("%Y-%m")
          pb_info[id]['date'] = formatted_date
        except ValueError:
          pass
        try:
          pb_info[id]['doi']  = next((item['value'] for item in ids_info['result'][id]['articleids'] if item['idtype'] == 'doi'), None)
        except KeyError:
          pass
        try:
          pb_info[id]['pmcid']  = next((item['value'] for item in ids_info['result'][id]['articleids'] if item['idtype'] == 'pmc'), None)
        except KeyError:
          pass  
      break


  with open('data/pb_info_to_ingest.json', 'w') as f:
    json.dump(pb_info, f)

  already_ingested = [
    r['pmid']
    for r in plpy.cursor(
      f'''
        select pmid
        from app_public_v2.pmid_info;
      '''
    )
  ]

  for pmid in tqdm(to_ingest):
    if pmid in pb_info and pmid not in already_ingested:
      plpy.execute(
        '''
          insert into app_public_v2.pmid_info (pmid, pmcid, title, pub_date, doi)
          values (%s, %s, %s, %s, %s)
          on conflict (pmid) do nothing;
        ''',
        (pmid, pb_info[pmid]['pmcid'], pb_info[pmid]['title'], pb_info[pmid]['date'], pb_info[pmid]['doi'])
      )

def import_gse_attrs(plpy, species='human'):
  import json
  from tqdm import tqdm

  to_ingest = [
    r['gse']
    for r in plpy.cursor(
      f'''
        select gse
        from app_public_v2.gse_info
        where species = '{species}'
          and gse not in (
            select gse
            from app_public_v2.gse_terms
            where species = '{species}'
          )
      '''
    )
  ]

  to_ingest = list(set(to_ingest))
  
  with open(f'data/keyterms_{species}.json') as f:
    gse_attrs = json.load(f)

  for gse in tqdm(to_ingest):
    sql = """
        INSERT INTO app_public_v2.gse_terms (gse, llm_attrs, pubmed_attrs, mesh_attrs, species)
        VALUES (%s, %s, %s, %s, %s)
        ON CONFLICT (gse, species) DO UPDATE SET
        llm_attrs = EXCLUDED.llm_attrs,
        pubmed_attrs = EXCLUDED.pubmed_attrs,
        mesh_attrs = EXCLUDED.mesh_attrs;
        """
    if gse not in gse_attrs:
      plpy.execute(sql, (gse, [], [], [], species))
    else:
      llm_attrs = gse_attrs[gse]
      pubmed_attrs = []
      mesh_attrs = []
      plpy.execute(sql, (gse, llm_attrs, pubmed_attrs, mesh_attrs, species))

  #plpy.execute('refresh materialized view app_public_v2.terms_count_combined', [])

  plpy.execute("""alter table app_public_v2.gse_info drop column gse_attrs cascade;
                  alter table app_public_v2.gse_info add column gse_attrs varchar;
                  UPDATE app_public_v2.gse_info gi
                  SET    gse_attrs = array_to_string(
                          ARRAY(
                            SELECT val
                            FROM   unnest(gt.llm_attrs) AS val
                            WHERE  val IS NOT NULL
                            UNION ALL
                            SELECT val
                            FROM   unnest(gt.mesh_attrs) AS val
                            WHERE  val IS NOT NULL
                            UNION ALL
                            SELECT val
                            FROM   unnest(gt.pubmed_attrs) AS val
                            WHERE  val IS NOT NULL
                          ), ' ')
                  FROM   app_public_v2.gse_terms gt
                  WHERE  gi.gse = gt.gse and gi.species = gt.species;""", [])


def import_term_categories(plpy):
  import pandas as pd
  from tqdm import tqdm

  with open('data/keyterm_categories.json') as f:
    new_cats = json.load(f)

  copy_from_records(
    plpy.conn, 'app_public_v2.term_categories', ('term_name', 'category',),
    (
      dict(term_name=term, category=new_cats[term])
      for term in tqdm(new_cats, desc='Inserting term categories...', total=len(new_cats))
    ),
  )

  plpy.execute('refresh materialized view app_public_v2.terms_count_combined;', [])
  plpy.execute('refresh materialized view app_public_v2.category_total_count;', [])

def import_gse_info(plpy, species='human'):
  import GEOparse
  from itertools import chain
  from tqdm import tqdm
  import pandas as pd
  import json
  from datetime import datetime

  plpy.execute('refresh materialized view app_public_v2.gene_set_gse;')

  # find subset to add info to
  to_ingest = [
    r['gse']
    for r in plpy.cursor(
      f'''
        select gse, species
        from app_public_v2.gene_set_gse
        where gse not in (
          select gse
          from app_public_v2.gse_info
          where gse_info.species = '{species}'
        ) and species = '{species}'
      '''
    )
  ]

  to_ingest = list(set(to_ingest))

  print(f'Found {len(to_ingest)} new GSEs to ingest')
  # fetch processed gse groupings/conditions titles info
  if not os.path.exists(f'data/gse_info_to_ingest_{species}.json'):
    try:
      with open(f'data/{species}-gse-processed-meta.json') as f:
        gse_info = json.load(f)
    except:
      raise RuntimeError('Missing metadata. Please ensure path is correct.')
    
    os.makedirs('data/geo', exist_ok=True)

    gse_info_to_ingest = {}
    samples_to_ingest = set()

    for gse in tqdm(to_ingest, desc='Fetching GSE info...'):
      if ',' in gse:
        gse_split = gse.split(',')[0]
      else:
        gse_split = gse
      for i in range(10):
        try:
          geo_meta = GEOparse.GEOparse.get_GEO(geo=gse_split, silent=True, destdir='data/geo')
        except:
          print(f'Failed to fetch {gse}')
          continue
        break
      gse_info_to_ingest[gse] = {}
      gse_info_to_ingest[gse]['title'] = geo_meta.get_metadata_attribute('title')
      gse_info_to_ingest[gse]['summary'] = geo_meta.get_metadata_attribute('summary')
      try:
          gse_info_to_ingest[gse]['pmid'] = geo_meta.get_metadata_attribute('pubmed_id')
      except:
          gse_info_to_ingest[gse]['pmid'] = None

      date_object = datetime.strptime(geo_meta.get_metadata_attribute('status'), "Public on %b %d %Y")
      formatted_date = date_object.strftime("%Y-%m-%d")
      gse_info_to_ingest[gse]['publication_date'] = formatted_date  # yyyy-mm-dd
      gse_info_to_ingest[gse]['platform'] = geo_meta.get_metadata_attribute('platform_id')
      gse_info_to_ingest[gse]['species'] = species
      gse_info_to_ingest[gse]['sample_groups'] = {"samples": gse_info[gse]['samples'], "titles": gse_info[gse]['titles']}
      gse_info_to_ingest[gse]['silhouette_score'] = gse_info[gse]['silhouette_score']

      samples = list(chain.from_iterable(gse_info[gse]['samples'].values()))
      samples_to_ingest.update(samples)

    
    with open(f'data/gse_info_to_ingest_{species}.json', 'w') as f:
      json.dump(gse_info_to_ingest, f)
    
    with open(f'data/samps_to_ingest_{species}.json', 'w') as f2:
      json.dump(list(samples_to_ingest), f2)
  else:
    with open(f'data/gse_info_to_ingest_{species}.json') as f:
      gse_info_to_ingest = json.load(f)
    
    with open(f'data/samps_to_ingest_{species}.json') as f2:
      samples_to_ingest = json.load(f2)


  gsms_ingested = [
    r['gsm']
    for r in plpy.cursor(
      f'''
        select gsm
        from app_public_v2.gsm_meta
      '''
    )
  ]


  samps_df = pd.read_csv(f'data/gse_gsm_meta_{species}.csv').set_index('gsm')

  samples_to_ingest = list(set(samples_to_ingest) - set(gsms_ingested))
  samps_df_to_ingest = samps_df.loc[list(set(samples_to_ingest).intersection(samps_df.index))]


  copy_from_records(
    plpy.conn, 'app_public_v2.gse_info', ('gse', 'pmid', 'title', 'summary', 'published_date', 'species', 'platform', 'sample_groups', 'silhouette_score'),
    tqdm((
      dict(
        gse=gse,
        pmid=gse_info_to_ingest[gse]['pmid'],
        title=gse_info_to_ingest[gse]['title'],
        summary=gse_info_to_ingest[gse]['summary'],
        published_date=gse_info_to_ingest[gse]['publication_date'],
        species=species,
        platform=gse_info_to_ingest[gse]['platform'],
        sample_groups=json.dumps(gse_info_to_ingest[gse]['sample_groups']),
        silhouette_score=gse_info_to_ingest[gse]['silhouette_score']
      )
      for gse in to_ingest
      if gse in gse_info_to_ingest
    ),
    total=len(to_ingest),
    desc='Inserting GSE info..')
  )

  copy_from_records(
    plpy.conn, 'app_public_v2.gsm_meta', ('gsm', 'gse', 'title', 'characteristics_ch1', 'source_name_ch1'),
    tqdm((
      dict(
        gsm=gsm,
        gse=row['gse'],
        title=row['title'],
        characteristics_ch1=row['characteristics_ch1'],
        source_name_ch1=row['source_name_ch1'],
      )
      for gsm, row in samps_df_to_ingest.iterrows()
    ),
    total=len(samps_df_to_ingest),
    desc='Inserting GSM info..')
  )

  plpy.execute('refresh materialized view app_public_v2.gene_set_pmid', [])
  plpy.execute('refresh materialized view app_public_v2.gene_set_gse;', [])


def replace_infinity_with_none(obj):
    import math
    if isinstance(obj, float) and math.isinf(obj):
        return None
    elif isinstance(obj, dict):
        return {k: replace_infinity_with_none(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [replace_infinity_with_none(elem) for elem in obj]
    else:
        return obj

def import_enrichr_terms(plpy, species='human'):
  import json
  from tqdm import tqdm

  with open(f'data/enrichr-terms-{species}.json') as f:
    enrichr_terms = json.load(f)

  enrichr_terms_ingested = set([
    r['sig']
    for r in plpy.cursor(
      f'''
        select sig
        from app_public_v2.enrichr_terms
      '''
    )
  ])

  for sig in tqdm(enrichr_terms):
    sql = """
        INSERT INTO app_public_v2.enrichr_terms (sig, organism, sig_terms, enrichr_stats)
        VALUES (%s, %s, %s, %s)
        ON CONFLICT (sig, organism) DO UPDATE SET
        sig_terms = EXCLUDED.sig_terms,
        enrichr_stats = EXCLUDED.enrichr_stats;
        """
    sig_name = list(sig.keys())[0]
    if sig_name in enrichr_terms_ingested:
      continue
    enrichr_stats=json.dumps(replace_infinity_with_none(sig[list(sig.keys())[0]]))

    sig_terms = []
    for lib in sig[sig_name]:
      sig_terms += [items[0] for items in sig[list(sig.keys())[0]][lib] if items[2] < 0.05]
    plpy.execute(sql, (sig_name, species, sig_terms, enrichr_stats))



@click.group()
def cli(): pass
@cli.command()
@click.option('-i', '--input', type=click.Path(exists=True, file_okay=True, path_type=Path), help='GMT file to ingest')
@click.option('--prefix', type=str, default='', help='Prefix to add to terms')
@click.option('--postfix', type=str, default='', help='Postfix to add to terms')
@click.option('--species', type=str, default='human', help='Terms species')
def ingest(input, prefix, postfix, species):
  from plpy import plpy
  print(f'Ingesting {input}...')
  print(f'species: {species}')
  try:
    import_gene_set_library(plpy, input, prefix=prefix, postfix=postfix, species=species)
  except:
    plpy.conn.rollback()
    raise
  else:
    plpy.conn.commit()

@cli.command()
@click.option('--species', type=str, default='human', help='Terms species')
def ingest_gse_info(species):
  from plpy import plpy
  try:
    import_gse_info(plpy, species)
  except:
    plpy.conn.rollback()
    raise
  else:
    plpy.conn.commit()

@cli.command()
@click.option('--species', type=str, default='human', help='Terms species')
def ingest_enrichr_terms(species):
  from plpy import plpy
  try:
    import_enrichr_terms(plpy, species)
  except:
    plpy.conn.rollback()
    raise
  else:
    plpy.conn.commit()

@cli.command()
@click.option('--species', type=str, default='human', help='Terms species')
def ingest_gse_attrs(species):
  from plpy import plpy
  try:
    import_gse_attrs(plpy, species)
  except:
    plpy.conn.rollback()
    raise
  else:
    plpy.conn.commit()

@cli.command()
def ingest_pb_info():
  from plpy import plpy
  try:
    import_pb_info(plpy)
  except:
    plpy.conn.rollback()
    raise
  else:
    plpy.conn.commit()

@cli.command()
def ingest_term_categories():
  from plpy import plpy
  try:
    import_term_categories(plpy)
  except:
    plpy.conn.rollback()
    raise
  else:
    plpy.conn.commit()


@cli.command()
@click.option('--enrich-url', envvar='ENRICH_URL', default='http://127.0.0.1:8000')
@click.option('--species', type=str, default='human')
def update_background(enrich_url, species):
  ''' A background is tied to a complete set of genes across all gene sets of a given species
  but also to a computed index in the enrich API. This function creates a
  new one, and drops the old one after ensuring the index is ready.
  '''
  import requests
  from plpy import plpy
  # record current backgrounds
  current_backgrounds = [row['id'] for row in plpy.cursor(f"select id, species from app_public_v2.background where species = '{species}'")]  # create updated background
  new_background, = plpy.cursor(f'''
    insert into app_public_v2.background (gene_ids, n_gene_ids, species)
    select
      jsonb_object_agg(distinct gsg.gene_id, null) as gene_ids,
      count(distinct gsg.gene_id) as n_gene_ids,
      '{species}' as species
    from app_public_v2.gene_set gs, jsonb_each(gs.gene_ids) gsg(gene_id, nil)
    where gs.species = '{species}'
    returning id;
  ''')
  plpy.conn.commit()
  # trigger index creation for the new background
  assert requests.get(f"{enrich_url}/{new_background['id']}").ok
  # remove old backgrounds
  plpy.execute(
    plpy.prepare('delete from app_public_v2.background where id = any($1::uuid[])', ['text[]']),
    [current_backgrounds]
  )
  plpy.conn.commit()
  # remove index for the old background
  for current_background in current_backgrounds:
    requests.delete(f"{enrich_url}/{current_background}")

if __name__ == '__main__':
  cli()
