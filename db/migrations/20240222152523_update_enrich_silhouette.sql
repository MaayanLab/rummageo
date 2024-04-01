-- migrate:up

create type app_public_v2.enriched_term_result as (
  term varchar,
  statistic double precision,
  pvalue double precision,
  adj_pvalue double precision,
  count int,
  total_term_count int
);


drop type app_public_v2.paginated_enrich_result cascade;

create type app_public_v2.paginated_enrich_result as (
  nodes app_public_v2.enrich_result[],
  total_count int,
  enriched_terms varchar[]
);


create materialized view app_public_v2.llm_attrs_term_counts_human as
select elts as term, count(elts) as term_count
from (
  select unnest(llm_attrs) as elts
  from app_public_v2.gse_terms
  where species = 'human'
) subquery
group by elts;

create materialized view app_public_v2.pubmed_attrs_term_counts_human as
select elts as term, count(elts) as term_count
from (
  select unnest(pubmed_attrs) as elts
  from app_public_v2.gse_terms
  where species = 'human'
) subquery
group by elts;

create materialized view app_public_v2.mesh_attrs_term_counts_human as
select elts as term, count(elts) as term_count
from (
  select unnest(mesh_attrs) as elts
  from app_public_v2.gse_terms
  where species = 'human'
) subquery
group by elts;

create materialized view app_public_v2.llm_attrs_term_counts_mouse as
select elts as term, count(elts) as term_count
from (
  select unnest(llm_attrs) as elts
  from app_public_v2.gse_terms
  where species = 'mouse'
) subquery
group by elts;

create materialized view app_public_v2.pubmed_attrs_term_counts_mouse as
select elts as term, count(elts) as term_count
from (
  select unnest(pubmed_attrs) as elts
  from app_public_v2.gse_terms
  where species = 'mouse'
) subquery
group by elts;

create materialized view app_public_v2.mesh_attrs_term_counts_mouse as
select elts as term, count(elts) as term_count
from (
  select unnest(mesh_attrs) as elts
  from app_public_v2.gse_terms
  where species = 'mouse'
) subquery
group by elts;

create or replace function app_private_v2.enrich_functional_terms(
    concat_terms varchar[], source_type varchar, species varchar) returns app_public_v2.enriched_term_result[] as $$
    from scipy.stats import kstest
    from statsmodels.stats.multitest import multipletests
    from collections import Counter
    import math
    rank_dict = {}

    if concat_terms is None:
        return []

    for i, term in enumerate(concat_terms):
        if term not in rank_dict:
            rank_dict[term] = []
        rank_dict[term].append(i)
      

    results = []
    p_values = []
    for term in rank_dict:
      if len(rank_dict[term]) < 5:
          continue
      try:
        statistic, pvalue = kstest(rank_dict[term], 'uniform', args=(0, 9999), alternative='greater')
        p_values.append(pvalue)
        escaped_term = term.replace("'", "''")
        term_count_result = plpy.execute(f"SELECT term_count FROM app_public_v2.{source_type}_term_counts_{species} WHERE term = '{escaped_term}'")
        total_term_count = term_count_result[0]['term_count'] if term_count_result else 0
        results.append({
            'term': term,
            'statistic': statistic,
            'pvalue': pvalue,
            'count': len(rank_dict[term]),
            'total_term_count': total_term_count
        })
      except Exception as e:
        print(e)
        continue
    try:
        adjusted_p_values = multipletests(p_values, method='fdr_bh')[1]
        for i in range(len(results)):
            results[i]['adj_pvalue'] = adjusted_p_values[i]
    except Exception as e:
        print(e)
        for i in range(len(results)):
            results[i]['adj_pvalue'] = 1
    sig_results = list(filter(lambda r: r['adj_pvalue'] < .05, results))

    if len(sig_results) > 50:
        results = list(sorted(sig_results, key=lambda r: r['adj_pvalue']))[:50]
    else:
        results = list(sorted(results, key=lambda r: r['pvalue']))[:100]
    return results
$$ language plpython3u immutable parallel unsafe;


create or replace function app_public_v2.enriched_functional_terms(
  enriched_terms varchar[], source_type varchar, organism varchar) returns app_public_v2.enriched_term_result[] as $$ 
  with gse_ids as (
    -- Split enriched_terms vector to get GSE ids
    select regexp_replace(t, '\mGSE([^-]+)\M.*', 'GSE\1') as id
    from unnest(enriched_terms) as t
  ), gse_terms as (
    -- Get list of functional terms based on the source_type
    select array_agg(terms) as terms
    from (
      select unnest(
        case 
          when source_type = 'llm_attrs' then llm_attrs
          when source_type = 'pubmed_attrs' then pubmed_attrs
          when source_type = 'mesh_attrs' then mesh_attrs
        end
      ) as terms
      from app_public_v2.gse_terms
      where gse in (select id from gse_ids) and species = organism
    ) subquery
    )
    -- Call the enrich_functional_terms function
    select * from app_private_v2.enrich_functional_terms(
      (select terms from gse_terms), 
      source_type, organism)
$$ language sql immutable parallel unsafe security definer;
grant execute on function app_public_v2.enriched_functional_terms to guest, authenticated;

create or replace function app_private_v2.indexed_enrich(
  background app_public_v2.background,
  gene_ids uuid[],
  filter_term varchar default null,
  overlap_ge int default 1,
  pvalue_le double precision default 0.05,
  adj_pvalue_le double precision default 0.05,
  "offset" int default null,
  "first" int default null,
  filter_score_le double precision default -1
) returns app_public_v2.paginated_enrich_result as $$
  import requests
  import json
  
  params = dict(
    overlap_ge=overlap_ge,
    pvalue_le=pvalue_le,
    adj_pvalue_le=adj_pvalue_le,
    score_filter=filter_score_le
  )
  if filter_term: params['filter_term'] = filter_term
  if offset: params['offset'] = offset
  if first: params['limit'] = first
  req = requests.post(
    f"http://rummageo-enrich:8000/{background['id']}",
    params=params,
    json=gene_ids,
  )
  print(filter_score_le)
  total_count = req.headers.get('Content-Range').partition('/')[-1]
  req_json = req.json()
  enriched_terms = req_json.pop('terms')
  gses = set()
  enriched_terms_top_gses = []
  for term in enriched_terms:
    gse = term.split('-')[0]
    if gse not in gses:
      gses.add(gse)
      enriched_terms_top_gses.append(gse)
    if len(enriched_terms_top_gses) >= 10000:
      break
  return dict(nodes=req_json['results'], total_count=total_count, enriched_terms=enriched_terms_top_gses)
$$ language plpython3u immutable parallel safe;

create or replace function app_public_v2.background_enrich(
  background app_public_v2.background,
  genes varchar[],
  filter_term varchar default null,
  overlap_ge int default 1,
  pvalue_le double precision default 0.05,
  adj_pvalue_le double precision default 0.05,
  "offset" int default null,
  "first" int default null,
  filter_score_le double precision default -1
) returns app_public_v2.paginated_enrich_result
as $$
  select r.*
  from app_private_v2.indexed_enrich(
    background_enrich.background,
    (select array_agg(gene_id) from app_public_v2.gene_map(genes) gm),
    background_enrich.filter_term,
    background_enrich.overlap_ge,
    background_enrich.pvalue_le,
    background_enrich.adj_pvalue_le,
    background_enrich."offset",
    background_enrich."first",
    background_enrich.filter_score_le
  ) r;
$$ language sql immutable parallel safe security definer;
grant execute on function app_public_v2.background_enrich to guest, authenticated;

alter table app_public_v2.gse_info drop column gse_attrs cascade;
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
WHERE  gi.gse = gt.gse and gi.species = gt.species;

create materialized view app_public_v2.gene_set_pmid as
select
    gs.id,
    gs.term,
    gse_info.id as gse_id,
    gse_info.gse,
    gse_info.pmid,
    gse_info.title,
    gse_info.sample_groups,
    gse_info.platform,
    gse_info.published_date,
    gse_info.gse_attrs,
    gse_info.silhouette_score
from
    app_public_v2.gene_set gs
join
    app_public_v2.gse_info gse_info ON regexp_replace(gs.term, '\mGSE([^-]+)\M.*', 'GSE\1') = gse_info.gse;


comment on materialized view app_public_v2.gene_set_pmid is E'@foreignKey (id) references app_public_v2.gene_set (id)';

create index gene_set_pmid_id_idx on app_public_v2.gene_set_pmid (id);
create index gene_set_pmid_pmid_idx on app_public_v2.gene_set_pmid (pmid);

grant select on app_public_v2.gene_set_pmid to guest;
grant all privileges on app_public_v2.gene_set_pmid to authenticated;


create function app_public_v2.get_pb_info_by_ids(pmids varchar[])
returns setof app_public_v2.gene_set_pmid as
$$
  select *
  from app_public_v2.gene_set_pmid
  where pmid = ANY(pmids)
$$ language sql immutable strict parallel safe;
grant execute on function app_public_v2.get_pb_info_by_ids to guest, authenticated;

create or replace function app_public_v2.gene_set_term_search(terms varchar[]) returns setof app_public_v2.gene_set_pmid
as $$
  select distinct gs.*
  from app_public_v2.gene_set_pmid gs
  inner join unnest(terms) ut(term) on gs.title ilike ('%' || ut.term || '%');
$$ language sql immutable strict parallel safe;
grant execute on function app_public_v2.gene_set_term_search to guest, authenticated;

-- migrate:down
drop function app_public_v2.background_enrich;
drop function app_private_v2.indexed_enrich;
drop type app_public_v2.enriched_term_result cascade;
drop materialized view app_public_v2.gene_set_pmid cascade;

drop materialized view app_public_v2.llm_attrs_term_counts_human;
drop materialized view app_public_v2.pubmed_attrs_term_counts_human;
drop materialized view app_public_v2.mesh_attrs_term_counts_human;

drop materialized view app_public_v2.llm_attrs_term_counts_mouse;
drop materialized view app_public_v2.pubmed_attrs_term_counts_mouse;
drop materialized view app_public_v2.mesh_attrs_term_counts_mouse;

