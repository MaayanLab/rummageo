import re
import io
import os
import csv
import sys
import uuid
import json
import queue
import shutil
import tarfile
import tempfile
import traceback
import contextlib
import subprocess
import multiprocessing as mp
from multiprocessing.pool import ThreadPool
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path, PurePosixPath

import tqdm
import numpy as np
import pandas as pd

from docx import Document
from maayanlab_bioinformatics.harmonization.ncbi_genes import ncbi_genes_lookup

# import tabula
# java = shutil.which('java')
# assert java, 'Missing java, necessary for tabula-py'

soffice = shutil.which('soffice', path=':'.join(filter(None, [os.environ.get('PATH'), '/Applications/LibreOffice.app/Contents/MacOS/'])))
assert soffice, 'Missing `soffice` binary for converting doc to docx'

class _DevNull:
  ''' File handle that does nothing
  '''
  def write(self, *args, **kwargs): pass
  def flush(self, *args, **kwargs): pass

def _run_with_timeout(send, fn, *args):
  try:
    send.put((None, fn(*args)))
  except Exception as e:
    send.put((e, None))

def run_with_timeout(fn, *args, timeout: int = 60):
  mp_spawn = mp.get_context('spawn')
  recv = mp_spawn.Queue()
  proc = mp_spawn.Process(target=_run_with_timeout, args=(recv, fn, *args))
  proc.start()
  try:
    err, res = recv.get(timeout=timeout)
  except queue.Empty:
    raise TimeoutError()
  else:
    if err is not None:
      raise err
    else:
      return res
  finally:
    proc.join(1)
    if proc.exitcode is None:
      proc.terminate()
      proc.join(1)
      if proc.exitcode is None:
        proc.kill()
        proc.join(1)

ext_handlers = {}
def register_ext_handler(*exts):
  ''' We create a dictionary with functions capable of extracting tables for each extension type.
  Each function is a generator of (name, pandas data frame) tuples
  '''
  def decorator(func):
    for ext in exts:
      ext_handlers[ext] = func
    return func
  return decorator

def _read_docx_tab(tab):
  '''  This converts from a docx table object into a pandas dataframe
  '''
  vf = io.StringIO()
  writer = csv.writer(vf)
  for row in tab.rows:
    writer.writerow(cell.text for cell in row.cells)
  vf.seek(0)
  return pd.read_csv(vf)

def read_docx_tables(f):
  ''' This reads tables out of a docx file
  '''
  with contextlib.redirect_stderr(_DevNull()):
    doc = Document(f)
  for i, tab in enumerate(doc.tables):
    yield dict(label=str(i), df=_read_docx_tab(tab))

@register_ext_handler('.docx')
def prepare_docx(fr):
  ''' This calls read_docx_tables, first copying the reader into a ByteIO since
  tarfile reader doesn't support seeks.
  '''
  fh = io.BytesIO()
  shutil.copyfileobj(fr, fh)
  fh.seek(0)
  yield from read_docx_tables(fh)

@register_ext_handler('.doc')
def read_doc_as_docx(fr):
  ''' For doc support, convert .doc to .docx in a temporary directory and call read_docx_tables
  '''
  with tempfile.TemporaryDirectory() as tmpdir:
    tmpdir = Path(tmpdir)
    with (tmpdir / 'table.doc').open('wb') as fw:
      shutil.copyfileobj(fr, fw)
    subprocess.call(
      [soffice, '--headless', '--convert-to', 'docx', str(tmpdir/'table.doc')],
      cwd=tmpdir,
      stderr=subprocess.DEVNULL,
      stdout=subprocess.DEVNULL,
      timeout=60,
    )
    yield from read_docx_tables(tmpdir/'table.docx')

@register_ext_handler('.xls', '.xlsb', '.xlsm','.odf','.ods','.odt')
def read_excel_tables(f, engine=None):
  ''' Use pandas read_excel function for these files, return all tables from all sheets
  '''
  for sheet, df in pd.read_excel(f, sheet_name=None, engine=engine).items():
    yield dict(label=sheet, df=df)

@register_ext_handler('.xlsx')
def read_xlsx_tables(f):
  yield from read_excel_tables(f, engine='openpyxl')

@register_ext_handler('.csv')
def read_csv_tables(f):
  yield dict(df=pd.read_csv(f))

@register_ext_handler('.tsv')
def read_tsv_tables(f):
  yield dict(df=pd.read_csv(f, sep='\t'))

@register_ext_handler('.txt')
def read_txt_tables(f):
  ''' Try to read txt as a table using pandas infer functionality
  '''
  yield dict(df=pd.read_csv(f, sep=None, engine='python'))

def _read_xml_text(node):
  ''' Read the text from an xml node (or the text from all it's children)
  '''
  return node.text or ''.join(filter(None, (
    el.text
    for el in node.findall('.//')
  )))

def _read_xml_text_with_exclusion(node, exclude={'table-wrap', 'fig'}):
  ''' Read the text from an xml node (or the text from all it's children)
  '''
  Q = [node]
  while Q:
    node = Q.pop(0)
    if node.text:
      yield node
    else:
      Q += [
        el
        for el in node
        if el.tag not in exclude
      ]

def _read_xml_table(tbl):
  ''' This reads a xml table as a pandas dataframe
  '''
  columns = [[_read_xml_text(td).replace('\t', '  ').replace('\n', ' ').strip() for td in tr.findall('./')] for tr in tbl.find('./thead').findall('./tr')]
  values = [[_read_xml_text(td).replace('\t', '  ').replace('\n', ' ').strip() for td in tr.findall('./')] for tr in tbl.find('./tbody').findall('./tr')]
  n_cols = max(map(len, columns + values))
  if n_cols > len(columns[0]):
    columns[0] += ['']*(n_cols - len(columns[0]))
  df = pd.read_csv(
    io.StringIO('\n'.join('\t'.join(el for el in row) for row in (columns + values))),
    sep='\t',
    on_bad_lines='warn',
  )
  return df

def _read_xml_tables(root):
  ''' Tables are embedded in the xml files, they can be parsed 
  '''
  for i, tblWrap in enumerate(root.findall('.//table-wrap'), start=1):
    # find the table label
    label = tblWrap.find('./label')
    if label:
      label = _read_xml_text(label)

    # find the table caption
    caption = tblWrap.find('./caption')
    if caption:
      caption = _read_xml_text(caption)

    # find any mention of the table in the article
    ref = tblWrap.attrib.get('id')
    if ref:
      href = f"#{ref}"
      mentions = '\n'.join(
        ''.join(_read_xml_text_with_exclusion(mention))
        for mention in root.findall(f'.//xref[rid="{ref}"]/..')
      )
    else:
      href = None
      mentions = None

    # get the table itself
    tbl = tblWrap.find('./table')
    if tbl and tbl.find('./thead') and tbl.find('./tbody'):
      try:
        tbl = _read_xml_table(tbl)
      except KeyboardInterrupt:
        raise
      except:
        traceback.print_exc()
        continue
      #
      yield dict(type='table', href=href, df=tbl)
      yield dict(type='context', href=href, label=label, caption=caption, mentions=mentions)

def _read_xml_supplemental_context(root):
  ''' Tables are embedded in the xml files, they can be parsed 
  '''
  for supplementary_material in root.findAll('.//supplementary-material'):
    # find the supplemental material href
    media = root.find('./media')
    if not media: continue
    href = media.attrib['xlink:href']

    # find the supplemental material caption
    caption = media.find('./caption')
    if caption:
      caption = _read_xml_text(caption)

    # find any mention of the table in the article
    ref = supplementary_material.attrib.get('id')
    if ref:
      mentions = '\n'.join(
        ''.join(_read_xml_text_with_exclusion(mention))
        for mention in root.findall(f'.//xref[rid="{ref}"]/..')
      )
    else:
      mentions = None

    yield dict(
      type='context',
      href=href,
      caption=caption,
      mentions=mentions,
    )

@register_ext_handler('.nxml', '.xml')
def read_xml(f):
  ''' Tables are embedded in the xml files, they can be parsed 
  '''
  parsed = ET.parse(f)
  root = parsed.getroot()
  yield from _read_xml_tables(root)
  yield from _read_xml_supplemental_context(root)

# it's unclear whether this really gives us any new information, and it is extremely expensive to compute
#  so for now I'm leaving it out.
# @register_ext_handler('.pdf')
# def read_pdf_tables(f):
#   ''' pdf tables read by tabula library
#   '''
#   results = tabula.read_pdf(f, pages='all', multiple_tables=True, silent=True)
#   if type(results) == list:
#     for i, df in enumerate(results):
#       yield f"{i}", df
#   elif type(results) == dict:
#     for key, df in results.items():
#       yield key, df
#   else:
#     raise NotImplementedError()

lookup = None
def gene_lookup(value):
  ''' Don't allow pure numbers or spaces--numbers can typically match entrez ids
  '''
  if type(value) != str: return None
  if re.search(r'\s', value): return None
  if re.match(r'\d+(\.\d+)?', value): return None
  global lookup
  if lookup is None:
    lookup = ncbi_genes_lookup(filters=lambda ncbi: ncbi)
  return lookup(value)

def extract_gene_set_columns(df):
  ''' Given a pandas dataframe, find columns containing mostly mappable genes
  '''
  for column in df.columns:
    if df[column].dtype != np.dtype('O'): continue
    unique_genes = pd.Series(df[column].dropna().unique())
    if unique_genes.shape[0] >= 5:
      unique_genes_mapped = unique_genes.apply(gene_lookup).dropna()
      ratio = unique_genes_mapped.shape[0] / unique_genes.shape[0]
      if ratio > 0.5:
        yield dict(
          column=column,
          raw_genes=unique_genes.apply(lambda gene: re.sub(r'\s+', ' ', gene) if type(gene) == str else gene).tolist(),
          mapped_genes=unique_genes_mapped.tolist(),
        ) 

def slugify(s):
  ''' Replace non-characters/numbers with _
  '''
  return re.sub(r'[^\w\d-]+', '_', s).strip('_')

def extract_tables_from_oa_package(oa_package):
  ''' Load all the tables from an oa_package archive
  '''
  with tarfile.open(oa_package) as tar:
    for member in tar.getmembers():
      if member.isfile():
        member_name_path = PurePosixPath(member.name)
        handler = ext_handlers.get(member_name_path.suffix.lower())
        if handler:
          try:
            for info in enumerate(handler(tar.extractfile(member))):
              yield dict(member_name_path=str(member_name_path), **info)
          except KeyboardInterrupt:
            raise
          except:
            traceback.print_exc()

def extract_oa_package(oa_package):
  for info in extract_tables_from_oa_package(oa_package):
    table_id = uuid.uuid4()
    df = info.pop('df')
    if df:
      yield dict(type='table', id=table_id, **info)
      for column in extract_gene_set_columns(df):
        yield dict(type='column', table_id=table_id, **column)
    else:
      yield dict(type='paper', id=oa_package, **info)

def fetch_oa_file_list(data_dir = Path()):
  ''' Fetch the PMCID, PMID, oa_file listing; we sort it newest first.
  ['File'] has the oa_package which is a relative path to a tar.gz archive containing
   the paper and all figures.
  '''
  oa_file_list = data_dir / 'oa_file_list.csv'
  if not oa_file_list.exists():
    df = pd.read_csv('https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_file_list.csv')
    ts_col = df.columns[-3]
    df[ts_col] = pd.to_datetime(df[ts_col])
    df.sort_values(ts_col, ascending=False, inplace=True)
    df.to_csv(oa_file_list, index=None)
  else:
    df = pd.read_csv(oa_file_list)
  return df

def find_pmc_ids(term):
  ''' Given a term, return all PMC ids matching that term
  '''
  import os, itertools
  from Bio import Entrez
  Entrez.email = os.environ['EMAIL']
  batch = 1000000
  for i in itertools.count():
    try:
      handle = Entrez.esearch(db="pmc", term=term, api_key=os.environ['API_KEY'], retstart=i*batch, retmax=batch)
      records = Entrez.read(handle)
      if not records['IdList']:
        break
      for id in records['IdList']:
        yield f"PMC{id}"
    except KeyboardInterrupt:
      raise
    except:
      import traceback
      traceback.print_exc()
      break

def filter_oa_file_list_by(oa_file_list, pmc_ids):
  ''' Filter oa_file_list by PMC IDs
  '''
  return oa_file_list[oa_file_list['Accession ID'].isin(list(pmc_ids))]

def fetch_extract_oa_package(oa_package):
  ''' Given the oa_package name from the oa_file_list, we'll download it temporarily and then extract tables/columns out of it
  '''
  with tempfile.NamedTemporaryFile(suffix=''.join(PurePosixPath(oa_package).suffixes)) as tmp:
    with urllib.request.urlopen(f"https://ftp.ncbi.nlm.nih.gov/pub/pmc/{oa_package}") as fr:
      shutil.copyfileobj(fr, tmp)
    tmp.flush()
    return list(extract_oa_package(tmp.name))

def task(record):
  try:
    return record, None, run_with_timeout(fetch_extract_oa_package, record['File'], timeout=60*5)
  except KeyboardInterrupt:
    raise
  except:
    return record, traceback.format_exc(), None

def main(data_dir = Path(), oa_file_list = None, progress = 'done.txt', progress_output = 'done.new.txt', output = 'output.jsonl'):
  '''
  Work through oa_file_list (see: fetch_oa_file_list)
    -- you can filter it and provide it to this function
  Track progress by storing oa_packages already processed in done.txt
  Write all results to output.jsonl
  '''
  data_dir.mkdir(parents=True, exist_ok=True)
  done_file = data_dir / progress
  new_done_file = data_dir / progress_output
  output_file = data_dir / output

  # find out what we've already processed
  if done_file.exists():
    with done_file.open('r') as fr:
      done = set(filter(None, map(str.strip, fr)))
  else:
    done = set()

  # find out what there remains to process
  if oa_file_list is None:
    oa_file_list = fetch_oa_file_list(data_dir)

  oa_file_list_size = oa_file_list.shape[0]
  oa_file_list = oa_file_list[~oa_file_list['File'].isin(list(done))]

  # fetch and extract oa_packages using a process pool
  #  append extracted records as they are ready into one jsonl file
  with new_done_file.open('a') as done_file_fh:
    with output_file.open('a') as output_fh:
      with ThreadPool() as pool:
        for record, error, res in tqdm.tqdm(
          pool.imap_unordered(
            task,
            (row for _, row in oa_file_list.iterrows())
          ),
          initial=oa_file_list_size - oa_file_list.shape[0],
          total=oa_file_list_size
        ):
          if error is None:
            for record in res:
              print(json.dumps(record), file=output_fh)
          else:
            print(json.dumps(dict(type='error', error=error)), file=sys.stderr)
          print(record['File'], file=done_file_fh)
          output_fh.flush()
          done_file_fh.flush()

if __name__ == '__main__':
  import os
  from dotenv import load_dotenv; load_dotenv()
  data_dir = Path(os.environ.get('PTH', 'data'))
  main(data_dir)
