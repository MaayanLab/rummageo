-- migrate:up
create or replace function app_private_v2.indexed_enrich(
  background app_public_v2.background,
  gene_ids uuid[],
  filter_term varchar default null,
  overlap_ge int default 1,
  pvalue_le double precision default 0.05,
  adj_pvalue_le double precision default 0.05,
  "offset" int default null,
  "first" int default null,
  filter_score_le double precision default -1,
  sort_by varchar default null,
  sort_by_dir varchar default null
) returns app_public_v2.paginated_enrich_result as $$
  import os
  import requests
  import json
  
  params = dict(
    overlap_ge=overlap_ge,
    pvalue_le=pvalue_le,
    adj_pvalue_le=adj_pvalue_le,
    score_filter=filter_score_le,
    sort_by=sort_by,
    sort_by_dir=sort_by_dir
  )
  if filter_term: params['filter_term'] = filter_term
  if offset: params['offset'] = offset
  if first: params['limit'] = first
  req = requests.post(
    f"{os.environ.get('ENRICH_URL', 'http://rummageo-enrich:8000')}/{background['id']}",
    params=params,
    json=gene_ids,
  )
  print(filter_score_le)
  total_count = req.headers.get('Content-Range').partition('/')[-1]
  req_json = req.json()
  enriched_terms = req_json.pop('terms')
  gses = set()
  enriched_terms_top_gses = []
  enriched_terms_top_sigs = []
  for i, term in enumerate(enriched_terms):
    if i < 500:
        enriched_terms_top_sigs.append(term.replace('.tsv', ''))
    gse = term.split('-')[0]
    if gse not in gses:
      gses.add(gse)
      enriched_terms_top_gses.append(term)
    if len(enriched_terms_top_gses) >= 5000:
      break
  return dict(nodes=req_json['results'], total_count=total_count, enriched_terms=enriched_terms_top_gses, top_enriched_sigs=enriched_terms_top_sigs)
$$ language plpython3u immutable parallel safe;

-- migrate:down
create or replace function app_private_v2.indexed_enrich(
  background app_public_v2.background,
  gene_ids uuid[],
  filter_term varchar default null,
  overlap_ge int default 1,
  pvalue_le double precision default 0.05,
  adj_pvalue_le double precision default 0.05,
  "offset" int default null,
  "first" int default null,
  filter_score_le double precision default -1,
  sort_by varchar default null,
  sort_by_dir varchar default null
) returns app_public_v2.paginated_enrich_result as $$
  import requests
  import json
  
  params = dict(
    overlap_ge=overlap_ge,
    pvalue_le=pvalue_le,
    adj_pvalue_le=adj_pvalue_le,
    score_filter=filter_score_le,
    sort_by=sort_by,
    sort_by_dir=sort_by_dir
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
  enriched_terms_top_sigs = []
  for i, term in enumerate(enriched_terms):
    if i < 500:
        enriched_terms_top_sigs.append(term.replace('.tsv', ''))
    gse = term.split('-')[0]
    if gse not in gses:
      gses.add(gse)
      enriched_terms_top_gses.append(term)
    if len(enriched_terms_top_gses) >= 5000:
      break
  return dict(nodes=req_json['results'], total_count=total_count, enriched_terms=enriched_terms_top_gses, top_enriched_sigs=enriched_terms_top_sigs)
$$ language plpython3u immutable parallel safe;
