-- migrate:up

drop type app_public_v2.paginated_enrich_result cascade;

create type app_public_v2.paginated_enrich_result as (
  nodes app_public_v2.enrich_result[],
  total_count int,
  enriched_terms varchar[],
  top_enriched_sigs varchar[]
);

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


-- migrate:down

