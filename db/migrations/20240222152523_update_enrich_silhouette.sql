-- migrate:up

create type app_public_v2.enriched_term_result as (
  term varchar,
  count int,
  odds_ratio double precision,
  pvalue double precision,
  adj_pvalue double precision,
  not_term_count int,
  total_term_count int,
  total_not_term_count int
);

drop type app_public_v2.paginated_enrich_result cascade;

create type app_public_v2.paginated_enrich_result as (
  nodes app_public_v2.enrich_result[],
  total_count int,
  enriched_terms app_public_v2.enriched_term_result[]
);

create materialized view app_public_v2.gse_attrs_terms_count as
select count(elts) as total_terms
from (
  select unnest(gse_attrs) as elts
  from app_public_v2.gse_info
) subquery;

create materialized view app_public_v2.gse_attrs_term_counts as
select elts as term, count(elts) as term_count
from (
  select unnest(gse_attrs) as elts
  from app_public_v2.gse_info
) subquery
group by elts;

create or replace function app_private_v2.enrich_functional_terms(
    concat_terms varchar[], total_count bigint) returns app_public_v2.enriched_term_result[] as $$
    from scipy.stats import fisher_exact
    from statsmodels.stats.multitest import multipletests
    from collections import Counter

    term_counter = Counter(concat_terms)
    term_counts = list(term_counter.items())
    total_enrich_term_count = sum(term_counter.values())

    results = []
    p_values = []
    for term, count in term_counts:
        term_count_result = plpy.execute(f"SELECT term_count FROM app_public_v2.gse_attrs_term_counts WHERE term = '{term}'")
        total_term_count = term_count_result[0]['term_count'] if term_count_result else 0
        contingency_table = [[count, total_enrich_term_count - count], [total_term_count, total_count - total_term_count]]
        odds_ratio, p_value = fisher_exact(contingency_table)
        p_values.append(p_value)
        results.append({
            'term': term,
            'count': count,
            'odds_ratio': odds_ratio,
            'pvalue': p_value,
            'not_term_count': total_enrich_term_count - count,
            'total_term_count': total_term_count,
            'total_not_term_count': total_count - total_term_count
        })
    try:
        adjusted_p_values = multipletests(p_values, method='fdr_bh')[1]
        for i in range(len(results)):
            results[i]['adj_pvalue'] = adjusted_p_values[i]
    except Exception as e:
        for i in range(len(results)):
            results[i]['adj_pvalue'] = 1
    results = list(sorted(filter(lambda r: r['adj_pvalue'] < .05, results), key=lambda r: r['adj_pvalue']))
    return results
$$ language plpython3u immutable parallel safe;


create or replace function app_private_v2.enriched_functional_terms(
    enriched_terms varchar[]) returns app_public_v2.enriched_term_result[] as $$ 
    with gse_ids as (
        -- Split enriched_terms vector to get GSE ids
        select regexp_replace(t, '\mGSE([^-]+)\M.*', 'GSE\1') as id
        from unnest(enriched_terms) as t
    ), gse_terms as (
    -- Get list of functional terms called gse_attrs from app_public_v2.gse_info table
        select array_agg(terms) as terms
        from (
            select unnest(gse_attrs) as terms
            from app_public_v2.gse_info
            where gse in (select id from gse_ids)
        ) subquery
    )
    -- Call the enrich_functional_terms function
    select * from app_private_v2.enrich_functional_terms((select terms from gse_terms), (select total_terms from app_public_v2.gse_attrs_terms_count));
$$ language sql immutable parallel safe security definer;
grant execute on function app_private_v2.enriched_functional_terms to guest, authenticated;

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
    filter_score=filter_score_le
  )
  if filter_term: params['filter_term'] = filter_term
  if offset: params['offset'] = offset
  if first: params['limit'] = first
  req = requests.post(
    f"http://rummageo-enrich:8000/{background['id']}",
    params=params,
    json=gene_ids,
  )
  total_count = req.headers.get('Content-Range').partition('/')[-1]
  req_json = req.json()
  enriched_terms = req_json.pop('terms')
  enriched_terms = ','.join(f"'{item}'" for item in enriched_terms)
  enriched_terms_result = plpy.execute(f"SELECT app_private_v2.enriched_functional_terms(ARRAY[{enriched_terms}::varchar])")[0]['enriched_functional_terms']
  return dict(nodes=req_json['results'], total_count=total_count, enriched_terms=enriched_terms_result)
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
drop function app_public_v2.background_enrich;
drop function app_private_v2.indexed_enrich;
drop type app_public_v2.enriched_term_result cascade;
drop materialized view app_public_v2.gse_attrs_terms_count;
drop materialized view gse_attrs_term_counts;

