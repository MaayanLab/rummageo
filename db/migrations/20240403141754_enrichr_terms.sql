-- migrate:up


create type app_public_v2.enrichr_result as (
  term varchar,
  pvalue double precision,
  adj_pvalue double precision,
  count int
);


/* UPDATE app_public_v2.gene_set */
/* SET term = replace(term, '.tsv', ''); */

grant select on app_public_v2.enrichr_terms to guest;
grant all privileges on app_public_v2.enrichr_terms to authenticated;

create or replace function app_private_v2.sig_enrichr_terms(
    concat_terms jsonb, species varchar) returns app_public_v2.enrichr_result[] as $$
    from scipy.stats import kstest
    from statsmodels.stats.multitest import multipletests
    from collections import Counter
    import json

    import math
    rank_dict = {}

    if concat_terms is None:
        return []

    concat_terms_json = json.loads(concat_terms)
    
    for rank in concat_terms_json:
        for term in concat_terms_json[rank]:
            if term not in rank_dict:
                rank_dict[term] = []
            rank_dict[term].append(int(rank))

    results = []
    p_values = []
    for term in rank_dict:
      if len(rank_dict[term]) < 5:
          continue
      try:
        statistic, pvalue = kstest(rank_dict[term], 'uniform', args=(0, 500), alternative='greater')
        p_values.append(pvalue)
        escaped_term = term.replace("'", "''")
        results.append({
            'term': term,
            'statistic': statistic,
            'pvalue': pvalue,
            'count': len(rank_dict[term]),
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
$$ language plpython3u immutable parallel safe;


create or replace function app_public_v2.enriched_enrichr_terms(
    enriched_sigs varchar[], organism varchar) RETURNS app_public_v2.enrichr_result[] AS $$ 
with enriched_sigs_with_row_number as (
    select enriched_sigs[i] as sig, i as rn
    from generate_subscripts(enriched_sigs, 1) as i
), enrichr_attrs AS (
    select jsonb_object_agg(es.rn::text, sig_terms) as attrs
    from app_public_v2.enrichr_terms et
    join enriched_sigs_with_row_number es ON et.sig = es.sig and et.organism = organism
    where cardinality(sig_terms) > 0
)
-- Call the enrich_functional_terms function and order by row number
select * from app_private_v2.sig_enrichr_terms((SELECT attrs FROM enrichr_attrs), organism)
$$ language sql immutable parallel safe security definer;
grant execute on function app_public_v2.enriched_enrichr_terms to guest, authenticated;



-- migrate:down

drop function app_public_v2.enriched_enrichr_terms(enriched_sigs varchar[], organism varchar);
drop function app_private_v2.sig_enrichr_terms(concat_terms jsonb, species varchar);
drop type app_public_v2.enrichr_result;

