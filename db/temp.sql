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



create or replace function app_private_v2.enrich_functional_terms(
    concat_terms varchar[], source_type varchar, species varchar) returns app_public_v2.enriched_term_result[] as $$
    from scipy.stats import fisher_exact
    from statsmodels.stats.multitest import multipletests
    from collections import Counter
    import math

    total_count = plpy.execute(f"SELECT {source_type}_total FROM app_public_v2.gse_attrs_terms_count_{species}")[0][f'{source_type}_total']
    total_documents = plpy.execute(f"SELECT count(gse) as gse_count FROM app_public_v2.gse_terms where species = '{species}'")[0][f'gse_count']

    term_counter = Counter(concat_terms)
    term_counts = list(term_counter.items())
    total_enrich_term_count = sum(term_counter.values())

    results = []
    p_values = []
    for term, count in term_counts:
      try:
        escaped_term = term.replace("'", "''")
        term_count_result = plpy.execute(f"SELECT term_count FROM app_public_v2.{source_type}_term_counts_{species} WHERE term = '{escaped_term}'")
        total_term_count = term_count_result[0]['term_count'] if term_count_result else 0
        contingency_table = [[count, total_enrich_term_count - count], [total_term_count, total_count - total_term_count]]
        odds_ratio, p_value = fisher_exact(contingency_table)
        #p_value = (1 / math.log((total_documents/total_count) + 1)) * p_value
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