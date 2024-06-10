SET statement_timeout = 0;
SET lock_timeout = 0;
SET idle_in_transaction_session_timeout = 0;
SET client_encoding = 'UTF8';
SET standard_conforming_strings = on;
SELECT pg_catalog.set_config('search_path', '', false);
SET check_function_bodies = false;
SET xmloption = content;
SET client_min_messages = warning;
SET row_security = off;

--
-- Name: app_private_v2; Type: SCHEMA; Schema: -; Owner: -
--

CREATE SCHEMA app_private_v2;


--
-- Name: app_public_v2; Type: SCHEMA; Schema: -; Owner: -
--

CREATE SCHEMA app_public_v2;


--
-- Name: plpython3u; Type: EXTENSION; Schema: -; Owner: -
--

CREATE EXTENSION IF NOT EXISTS plpython3u WITH SCHEMA pg_catalog;


--
-- Name: EXTENSION plpython3u; Type: COMMENT; Schema: -; Owner: -
--

COMMENT ON EXTENSION plpython3u IS 'PL/Python3U untrusted procedural language';


--
-- Name: pg_trgm; Type: EXTENSION; Schema: -; Owner: -
--

CREATE EXTENSION IF NOT EXISTS pg_trgm WITH SCHEMA public;


--
-- Name: EXTENSION pg_trgm; Type: COMMENT; Schema: -; Owner: -
--

COMMENT ON EXTENSION pg_trgm IS 'text similarity measurement and index searching based on trigrams';


--
-- Name: uuid-ossp; Type: EXTENSION; Schema: -; Owner: -
--

CREATE EXTENSION IF NOT EXISTS "uuid-ossp" WITH SCHEMA public;


--
-- Name: EXTENSION "uuid-ossp"; Type: COMMENT; Schema: -; Owner: -
--

COMMENT ON EXTENSION "uuid-ossp" IS 'generate universally unique identifiers (UUIDs)';


--
-- Name: enrich_result; Type: TYPE; Schema: app_public_v2; Owner: -
--

CREATE TYPE app_public_v2.enrich_result AS (
	gene_set_id uuid,
	n_overlap integer,
	odds_ratio double precision,
	pvalue double precision,
	adj_pvalue double precision
);


--
-- Name: TYPE enrich_result; Type: COMMENT; Schema: app_public_v2; Owner: -
--

COMMENT ON TYPE app_public_v2.enrich_result IS '@foreign key (gene_set_id) references app_public_v2.gene_set (id)';


--
-- Name: enriched_term_result; Type: TYPE; Schema: app_public_v2; Owner: -
--

CREATE TYPE app_public_v2.enriched_term_result AS (
	term character varying,
	count integer,
	odds_ratio double precision,
	pvalue double precision,
	adj_pvalue double precision,
	not_term_count integer,
	total_term_count integer,
	total_not_term_count integer
);


--
-- Name: enrichr_result; Type: TYPE; Schema: app_public_v2; Owner: -
--

CREATE TYPE app_public_v2.enrichr_result AS (
	term character varying,
	pvalue double precision,
	adj_pvalue double precision,
	count integer
);


--
-- Name: paginated_enrich_result; Type: TYPE; Schema: app_public_v2; Owner: -
--

CREATE TYPE app_public_v2.paginated_enrich_result AS (
	nodes app_public_v2.enrich_result[],
	total_count integer,
	enriched_terms character varying[],
	top_enriched_sigs character varying[]
);


--
-- Name: cell_line_results; Type: TYPE; Schema: public; Owner: -
--

CREATE TYPE public.cell_line_results AS (
	cell_line character varying,
	pcc double precision
);


--
-- Name: screen_target_results; Type: TYPE; Schema: public; Owner: -
--

CREATE TYPE public.screen_target_results AS (
	gene character varying,
	log2fc numeric,
	t double precision,
	p double precision
);


--
-- Name: welchs_t_test_results; Type: TYPE; Schema: public; Owner: -
--

CREATE TYPE public.welchs_t_test_results AS (
	t double precision,
	p double precision
);


--
-- Name: welchs_t_test_vectorized_results; Type: TYPE; Schema: public; Owner: -
--

CREATE TYPE public.welchs_t_test_vectorized_results AS (
	gene character varying,
	log2fc double precision,
	t double precision,
	p double precision,
	adj_p double precision
);


--
-- Name: welchs_t_test_vectorized_transcript_results; Type: TYPE; Schema: public; Owner: -
--

CREATE TYPE public.welchs_t_test_vectorized_transcript_results AS (
	gene character varying,
	transcript character varying,
	log2fc double precision,
	t double precision,
	p double precision,
	adj_p double precision
);


--
-- Name: enrich_functional_terms(character varying[], character varying, character varying, character varying[], bigint, jsonb); Type: FUNCTION; Schema: app_private_v2; Owner: -
--

CREATE FUNCTION app_private_v2.enrich_functional_terms(terms_concat character varying[], source_type character varying, species character varying, category_terms character varying[], total_count bigint, term_counts_json jsonb) RETURNS app_public_v2.enriched_term_result[]
    LANGUAGE plpython3u IMMUTABLE PARALLEL SAFE
    AS $$
    from scipy.stats import fisher_exact
    from statsmodels.stats.multitest import multipletests
    from collections import Counter
    import json
    import math

    if terms_concat is None or not terms_concat:
      return []

    concat_terms = list(filter(lambda t: t in category_terms, terms_concat))

    term_counter = Counter(concat_terms)
    term_counts = list(term_counter.items())
    total_enrich_term_count = sum(term_counter.values())
    term_counts_dict = json.loads(term_counts_json)

    results = []
    p_values = []
    for term, count in term_counts:
      if count < 5:
          continue
      try:
        escaped_term = term.replace("'", "''")
        total_term_count = term_counts_dict.get(term, 0)
        a = count
        b = total_enrich_term_count - count
        c = total_term_count
        d = total_count - total_term_count
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
$$;


SET default_tablespace = '';

SET default_table_access_method = heap;

--
-- Name: background; Type: TABLE; Schema: app_public_v2; Owner: -
--

CREATE TABLE app_public_v2.background (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    gene_ids jsonb NOT NULL,
    n_gene_ids integer NOT NULL,
    species character varying NOT NULL,
    created timestamp without time zone DEFAULT now()
);


--
-- Name: indexed_enrich(app_public_v2.background, uuid[], character varying, integer, double precision, double precision, integer, integer, double precision); Type: FUNCTION; Schema: app_private_v2; Owner: -
--

CREATE FUNCTION app_private_v2.indexed_enrich(background app_public_v2.background, gene_ids uuid[], filter_term character varying DEFAULT NULL::character varying, overlap_ge integer DEFAULT 1, pvalue_le double precision DEFAULT 0.05, adj_pvalue_le double precision DEFAULT 0.05, "offset" integer DEFAULT NULL::integer, first integer DEFAULT NULL::integer, filter_score_le double precision DEFAULT '-1'::integer) RETURNS app_public_v2.paginated_enrich_result
    LANGUAGE plpython3u IMMUTABLE PARALLEL SAFE
    AS $$
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
$$;


--
-- Name: sig_enrichr_terms(jsonb, character varying); Type: FUNCTION; Schema: app_private_v2; Owner: -
--

CREATE FUNCTION app_private_v2.sig_enrichr_terms(concat_terms jsonb, species character varying) RETURNS app_public_v2.enrichr_result[]
    LANGUAGE plpython3u IMMUTABLE PARALLEL SAFE
    AS $$
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
$$;


--
-- Name: user_gene_set; Type: TABLE; Schema: app_public_v2; Owner: -
--

CREATE TABLE app_public_v2.user_gene_set (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    genes character varying[],
    description character varying DEFAULT ''::character varying,
    created timestamp without time zone DEFAULT now() NOT NULL
);


--
-- Name: add_user_gene_set(character varying[], character varying); Type: FUNCTION; Schema: app_public_v2; Owner: -
--

CREATE FUNCTION app_public_v2.add_user_gene_set(genes character varying[], description character varying DEFAULT ''::character varying) RETURNS app_public_v2.user_gene_set
    LANGUAGE sql SECURITY DEFINER
    AS $$
  insert into app_public_v2.user_gene_set (genes, description)
  select
    (
      select array_agg(ug.gene order by ug.gene)
      from unnest(add_user_gene_set.genes) ug(gene)
    ) as genes,
    add_user_gene_set.description
  returning *;
$$;


--
-- Name: background_enrich(app_public_v2.background, character varying[], character varying, integer, double precision, double precision, integer, integer, double precision); Type: FUNCTION; Schema: app_public_v2; Owner: -
--

CREATE FUNCTION app_public_v2.background_enrich(background app_public_v2.background, genes character varying[], filter_term character varying DEFAULT NULL::character varying, overlap_ge integer DEFAULT 1, pvalue_le double precision DEFAULT 0.05, adj_pvalue_le double precision DEFAULT 0.05, "offset" integer DEFAULT NULL::integer, first integer DEFAULT NULL::integer, filter_score_le double precision DEFAULT '-1'::integer) RETURNS app_public_v2.paginated_enrich_result
    LANGUAGE sql IMMUTABLE SECURITY DEFINER PARALLEL SAFE
    AS $$
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
$$;


--
-- Name: background_overlap(app_public_v2.background, character varying[], integer); Type: FUNCTION; Schema: app_public_v2; Owner: -
--

CREATE FUNCTION app_public_v2.background_overlap(background app_public_v2.background, genes character varying[], overlap_greater_than integer DEFAULT 0) RETURNS TABLE(gene_set_id uuid, n_overlap_gene_ids integer, n_gs_gene_ids integer)
    LANGUAGE sql IMMUTABLE STRICT
    AS $$
  select
    gs.id as gene_set_id,
    count(ig.gene_id) as n_overlap_gene_ids,
    gs.n_gene_ids as n_gs_gene_ids
  from
    (
      select distinct g.gene_id::text
      from app_public_v2.gene_map(background_overlap.genes) g
    ) ig
    inner join app_public_v2.gene_set gs on gs.gene_ids ? ig.gene_id
  group by gs.id
  having count(ig.gene_id) > background_overlap.overlap_greater_than;
$$;


--
-- Name: current_background(); Type: FUNCTION; Schema: app_public_v2; Owner: -
--

CREATE FUNCTION app_public_v2.current_background() RETURNS app_public_v2.background
    LANGUAGE sql IMMUTABLE STRICT SECURITY DEFINER PARALLEL SAFE
    AS $$
  select *
  from app_public_v2.background
  order by created asc
  limit 1;
$$;


--
-- Name: gene_set; Type: TABLE; Schema: app_public_v2; Owner: -
--

CREATE TABLE app_public_v2.gene_set (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    term character varying NOT NULL,
    gene_ids jsonb NOT NULL,
    n_gene_ids integer NOT NULL,
    species character varying NOT NULL,
    created timestamp without time zone DEFAULT now() NOT NULL,
    created_at timestamp without time zone DEFAULT CURRENT_TIMESTAMP
);


--
-- Name: enrich_result_gene_set(app_public_v2.enrich_result); Type: FUNCTION; Schema: app_public_v2; Owner: -
--

CREATE FUNCTION app_public_v2.enrich_result_gene_set(enrich_result app_public_v2.enrich_result) RETURNS app_public_v2.gene_set
    LANGUAGE sql IMMUTABLE STRICT
    AS $$
  select gs.*
  from app_public_v2.gene_set gs
  where gs.id = enrich_result.gene_set_id;
$$;


--
-- Name: enriched_enrichr_terms(character varying[], character varying); Type: FUNCTION; Schema: app_public_v2; Owner: -
--

CREATE FUNCTION app_public_v2.enriched_enrichr_terms(enriched_sigs character varying[], organism character varying) RETURNS app_public_v2.enrichr_result[]
    LANGUAGE sql IMMUTABLE SECURITY DEFINER PARALLEL SAFE
    AS $$
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
$$;


--
-- Name: enriched_functional_terms(character varying[], character varying, character varying); Type: FUNCTION; Schema: app_public_v2; Owner: -
--

CREATE FUNCTION app_public_v2.enriched_functional_terms(enriched_terms character varying[], source_type character varying, spec character varying) RETURNS app_public_v2.enriched_term_result[]
    LANGUAGE sql IMMUTABLE SECURITY DEFINER
    AS $$
  with gse_ids as (
    -- Split enriched_terms vector to get GSE ids
    select regexp_replace(t, '\mGSE([^-]+)\M.*', 'GSE\1') as id
    from unnest(enriched_terms) as t
  ), gse_terms as (
    -- Get list of functional terms based on the source_type
    select array_agg(terms) as terms
    from (
      select unnest(llm_attrs) as terms
      from app_public_v2.gse_terms
      where gse in (select id from gse_ids) and species = spec
    ) subquery
    )
    -- Call the enrich_functional_terms function
    select * from app_private_v2.enrich_functional_terms(
      (select terms from gse_terms),
      source_type, spec,
      (SELECT array_agg(term_name )FROM app_public_v2.term_categories WHERE category = source_type),
      (SELECT term_total::bigint FROM app_public_v2.category_total_count WHERE category = source_type),
      (SELECT COALESCE(jsonb_object_agg(terms, term_count), '{}'::jsonb)
      FROM app_public_v2.terms_count_combined
      WHERE terms IN (SELECT unnest(terms) FROM gse_terms) AND organism = spec)
      );
$$;


--
-- Name: gene_map(character varying[]); Type: FUNCTION; Schema: app_public_v2; Owner: -
--

CREATE FUNCTION app_public_v2.gene_map(genes character varying[]) RETURNS TABLE(gene_id uuid, gene character varying)
    LANGUAGE sql IMMUTABLE STRICT PARALLEL SAFE
    AS $$
  select g.id as gene_id, ug.gene as gene
  from unnest(gene_map.genes) ug(gene)
  inner join app_public_v2.gene g on g.symbol = ug.gene or g.synonyms ? ug.gene;
$$;


--
-- Name: gene_set_gene_search(character varying[]); Type: FUNCTION; Schema: app_public_v2; Owner: -
--

CREATE FUNCTION app_public_v2.gene_set_gene_search(genes character varying[]) RETURNS SETOF app_public_v2.gene_set
    LANGUAGE sql IMMUTABLE STRICT PARALLEL SAFE
    AS $$
  select distinct gs.*
  from
    app_public_v2.gene_map(genes) g
    inner join app_public_v2.gene_set gs on gs.gene_ids ? g.gene_id::text;
$$;


--
-- Name: gene; Type: TABLE; Schema: app_public_v2; Owner: -
--

CREATE TABLE app_public_v2.gene (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    symbol character varying NOT NULL,
    synonyms jsonb DEFAULT '{}'::jsonb
);


--
-- Name: gene_set_genes(app_public_v2.gene_set); Type: FUNCTION; Schema: app_public_v2; Owner: -
--

CREATE FUNCTION app_public_v2.gene_set_genes(gene_set app_public_v2.gene_set) RETURNS SETOF app_public_v2.gene
    LANGUAGE sql IMMUTABLE STRICT PARALLEL SAFE
    AS $$
  select g.*
  from app_public_v2.gene g
  where gene_set_genes.gene_set.gene_ids ? g.id::text;
$$;


--
-- Name: gene_set_overlap(app_public_v2.gene_set, character varying[]); Type: FUNCTION; Schema: app_public_v2; Owner: -
--

CREATE FUNCTION app_public_v2.gene_set_overlap(gene_set app_public_v2.gene_set, genes character varying[]) RETURNS SETOF app_public_v2.gene
    LANGUAGE sql IMMUTABLE STRICT
    AS $$
  select distinct g.*
  from app_public_v2.gene_map(gene_set_overlap.genes) gm
  inner join app_public_v2.gene g on g.id = gm.gene_id
  where gene_set.gene_ids ? gm.gene_id::text;
$$;


--
-- Name: gse_info; Type: TABLE; Schema: app_public_v2; Owner: -
--

CREATE TABLE app_public_v2.gse_info (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    gse character varying,
    pmid character varying,
    title character varying,
    summary character varying,
    published_date date,
    species character varying,
    platform character varying,
    sample_groups jsonb,
    silhouette_score double precision,
    gse_attrs character varying
);


--
-- Name: gene_set_pmid; Type: MATERIALIZED VIEW; Schema: app_public_v2; Owner: -
--

CREATE MATERIALIZED VIEW app_public_v2.gene_set_pmid AS
 SELECT gs.id,
    gs.term,
    gse_info.id AS gse_id,
    gse_info.gse,
    gse_info.pmid,
    gse_info.title,
    gse_info.sample_groups,
    gse_info.platform,
    gse_info.published_date,
    gse_info.gse_attrs,
    gse_info.silhouette_score
   FROM (app_public_v2.gene_set gs
     JOIN app_public_v2.gse_info gse_info ON ((regexp_replace((gs.term)::text, '\mGSE([^-]+)\M.*'::text, 'GSE\1'::text) = (gse_info.gse)::text)))
  WITH NO DATA;


--
-- Name: MATERIALIZED VIEW gene_set_pmid; Type: COMMENT; Schema: app_public_v2; Owner: -
--

COMMENT ON MATERIALIZED VIEW app_public_v2.gene_set_pmid IS '@foreignKey (id) references app_public_v2.gene_set (id)';


--
-- Name: gene_set_term_search(character varying[]); Type: FUNCTION; Schema: app_public_v2; Owner: -
--

CREATE FUNCTION app_public_v2.gene_set_term_search(terms character varying[]) RETURNS SETOF app_public_v2.gene_set_pmid
    LANGUAGE sql IMMUTABLE STRICT PARALLEL SAFE
    AS $$
  select distinct gs.*
  from app_public_v2.gene_set_pmid gs
  inner join unnest(terms) ut(term) on gs.title ilike ('%' || ut.term || '%');
$$;


--
-- Name: gsm_meta; Type: TABLE; Schema: app_public_v2; Owner: -
--

CREATE TABLE app_public_v2.gsm_meta (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    gsm character varying NOT NULL,
    gse character varying,
    title character varying,
    characteristics_ch1 character varying,
    source_name_ch1 character varying
);


--
-- Name: get_gsm_meta(character varying[]); Type: FUNCTION; Schema: app_public_v2; Owner: -
--

CREATE FUNCTION app_public_v2.get_gsm_meta(gsms character varying[]) RETURNS SETOF app_public_v2.gsm_meta
    LANGUAGE sql IMMUTABLE STRICT PARALLEL SAFE
    AS $$
  select *
  from app_public_v2.gsm_meta
  where gsm = ANY (gsms);
$$;


--
-- Name: get_pb_info_by_ids(character varying[]); Type: FUNCTION; Schema: app_public_v2; Owner: -
--

CREATE FUNCTION app_public_v2.get_pb_info_by_ids(pmids character varying[]) RETURNS SETOF app_public_v2.gene_set_pmid
    LANGUAGE sql IMMUTABLE STRICT PARALLEL SAFE
    AS $$
  select *
  from app_public_v2.gene_set_pmid
  where pmid = ANY(pmids)
$$;


--
-- Name: pmid_info; Type: TABLE; Schema: app_public_v2; Owner: -
--

CREATE TABLE app_public_v2.pmid_info (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    pmid character varying NOT NULL,
    pmcid character varying,
    title character varying,
    pub_date character varying,
    doi character varying
);


--
-- Name: get_pb_meta_by_ids(character varying[]); Type: FUNCTION; Schema: app_public_v2; Owner: -
--

CREATE FUNCTION app_public_v2.get_pb_meta_by_ids(pmids character varying[]) RETURNS SETOF app_public_v2.pmid_info
    LANGUAGE sql IMMUTABLE STRICT PARALLEL SAFE
    AS $$
  select *
  from app_public_v2.pmid_info
  where pmid = ANY(pmids)
$$;


--
-- Name: pmc_info; Type: TABLE; Schema: app_public_v2; Owner: -
--

CREATE TABLE app_public_v2.pmc_info (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    pmcid character varying NOT NULL,
    title character varying,
    yr integer,
    doi character varying
);


--
-- Name: TABLE pmc_info; Type: COMMENT; Schema: app_public_v2; Owner: -
--

COMMENT ON TABLE app_public_v2.pmc_info IS '@foreignKey (pmcid) references app_public_v2.gene_set_pmc (pmc)';


--
-- Name: get_pmc_info_by_ids(character varying[]); Type: FUNCTION; Schema: app_public_v2; Owner: -
--

CREATE FUNCTION app_public_v2.get_pmc_info_by_ids(pmcids character varying[]) RETURNS SETOF app_public_v2.pmc_info
    LANGUAGE sql IMMUTABLE STRICT PARALLEL SAFE
    AS $$
  select *
  from app_public_v2.pmc_info
  where pmcid = ANY (pmcIds);
$$;


--
-- Name: release; Type: TABLE; Schema: app_public_v2; Owner: -
--

CREATE TABLE app_public_v2.release (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    n_publications_processed bigint,
    created timestamp without time zone DEFAULT now()
);


--
-- Name: pmc_stats; Type: MATERIALIZED VIEW; Schema: app_private_v2; Owner: -
--

CREATE MATERIALIZED VIEW app_private_v2.pmc_stats AS
 SELECT sum(release.n_publications_processed) AS n_publications_processed
   FROM app_public_v2.release
  WITH NO DATA;


--
-- Name: pmc_stats(); Type: FUNCTION; Schema: app_public_v2; Owner: -
--

CREATE FUNCTION app_public_v2.pmc_stats() RETURNS app_private_v2.pmc_stats
    LANGUAGE sql IMMUTABLE STRICT SECURITY DEFINER PARALLEL SAFE
    AS $$
  select * from app_private_v2.pmc_stats;
$$;


--
-- Name: terms_pmcs_count(character varying[]); Type: FUNCTION; Schema: app_public_v2; Owner: -
--

CREATE FUNCTION app_public_v2.terms_pmcs_count(pmcids character varying[]) RETURNS TABLE(pmc character varying, term character varying, id uuid, count integer)
    LANGUAGE sql IMMUTABLE STRICT PARALLEL SAFE
    AS $$
  select gsp.pmc, gs.term, gs.id, gs.n_gene_ids as count
  from
    app_public_v2.gene_set_pmc as gsp
    inner join app_public_v2.gene_set as gs on gs.id = gsp.id
  where gsp.pmc = ANY (pmcids);
$$;


--
-- Name: agg_mean_std_count(jsonb); Type: FUNCTION; Schema: public; Owner: -
--

CREATE FUNCTION public.agg_mean_std_count(stats jsonb) RETURNS jsonb
    LANGUAGE plpython3u IMMUTABLE
    AS $$
  if not stats: return None
  import numpy as np
  import pandas as pd
  import json

  def combined_mean(means, ns):
    if ns.sum() == 0: return float('nan')
    return (means * ns).sum() / ns.sum()

  def combined_std(stds, ns, means):
    if ns.sum() == 0: return float('nan')
    mean = combined_mean(means, ns)
    ds = means - mean
    return (((ns*stds**2).sum() + (ns*ds**2).sum()) / ns.sum())**(1/2)

  stats_df = pd.DataFrame(json.loads(stats)).dropna()
  return pd.Series({
    'mean': combined_mean(stats_df['mean'], stats_df['count']),
    'std': combined_std(stats_df['std'], stats_df['count'], stats_df['mean']),
    'count': stats_df['count'].sum(),
  }).to_json()
$$;


--
-- Name: aggregate_stats(jsonb); Type: FUNCTION; Schema: public; Owner: -
--

CREATE FUNCTION public.aggregate_stats(stats jsonb) RETURNS jsonb
    LANGUAGE sql IMMUTABLE
    AS $$
  select jsonb_build_object(
    'mean',
    (
      select to_jsonb(avg(j.value::text::double precision))
      from jsonb_each(stats->'mean') j
    ),
    'std',
    (
      select to_jsonb(sqrt(sum(power(j.value::text::double precision,2))))
      from jsonb_each(stats->'std') j
    ),
    'count',
    (
      select to_jsonb(sum(j.value::text::double precision))
      from jsonb_each(stats->'count') j
    )
  )
$$;


--
-- Name: data; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.data (
    id integer NOT NULL,
    database uuid NOT NULL,
    gene uuid NOT NULL,
    "values" jsonb NOT NULL
);


--
-- Name: data_transcript; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.data_transcript (
    id integer NOT NULL,
    database uuid NOT NULL,
    transcript uuid NOT NULL,
    "values" jsonb NOT NULL
);


--
-- Name: database; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.database (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    dbname character varying,
    created timestamp without time zone DEFAULT now()
);


--
-- Name: gene; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.gene (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    gene character varying
);


--
-- Name: transcript; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.transcript (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    transcript character varying
);


--
-- Name: data_complete; Type: VIEW; Schema: public; Owner: -
--

CREATE VIEW public.data_complete AS
 SELECT database.dbname,
    gene.gene,
    COALESCE(data."values", data_transcript."values") AS "values",
    transcript.transcript
   FROM ((((public.database
     LEFT JOIN public.data ON ((data.database = database.id)))
     LEFT JOIN public.gene ON ((data.gene = gene.id)))
     LEFT JOIN public.data_transcript ON ((data_transcript.database = database.id)))
     LEFT JOIN public.transcript ON ((data_transcript.transcript = transcript.id)));


--
-- Name: cell_line_ppc_vectorized(double precision[], public.data_complete[]); Type: FUNCTION; Schema: public; Owner: -
--

CREATE FUNCTION public.cell_line_ppc_vectorized(a_mean double precision[], input_background public.data_complete[]) RETURNS public.cell_line_results[]
    LANGUAGE plpython3u IMMUTABLE
    AS $$
  import numpy as np
  import json

  cell_lines_dict = {}

  gene_order = []
  print(input_background[0].keys())
  for d in input_background:
    input_bg_data = json.loads(d['values'])
    expr = input_bg_data['value']
    gene = d['gene']
    for cl in expr:
      if cl not in cell_lines_dict:
        cell_lines_dict[cl] = []
      cell_lines_dict[cl].append(expr[cl])

  result = []
  for cell_line in cell_lines_dict:
      result.append((cell_line, np.corrcoef(np.log(np.array(a_mean) + 1), np.log(np.array(cell_lines_dict[cell_line]) + 1))[0][1]))

  return result
$$;


--
-- Name: cell_line_ppc_vectorized(double precision[], public.data_complete[], character varying); Type: FUNCTION; Schema: public; Owner: -
--

CREATE FUNCTION public.cell_line_ppc_vectorized(a_mean double precision[], input_background public.data_complete[], cell_line_target character varying) RETURNS public.cell_line_results[]
    LANGUAGE plpython3u IMMUTABLE
    AS $$
  import numpy as np
  import json

  cell_lines_dict = {}

  gene_order = []
  print(input_background[0].keys())
  for d in input_background:
    input_bg_data = json.loads(d['values'])
    expr = input_bg_data['value']
    gene = d['gene']
    for cl in expr:
      if cell_line_target == 'All':
        if cl not in cell_lines_dict:
            cell_lines_dict[cl] = []
        cell_lines_dict[cl].append(expr[cl])
      else:
        if cl.split(' - ')[0] == cell_line_target:
          if cl not in cell_lines_dict:
            cell_lines_dict[cl] = []
          cell_lines_dict[cl].append(expr[cl])

  result = []
  for cell_line in cell_lines_dict:
      result.append((cell_line, np.corrcoef(np.log(np.array(a_mean) + 1), np.log(np.array(cell_lines_dict[cell_line]) + 1))[0][1]))

  return result
$$;


--
-- Name: safe_log2fc(numeric, numeric); Type: FUNCTION; Schema: public; Owner: -
--

CREATE FUNCTION public.safe_log2fc(b numeric, a numeric) RETURNS numeric
    LANGUAGE sql IMMUTABLE
    AS $$
  select case
    when b = 0 and a = 0 then 0.0
    when b = 0 and a <> 0 then 'infinity'::numeric
    when b <> 0 and a = 0 then '-infinity'::numeric
    else log(2, b) - log(2, a)
  end
$$;


--
-- Name: screen_cell_lines_vectorized(jsonb, character varying, character varying); Type: FUNCTION; Schema: public; Owner: -
--

CREATE FUNCTION public.screen_cell_lines_vectorized(input_data jsonb, background character varying, cell_line_target character varying) RETURNS SETOF public.cell_line_results
    LANGUAGE sql IMMUTABLE
    AS $$
  with
  input_data_each as (
    select
      j.key as gene,
      j.value || jsonb_build_object('count', input_data->'n') as values
    from jsonb_each(input_data->'genes') as j
  ),
  input_background as (
    select
      input_data_each.gene,
      input_data_each.values as input_data,
      data_complete as background_data,
      row_number() over () / 1000 as group_number  -- Change 1000 to your desired chunk size
    from input_data_each
    inner join gene on gene.gene = input_data_each.gene
    inner join data_complete on data_complete.gene = gene.gene
    where data_complete.dbname = background
  ),
  vectorized_stats as (
    select
      group_number,
      cell_line_ppc_vectorized(
        array_agg((input_data->>'mean')::double precision),
        array_agg((background_data)::data_complete),
        cell_line_target
      ) as value
    from input_background
    group by group_number
  ),
  stats as (
    select r.*
    from
      vectorized_stats,
      unnest(vectorized_stats.value) r
    where
      r.pcc >= 0
  )
select *
from stats
$$;


--
-- Name: screen_targets(jsonb, uuid); Type: FUNCTION; Schema: public; Owner: -
--

CREATE FUNCTION public.screen_targets(input_data jsonb, background uuid) RETURNS SETOF public.screen_target_results
    LANGUAGE sql IMMUTABLE
    AS $$
  with
  input_data_each as (
    select
      j.key as gene,
      j.value || jsonb_build_object('count', input_data->'n') as values
    from jsonb_each(input_data->'genes') as j
  ),
  input_background as (
    select
      input_data_each.gene,
      input_data_each.values as input_data,
      database_agg.values as background_data
    from input_data_each
    inner join gene on gene.gene = input_data_each.gene
    inner join database_agg on database_agg.gene = gene.id
    where database_agg.database = background
  ),
  stats as (
    select
      gene,
      safe_log2fc((input_data->>'mean')::numeric, (background_data->>'mean')::numeric) as log2fc,
      r.t as t,
      r.p as p
    from input_background
    cross join lateral welchs_t_test(
        (input_data->>'mean')::double precision,
        (input_data->>'std')::double precision,
        (input_data->>'count')::double precision,
        (background_data->>'mean')::double precision,
        (background_data->>'std')::double precision,
        (background_data->>'count')::double precision,
        false,
        'two-sided'
    ) as r
    where
      r.t > 0
  )
select *
from stats
$$;


--
-- Name: screen_targets_transcript_vectorized(jsonb, uuid); Type: FUNCTION; Schema: public; Owner: -
--

CREATE FUNCTION public.screen_targets_transcript_vectorized(input_data jsonb, background uuid) RETURNS SETOF public.welchs_t_test_vectorized_transcript_results
    LANGUAGE sql IMMUTABLE
    AS $$
  with
  input_data_each as (
    select
      j.key as transcript,
      j.value || jsonb_build_object('count', input_data->'n') as values
    from jsonb_each(input_data->'transcripts') as j
  ),
  input_background as (
    select
      input_data_each.transcript,
      input_data_each.values as input_data,
      database_agg_transcript.values as background_data
    from input_data_each
    inner join transcript on transcript.transcript = input_data_each.transcript
    inner join database_agg_transcript on database_agg_transcript.transcript = transcript.id
    where database_agg_transcript.database = background
  ),
  vectorized_stats as (
    select
      welchs_t_test_vectorized(
        array_agg(transcript),
        array_agg((input_data->>'mean')::double precision),
        array_agg((input_data->>'std')::double precision),
        array_agg((input_data->>'count')::double precision),
        array_agg((background_data->>'mean')::double precision),
        array_agg((background_data->>'std')::double precision),
        array_agg((background_data->>'count')::double precision),
        false,
        'two-sided',
        0.05,
        'fdr_bh'
      ) as value
    from input_background
  ),
  stats as (
    select r.*
    from
      vectorized_stats,
      unnest(vectorized_stats.value) r
    where
      r.t > 0
  ),
  stats_mapped as (
    select
      mapper.gene,
      stats.gene as transcript,
      stats.log2fc,
      stats.t,
      stats.p,
      stats.adj_p
    from stats
    left join mapper on stats.gene = mapper.transcript
  )
select *
from stats_mapped
$$;


--
-- Name: screen_targets_vectorized(jsonb, uuid); Type: FUNCTION; Schema: public; Owner: -
--

CREATE FUNCTION public.screen_targets_vectorized(input_data jsonb, background uuid) RETURNS SETOF public.welchs_t_test_vectorized_results
    LANGUAGE sql IMMUTABLE
    AS $$
  with
  input_data_each as (
    select
      j.key as gene,
      j.value || jsonb_build_object('count', input_data->'n') as values
    from jsonb_each(input_data->'genes') as j
  ),
  input_background as (
    select
      input_data_each.gene,
      input_data_each.values as input_data,
      database_agg.values as background_data
    from input_data_each
    inner join gene on gene.gene = input_data_each.gene
    inner join database_agg on database_agg.gene = gene.id
    where database_agg.database = background
  ),
  vectorized_stats as (
    select
      welchs_t_test_vectorized(
        array_agg(gene),
        array_agg((input_data->>'mean')::double precision),
        array_agg((input_data->>'std')::double precision),
        array_agg((input_data->>'count')::double precision),
        array_agg((background_data->>'mean')::double precision),
        array_agg((background_data->>'std')::double precision),
        array_agg((background_data->>'count')::double precision),
        false,
        'two-sided',
        0.05,
        'fdr_bh'
      ) as value
    from input_background
  ),
  stats as (
    select r.*
    from
      vectorized_stats,
      unnest(vectorized_stats.value) r
    where
      r.t > 0
  )
select *
from stats
$$;


--
-- Name: unified_data_complete(character varying, character varying[]); Type: FUNCTION; Schema: public; Owner: -
--

CREATE FUNCTION public.unified_data_complete(gene character varying, dbnames character varying[]) RETURNS jsonb
    LANGUAGE sql IMMUTABLE
    AS $_$
with cte as (
  select
    db.dbname,
    dtv."desc" as "desc",
    jsonb_object_agg(
      utm."term",
      ((tv."value"::double precision) - (dd.dist->>'mean')::double precision) / (dd.dist->>'std')::double precision
    ) as "term_value"
  from
    unnest(dbnames) dbn (dbname)
    inner join database db on db.dbname = dbn.dbname
    inner join database_dist dd on dd.database = db.id
    inner join data d on d.database = db.id
    inner join gene g on d.gene = g.id,
    jsonb_each(d."values") as dtv("desc", "term_value"),
    jsonb_each_text(dtv."term_value") as tv("term", "value")
    inner join unified_term_map utm on tv."term" = utm."original_term"
  where g.gene = $1
  group by db.dbname, dtv."desc"
), cte2 as (
  select cte."dbname", jsonb_object_agg(cte."desc", cte."term_value") as "values"
  from cte
  group by dbname
) select jsonb_object_agg("dbname", "values") from cte2;
$_$;


--
-- Name: welchs_t_test(double precision, double precision, double precision, double precision, double precision, double precision, boolean, character varying); Type: FUNCTION; Schema: public; Owner: -
--

CREATE FUNCTION public.welchs_t_test(a_mean double precision, a_std double precision, a_n double precision, b_mean double precision, b_std double precision, b_n double precision, equal_var boolean, alternative character varying) RETURNS public.welchs_t_test_results
    LANGUAGE plpython3u IMMUTABLE
    AS $$
  import scipy.stats
  if a_mean is None or not a_std or not a_n or b_mean is None or not b_std or not b_n:
    return (None, None)
  result = scipy.stats.ttest_ind_from_stats(a_mean, a_std, a_n, b_mean, b_std, b_n, equal_var=equal_var, alternative=alternative)
  return (result.statistic, result.pvalue)
$$;


--
-- Name: welchs_t_test_vectorized(character varying[], double precision[], double precision[], double precision[], double precision[], double precision[], double precision[], boolean, character varying, double precision, character varying); Type: FUNCTION; Schema: public; Owner: -
--

CREATE FUNCTION public.welchs_t_test_vectorized(genes character varying[], a_mean double precision[], a_std double precision[], a_n double precision[], b_mean double precision[], b_std double precision[], b_n double precision[], ttest_equal_var boolean, ttest_alternative character varying, p_adjust_alpha double precision, p_adjust_method character varying) RETURNS public.welchs_t_test_vectorized_results[]
    LANGUAGE plpython3u IMMUTABLE
    AS $$
  import numpy as np
  import scipy.stats
  import statsmodels.stats.multitest

  n = len(genes)
  np_a_mean = np.array(a_mean, dtype=np.float64)
  np_a_std = np.array(a_std, dtype=np.float64)
  np_a_n = np.array(a_n, dtype=np.float64)
  np_b_mean = np.array(b_mean, dtype=np.float64)
  np_b_std = np.array(b_std, dtype=np.float64)
  np_b_n = np.array(b_n, dtype=np.float64)

  np_a_mean_is_zero = np.isclose(np_a_mean, 0.)
  np_b_mean_is_zero = np.isclose(np_b_mean, 0.)
  np_a_and_b_are_nonzero = ~np_a_mean_is_zero & ~np_b_mean_is_zero
  log2fc = np.zeros(n)
  log2fc[np_a_mean_is_zero & np_b_mean_is_zero] = 0.
  log2fc[np_a_mean_is_zero & ~np_b_mean_is_zero] = -np.inf
  log2fc[~np_a_mean_is_zero & np_b_mean_is_zero] = np.inf
  log2fc[np_a_and_b_are_nonzero] = (
    np.log2(np_a_mean[np_a_and_b_are_nonzero])
    - np.log2(np_b_mean[np_a_and_b_are_nonzero])
  )

  mask = ~(
    np.isnan(np_a_mean)
    |np.isnan(np_a_std)|np.isclose(np_a_std, 0.)
    |np.isnan(np_a_n)|np.isclose(np_a_n, 0.)
    |np.isnan(np_b_mean)
    |np.isnan(np_b_std)|np.isclose(np_b_std, 0.)
    |np.isnan(np_b_n)|np.isclose(np_b_n, 0.)
  )
  result = scipy.stats.ttest_ind_from_stats(
    np_a_mean[mask], np_a_std[mask], np_a_n[mask],
    np_b_mean[mask], np_b_std[mask], np_b_n[mask],
    equal_var=ttest_equal_var, alternative=ttest_alternative,
  )
  tstats = result.statistic
  pvals = result.pvalue

  try:
    reject, adj_pvals, alphacSidak, alphacBonf = statsmodels.stats.multitest.multipletests(
      pvals,
      p_adjust_alpha,
      p_adjust_method,
    )
  except:
    adj_pvals = np.nan

  results = np.zeros((n, 4))
  results[~mask, :] = np.nan
  results[:, 0] = log2fc
  results[mask, 1] = tstats
  results[mask, 2] = pvals
  results[mask, 3] = adj_pvals
  return [
    (gene, log2fc, t, p, adj_p)
    for gene, (log2fc, t, p, adj_p) in zip(genes, results)
  ]
$$;


--
-- Name: gse_terms; Type: TABLE; Schema: app_public_v2; Owner: -
--

CREATE TABLE app_public_v2.gse_terms (
    gse character varying NOT NULL,
    species character varying NOT NULL,
    llm_attrs character varying[],
    pubmed_attrs character varying[],
    mesh_attrs character varying[]
);


--
-- Name: term_categories; Type: TABLE; Schema: app_public_v2; Owner: -
--

CREATE TABLE app_public_v2.term_categories (
    term_name character varying NOT NULL,
    category character varying NOT NULL
);


--
-- Name: terms_count_combined; Type: MATERIALIZED VIEW; Schema: app_public_v2; Owner: -
--

CREATE MATERIALIZED VIEW app_public_v2.terms_count_combined AS
 SELECT subquery.terms,
    count(subquery.terms) AS term_count,
    'mouse'::text AS organism
   FROM ( SELECT unnest(gse_terms.llm_attrs) AS terms
           FROM app_public_v2.gse_terms
          WHERE ((gse_terms.species)::text = 'mouse'::text)) subquery
  GROUP BY subquery.terms
UNION ALL
 SELECT subquery.terms,
    count(subquery.terms) AS term_count,
    'human'::text AS organism
   FROM ( SELECT unnest(gse_terms.llm_attrs) AS terms
           FROM app_public_v2.gse_terms
          WHERE ((gse_terms.species)::text = 'human'::text)) subquery
  GROUP BY subquery.terms
  WITH NO DATA;


--
-- Name: category_total_count; Type: MATERIALIZED VIEW; Schema: app_public_v2; Owner: -
--

CREATE MATERIALIZED VIEW app_public_v2.category_total_count AS
 SELECT tc.category,
    sum(tcc.term_count) AS term_total
   FROM (app_public_v2.term_categories tc
     JOIN app_public_v2.terms_count_combined tcc ON (((tc.term_name)::text = (tcc.terms)::text)))
  GROUP BY tc.category
  WITH NO DATA;


--
-- Name: enrichr_terms; Type: TABLE; Schema: app_public_v2; Owner: -
--

CREATE TABLE app_public_v2.enrichr_terms (
    sig character varying NOT NULL,
    organism character varying NOT NULL,
    sig_terms character varying[],
    enrichr_stats jsonb
);


--
-- Name: gene_set_gse; Type: MATERIALIZED VIEW; Schema: app_public_v2; Owner: -
--

CREATE MATERIALIZED VIEW app_public_v2.gene_set_gse AS
 SELECT gs.id,
    regexp_replace((gs.term)::text, '\mGSE([^-]+)\M.*'::text, 'GSE\1'::text) AS gse,
    gs.species
   FROM app_public_v2.gene_set gs
  WITH NO DATA;


--
-- Name: MATERIALIZED VIEW gene_set_gse; Type: COMMENT; Schema: app_public_v2; Owner: -
--

COMMENT ON MATERIALIZED VIEW app_public_v2.gene_set_gse IS '@foreignKey (id) references app_public_v2.gene_set (id)';


--
-- Name: gene_set_pmc; Type: MATERIALIZED VIEW; Schema: app_public_v2; Owner: -
--

CREATE MATERIALIZED VIEW app_public_v2.gene_set_pmc AS
 SELECT gs.id,
    regexp_replace((gs.term)::text, '^(^PMC\d+)(.*)$'::text, '\1'::text) AS pmc
   FROM app_public_v2.gene_set gs
  WITH NO DATA;


--
-- Name: MATERIALIZED VIEW gene_set_pmc; Type: COMMENT; Schema: app_public_v2; Owner: -
--

COMMENT ON MATERIALIZED VIEW app_public_v2.gene_set_pmc IS '@foreignKey (id) references app_public_v2.gene_set (id)';


--
-- Name: gse; Type: VIEW; Schema: app_public_v2; Owner: -
--

CREATE VIEW app_public_v2.gse AS
 SELECT DISTINCT gene_set_gse.gse
   FROM app_public_v2.gene_set_gse;


--
-- Name: VIEW gse; Type: COMMENT; Schema: app_public_v2; Owner: -
--

COMMENT ON VIEW app_public_v2.gse IS '@foreignKey (gse) references app_public_v2.gene_set_gse (gse)';


--
-- Name: pmc; Type: VIEW; Schema: app_public_v2; Owner: -
--

CREATE VIEW app_public_v2.pmc AS
 SELECT DISTINCT gene_set_pmc.pmc
   FROM app_public_v2.gene_set_pmc;


--
-- Name: VIEW pmc; Type: COMMENT; Schema: app_public_v2; Owner: -
--

COMMENT ON VIEW app_public_v2.pmc IS '@foreignKey (pmc) references app_public_v2.gene_set_pmc (pmc)';


--
-- Name: data_id_seq; Type: SEQUENCE; Schema: public; Owner: -
--

CREATE SEQUENCE public.data_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


--
-- Name: data_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: -
--

ALTER SEQUENCE public.data_id_seq OWNED BY public.data.id;


--
-- Name: data_transcript_id_seq; Type: SEQUENCE; Schema: public; Owner: -
--

CREATE SEQUENCE public.data_transcript_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


--
-- Name: data_transcript_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: -
--

ALTER SEQUENCE public.data_transcript_id_seq OWNED BY public.data_transcript.id;


--
-- Name: database_agg; Type: MATERIALIZED VIEW; Schema: public; Owner: -
--

CREATE MATERIALIZED VIEW public.database_agg AS
 SELECT data.id,
    data.database,
    data.gene,
    public.aggregate_stats(data."values") AS "values"
   FROM public.data
  WITH NO DATA;


--
-- Name: database_agg_transcript; Type: MATERIALIZED VIEW; Schema: public; Owner: -
--

CREATE MATERIALIZED VIEW public.database_agg_transcript AS
 SELECT data_transcript.id,
    data_transcript.database,
    data_transcript.transcript,
    public.aggregate_stats(data_transcript."values") AS "values"
   FROM public.data_transcript
  WITH NO DATA;


--
-- Name: database_dist; Type: MATERIALIZED VIEW; Schema: public; Owner: -
--

CREATE MATERIALIZED VIEW public.database_dist AS
 SELECT database_agg.database,
    public.agg_mean_std_count(jsonb_agg(database_agg."values")) AS dist
   FROM public.database_agg
  GROUP BY database_agg.database
  WITH NO DATA;


--
-- Name: depmap_cell_line_mapping; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.depmap_cell_line_mapping (
    cell_line character varying,
    depmap_id character varying
);


--
-- Name: gene_info; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.gene_info (
    id integer NOT NULL,
    symbol character varying,
    synonyms jsonb,
    chromosome character varying,
    map_location character varying,
    description character varying,
    type_of_gene character varying,
    summary text
);


--
-- Name: gene_transcript; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.gene_transcript (
    gene uuid NOT NULL,
    transcript uuid NOT NULL
);


--
-- Name: mapper; Type: VIEW; Schema: public; Owner: -
--

CREATE VIEW public.mapper AS
 SELECT gene.gene,
    transcript.transcript
   FROM ((public.gene_transcript
     LEFT JOIN public.gene ON ((gene_transcript.gene = gene.id)))
     LEFT JOIN public.transcript ON ((gene_transcript.transcript = transcript.id)));


--
-- Name: schema_migrations; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.schema_migrations (
    version character varying(255) NOT NULL
);


--
-- Name: unified_term_map; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.unified_term_map (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    original_term character varying,
    term character varying
);


--
-- Name: data id; Type: DEFAULT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.data ALTER COLUMN id SET DEFAULT nextval('public.data_id_seq'::regclass);


--
-- Name: data_transcript id; Type: DEFAULT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.data_transcript ALTER COLUMN id SET DEFAULT nextval('public.data_transcript_id_seq'::regclass);


--
-- Name: background background_pkey; Type: CONSTRAINT; Schema: app_public_v2; Owner: -
--

ALTER TABLE ONLY app_public_v2.background
    ADD CONSTRAINT background_pkey PRIMARY KEY (id);


--
-- Name: enrichr_terms enrichr_terms_sig_organism_key; Type: CONSTRAINT; Schema: app_public_v2; Owner: -
--

ALTER TABLE ONLY app_public_v2.enrichr_terms
    ADD CONSTRAINT enrichr_terms_sig_organism_key UNIQUE (sig, organism);


--
-- Name: gene gene_pkey; Type: CONSTRAINT; Schema: app_public_v2; Owner: -
--

ALTER TABLE ONLY app_public_v2.gene
    ADD CONSTRAINT gene_pkey PRIMARY KEY (id);


--
-- Name: gene_set gene_set_pkey; Type: CONSTRAINT; Schema: app_public_v2; Owner: -
--

ALTER TABLE ONLY app_public_v2.gene_set
    ADD CONSTRAINT gene_set_pkey PRIMARY KEY (id);


--
-- Name: gene_set gene_set_term_key; Type: CONSTRAINT; Schema: app_public_v2; Owner: -
--

ALTER TABLE ONLY app_public_v2.gene_set
    ADD CONSTRAINT gene_set_term_key UNIQUE (term);


--
-- Name: gene gene_symbol_key; Type: CONSTRAINT; Schema: app_public_v2; Owner: -
--

ALTER TABLE ONLY app_public_v2.gene
    ADD CONSTRAINT gene_symbol_key UNIQUE (symbol);


--
-- Name: gse_info gse_info_pkey; Type: CONSTRAINT; Schema: app_public_v2; Owner: -
--

ALTER TABLE ONLY app_public_v2.gse_info
    ADD CONSTRAINT gse_info_pkey PRIMARY KEY (id);


--
-- Name: gse_terms gse_terms_gse_species_key; Type: CONSTRAINT; Schema: app_public_v2; Owner: -
--

ALTER TABLE ONLY app_public_v2.gse_terms
    ADD CONSTRAINT gse_terms_gse_species_key UNIQUE (gse, species);


--
-- Name: gsm_meta gsm_meta_gsm_key; Type: CONSTRAINT; Schema: app_public_v2; Owner: -
--

ALTER TABLE ONLY app_public_v2.gsm_meta
    ADD CONSTRAINT gsm_meta_gsm_key UNIQUE (gsm);


--
-- Name: gsm_meta gsm_meta_pkey; Type: CONSTRAINT; Schema: app_public_v2; Owner: -
--

ALTER TABLE ONLY app_public_v2.gsm_meta
    ADD CONSTRAINT gsm_meta_pkey PRIMARY KEY (id);


--
-- Name: pmc_info pmc_info_pkey; Type: CONSTRAINT; Schema: app_public_v2; Owner: -
--

ALTER TABLE ONLY app_public_v2.pmc_info
    ADD CONSTRAINT pmc_info_pkey PRIMARY KEY (id);


--
-- Name: pmc_info pmc_info_pmcid_key; Type: CONSTRAINT; Schema: app_public_v2; Owner: -
--

ALTER TABLE ONLY app_public_v2.pmc_info
    ADD CONSTRAINT pmc_info_pmcid_key UNIQUE (pmcid);


--
-- Name: pmid_info pmid_info_pkey; Type: CONSTRAINT; Schema: app_public_v2; Owner: -
--

ALTER TABLE ONLY app_public_v2.pmid_info
    ADD CONSTRAINT pmid_info_pkey PRIMARY KEY (id);


--
-- Name: pmid_info pmid_info_pmid_key; Type: CONSTRAINT; Schema: app_public_v2; Owner: -
--

ALTER TABLE ONLY app_public_v2.pmid_info
    ADD CONSTRAINT pmid_info_pmid_key UNIQUE (pmid);


--
-- Name: release release_pkey; Type: CONSTRAINT; Schema: app_public_v2; Owner: -
--

ALTER TABLE ONLY app_public_v2.release
    ADD CONSTRAINT release_pkey PRIMARY KEY (id);


--
-- Name: term_categories term_categories_term_name_category_key; Type: CONSTRAINT; Schema: app_public_v2; Owner: -
--

ALTER TABLE ONLY app_public_v2.term_categories
    ADD CONSTRAINT term_categories_term_name_category_key UNIQUE (term_name, category);


--
-- Name: user_gene_set user_gene_set_pkey; Type: CONSTRAINT; Schema: app_public_v2; Owner: -
--

ALTER TABLE ONLY app_public_v2.user_gene_set
    ADD CONSTRAINT user_gene_set_pkey PRIMARY KEY (id);


--
-- Name: data data_database_gene_key; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.data
    ADD CONSTRAINT data_database_gene_key UNIQUE (database, gene);


--
-- Name: data data_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.data
    ADD CONSTRAINT data_pkey PRIMARY KEY (id);


--
-- Name: data_transcript data_transcript_database_transcript_key; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.data_transcript
    ADD CONSTRAINT data_transcript_database_transcript_key UNIQUE (database, transcript);


--
-- Name: data_transcript data_transcript_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.data_transcript
    ADD CONSTRAINT data_transcript_pkey PRIMARY KEY (id);


--
-- Name: database database_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.database
    ADD CONSTRAINT database_pkey PRIMARY KEY (id);


--
-- Name: gene_info gene_info_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.gene_info
    ADD CONSTRAINT gene_info_pkey PRIMARY KEY (id);


--
-- Name: gene gene_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.gene
    ADD CONSTRAINT gene_pkey PRIMARY KEY (id);


--
-- Name: gene_transcript gene_transcript_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.gene_transcript
    ADD CONSTRAINT gene_transcript_pkey PRIMARY KEY (gene, transcript);


--
-- Name: schema_migrations schema_migrations_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.schema_migrations
    ADD CONSTRAINT schema_migrations_pkey PRIMARY KEY (version);


--
-- Name: transcript transcript_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.transcript
    ADD CONSTRAINT transcript_pkey PRIMARY KEY (id);


--
-- Name: unified_term_map unified_term_map_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.unified_term_map
    ADD CONSTRAINT unified_term_map_pkey PRIMARY KEY (id);


--
-- Name: background_gene_ids_idx; Type: INDEX; Schema: app_public_v2; Owner: -
--

CREATE INDEX background_gene_ids_idx ON app_public_v2.background USING gin (gene_ids);


--
-- Name: gene_set_gene_ids_idx; Type: INDEX; Schema: app_public_v2; Owner: -
--

CREATE INDEX gene_set_gene_ids_idx ON app_public_v2.gene_set USING gin (gene_ids);


--
-- Name: gene_set_gse_gse_idx; Type: INDEX; Schema: app_public_v2; Owner: -
--

CREATE INDEX gene_set_gse_gse_idx ON app_public_v2.gene_set_gse USING btree (gse);


--
-- Name: gene_set_gse_id_gse_idx; Type: INDEX; Schema: app_public_v2; Owner: -
--

CREATE UNIQUE INDEX gene_set_gse_id_gse_idx ON app_public_v2.gene_set_gse USING btree (id, gse);


--
-- Name: gene_set_gse_id_idx; Type: INDEX; Schema: app_public_v2; Owner: -
--

CREATE INDEX gene_set_gse_id_idx ON app_public_v2.gene_set_gse USING btree (id);


--
-- Name: gene_set_pmc_id_idx; Type: INDEX; Schema: app_public_v2; Owner: -
--

CREATE INDEX gene_set_pmc_id_idx ON app_public_v2.gene_set_pmc USING btree (id);


--
-- Name: gene_set_pmc_id_pmc_idx; Type: INDEX; Schema: app_public_v2; Owner: -
--

CREATE UNIQUE INDEX gene_set_pmc_id_pmc_idx ON app_public_v2.gene_set_pmc USING btree (id, pmc);


--
-- Name: gene_set_pmc_pmc_idx; Type: INDEX; Schema: app_public_v2; Owner: -
--

CREATE INDEX gene_set_pmc_pmc_idx ON app_public_v2.gene_set_pmc USING btree (pmc);


--
-- Name: gene_set_pmid_id_idx; Type: INDEX; Schema: app_public_v2; Owner: -
--

CREATE INDEX gene_set_pmid_id_idx ON app_public_v2.gene_set_pmid USING btree (id);


--
-- Name: gene_set_pmid_pmid_idx; Type: INDEX; Schema: app_public_v2; Owner: -
--

CREATE INDEX gene_set_pmid_pmid_idx ON app_public_v2.gene_set_pmid USING btree (pmid);


--
-- Name: gene_set_term_idx; Type: INDEX; Schema: app_public_v2; Owner: -
--

CREATE INDEX gene_set_term_idx ON app_public_v2.gene_set USING gin (term public.gin_trgm_ops);


--
-- Name: gene_synonyms_idx; Type: INDEX; Schema: app_public_v2; Owner: -
--

CREATE INDEX gene_synonyms_idx ON app_public_v2.gene USING gin (synonyms);


--
-- Name: idx_gene_set_species; Type: INDEX; Schema: app_public_v2; Owner: -
--

CREATE INDEX idx_gene_set_species ON app_public_v2.gene_set USING btree (species);


--
-- Name: idx_gsm_meta_gsm; Type: INDEX; Schema: app_public_v2; Owner: -
--

CREATE INDEX idx_gsm_meta_gsm ON app_public_v2.gsm_meta USING btree (gsm);


--
-- Name: pmid_info_pmid_idx; Type: INDEX; Schema: app_public_v2; Owner: -
--

CREATE INDEX pmid_info_pmid_idx ON app_public_v2.pmid_info USING btree (pmid);


--
-- Name: release_created_idx; Type: INDEX; Schema: app_public_v2; Owner: -
--

CREATE INDEX release_created_idx ON app_public_v2.release USING btree (created);


--
-- Name: data_database_fkey; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX data_database_fkey ON public.data USING btree (database);


--
-- Name: data_gene_fkey; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX data_gene_fkey ON public.data USING btree (gene);


--
-- Name: database_agg_database_idx; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX database_agg_database_idx ON public.database_agg USING btree (database);


--
-- Name: database_agg_id_idx; Type: INDEX; Schema: public; Owner: -
--

CREATE UNIQUE INDEX database_agg_id_idx ON public.database_agg USING btree (id);


--
-- Name: database_agg_trancript_database_idx; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX database_agg_trancript_database_idx ON public.database_agg_transcript USING btree (database);


--
-- Name: database_agg_trancript_id_idx; Type: INDEX; Schema: public; Owner: -
--

CREATE UNIQUE INDEX database_agg_trancript_id_idx ON public.database_agg_transcript USING btree (id);


--
-- Name: depmap_cell_line_mapping_cell_line_idx; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX depmap_cell_line_mapping_cell_line_idx ON public.depmap_cell_line_mapping USING btree (cell_line);


--
-- Name: gene_gene_idx; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX gene_gene_idx ON public.gene USING btree (gene);


--
-- Name: gene_gene_trgm_idx; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX gene_gene_trgm_idx ON public.gene USING gin (gene public.gin_trgm_ops);


--
-- Name: gene_info_symbol_idx; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX gene_info_symbol_idx ON public.gene_info USING btree (symbol);


--
-- Name: unified_term_map_original_term_idx; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX unified_term_map_original_term_idx ON public.unified_term_map USING btree (original_term);


--
-- Name: data data_database_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.data
    ADD CONSTRAINT data_database_fkey FOREIGN KEY (database) REFERENCES public.database(id) ON DELETE CASCADE;


--
-- Name: data data_gene_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.data
    ADD CONSTRAINT data_gene_fkey FOREIGN KEY (gene) REFERENCES public.gene(id) ON DELETE CASCADE;


--
-- Name: data_transcript data_transcript_database_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.data_transcript
    ADD CONSTRAINT data_transcript_database_fkey FOREIGN KEY (database) REFERENCES public.database(id) ON DELETE CASCADE;


--
-- Name: data_transcript data_transcript_transcript_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.data_transcript
    ADD CONSTRAINT data_transcript_transcript_fkey FOREIGN KEY (transcript) REFERENCES public.transcript(id) ON DELETE CASCADE;


--
-- Name: gene_transcript gene_transcript_gene_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.gene_transcript
    ADD CONSTRAINT gene_transcript_gene_fkey FOREIGN KEY (gene) REFERENCES public.gene(id) ON DELETE CASCADE;


--
-- Name: gene_transcript gene_transcript_transcript_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.gene_transcript
    ADD CONSTRAINT gene_transcript_transcript_fkey FOREIGN KEY (transcript) REFERENCES public.transcript(id) ON DELETE CASCADE;


--
-- PostgreSQL database dump complete
--


--
-- Dbmate schema migrations
--

INSERT INTO public.schema_migrations (version) VALUES
    ('20221123142705'),
    ('20221212153155'),
    ('20230117163145'),
    ('20230117191156'),
    ('20230221183030'),
    ('20230223214101'),
    ('20230313174842'),
    ('20230403155038'),
    ('20230406190001'),
    ('20230419203711'),
    ('20230419204219'),
    ('20230727174158'),
    ('20230906154745'),
    ('20230918153613'),
    ('20230920195024'),
    ('20230920201419'),
    ('20230925141013'),
    ('20230925165804'),
    ('20230925181844'),
    ('20231129200027'),
    ('20231205200001'),
    ('20231205220323'),
    ('20231206002401'),
    ('20231207235614'),
    ('20240202144258'),
    ('20240202163120'),
    ('20240222152523'),
    ('20240403141754'),
    ('20240403160325'),
    ('20240417143956'),
    ('20240418212539'),
    ('20240419161005'),
    ('20240610183415');
