-- migrate:up

-- Create the materialized view gene_set_pmid

drop materialized view if exists app_public_v2.gene_set_pmid cascade;

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

create table app_public_v2.gse_terms (
  gse varchar not null,
  species varchar not null,
  llm_attrs varchar[],
  pubmed_attrs varchar[],
  mesh_attrs varchar[],
  unique (gse, species)
);

-- migrate:down

drop materialized view if exists app_public_v2.gene_set_pmid cascade;