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

-- Add a comment for the foreign key
comment on materialized view app_public_v2.gene_set_pmid is E'@foreignKey (id) references app_public_v2.gene_set (id)';

-- Create unique and regular indexes
create index gene_set_pmid_id_idx on app_public_v2.gene_set_pmid (id);
create index gene_set_pmid_pmid_idx on app_public_v2.gene_set_pmid (pmid);

-- Grant permissions
grant select on app_public_v2.gene_set_pmid to guest;
grant all privileges on app_public_v2.gene_set_pmid to authenticated;
-- migrate:down

