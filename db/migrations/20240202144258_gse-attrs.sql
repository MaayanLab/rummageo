-- migrate:up

alter table app_public_v2.gse_info
add column gse_attrs varchar[],
add column silhouette_score float;


-- migrate:down

alter table app_public_v2.gse_info
drop column gse_attrs,
drop column silhouette_score;
