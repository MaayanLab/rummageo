-- migrate:up

ALTER TABLE app_public_v2.gene_set ADD COLUMN created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP;

-- Update existing records to set the created_at column to December 8, 2023 using explicit TIMESTAMP format
UPDATE app_public_v2.gene_set
SET created_at = TO_TIMESTAMP('2023-12-08 00:00:00', 'YYYY-MM-DD HH24:MI:SS')
WHERE created_at IS NULL;

-- migrate:down

drop column app_public_v2.gene_set.created_at;
