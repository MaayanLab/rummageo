# CONFIG
PYTHON=python3

.PHONY: pg-shell
pg-shell: check-deps ensure-db
	psql $(shell dotenv get DATABASE_URL)

.PHONY: ensure-db
ensure-db: check-deps ensure-env
	dbmate wait
	dbmate up

.PHONY: pg-dump
pg-dump: check-deps ensure-db
	mkdir -p ETL/out
	pg_dump -Fc --no-acl --no-owner $(shell dotenv get DATABASE_URL) > ETL/out/db.dump

.PHONY: check-deps
check-deps:
	@which dbmate > /dev/null || (echo "Install dbmate, see https://github.com/amacneil/dbmate" && exit 1)
	@which dotenv > /dev/null || (echo "Install python-dotenv, see https://pypi.org/project/python-dotenv/" && exit 1)
	@which curl > /dev/null || (echo "Install curl, see https://curl.se/" && exit 1)
	@bash -c "[[ \"$$($(PYTHON) -V)\" =~ \"Python 3\" ]] || (echo 'python3 is missing, maybe you need to override PYTHON' && exit 1)"

.PHONY: ensure-env
ensure-env: .env
	@dotenv get POSTGRES_PASSWORD || dotenv -qauto set POSTGRES_PASSWORD $(shell openssl rand -hex 16)
	dotenv -qauto set DATABASE_URL postgres://$(shell dotenv get POSTGRES_USER):$(shell dotenv get POSTGRES_PASSWORD)@$(shell dotenv get POSTGRES_HOST):5432/$(shell dotenv get POSTGRES_DB)?sslmode=disable

data/human-geo-auto.gmt:
	echo "Downloading human-geo-auto.gmt"
	curl -s -o data/human-geo-auto.gmt.gz https://s3.amazonaws.com/maayanlab-public/rummageo/human-geo-auto.gmt.gz
	gunzip data/human-geo-auto.gmt.gz

data/mouse-geo-auto.gmt:
	echo "Downloading mouse-geo-auto.gmt"
	curl -s -o data/mouse-geo-auto.gmt.gz https://s3.amazonaws.com/maayanlab-public/rummageo/mouse-geo-auto.gmt.gz
	gunzip data/mouse-geo-auto.gmt.gz

data/human-gse-processed-meta.json:
	echo "Downloading human-gse-processed-meta.json"
	curl -s -o data/human-gse-processed-meta.json https://s3.amazonaws.com/maayanlab-public/rummageo/gse_processed_meta_human.json

data/mouse-gse-processed-meta.json:
	echo "Downloading mouse-gse-processed-meta.json"
	curl -s -o data/mouse-gse-processed-meta.json https://s3.amazonaws.com/maayanlab-public/rummageo/gse_processed_meta_mouse.json

data/enrichr-terms-mouse.json:
	echo "Downloading enrichr-terms-mouse.json"
	curl -s -o data/enrichr-terms-mouse.json https://minio.dev.maayanlab.cloud/rummageo/enrichr-terms-mouse.json

data/enrichr-terms-human.json:
	echo "Downloading enrichr-terms-human.json"
	curl -s -o data/enrichr-terms-human.json https://minio.dev.maayanlab.cloud/rummageo/enrichr-terms-human.json

data/gse_gsm_meta_human.csv:
	echo "Downloading gse_gsm_meta_human.csv"
	curl -s -o data/gse_gsm_meta_human.csv https://minio.dev.maayanlab.cloud/rummageo/gse_gsm_meta_human.csv

data/gse_gsm_meta_mouse.csv:
	echo "Downloading gse_gsm_meta_mouse.csv"
	curl -s -o data/gse_gsm_meta_mouse.csv https://minio.dev.maayanlab.cloud/rummageo/gse_gsm_meta_mouse.csv

data/keyterms_human.json:
	echo "Downloading keyterms_human_2.4.json"
	curl -s -o data/keyterms_human.json https://minio.dev.maayanlab.cloud/rummageo/keyterms_human_2.4.json

data/keyterms_mouse.json:
	echo "Downloading keyterms_mouse_2.4.json"
	curl -s -o data/keyterms_mouse.json https://minio.dev.maayanlab.cloud/rummageo/keyterms_mouse_2.4.json

data/keyterm_categories.json:
	echo "Downloading keyterm_categories.json"
	curl -s -o data/keyterm_categories.json https://minio.dev.maayanlab.cloud/rummageo/keyterm_categories_2.4.json

.PHONY: ingest-db
ingest-db: data/human-geo-auto.gmt data/mouse-geo-auto.gmt data/human-gse-processed-meta.json data/mouse-gse-processed-meta.json data/enrichr-terms-mouse.json data/enrichr-terms-human.json data/gse_gsm_meta_human.csv data/gse_gsm_meta_mouse.csv data/keyterms_human.json data/keyterms_mouse.json data/keyterm_categories.json
	$(PYTHON) ETL/helper.py ingest -i data/human-geo-auto.gmt --species human
	$(PYTHON) ETL/helper.py ingest -i data/mouse-geo-auto.gmt --species mouse
	$(PYTHON) ETL/helper.py ingest-gse-info --species human
	$(PYTHON) ETL/helper.py ingest-gse-info --species mouse
	$(PYTHON) ETL/helper.py ingest-gse-attrs --species human
	$(PYTHON) ETL/helper.py ingest-gse-attrs --species mouse
	$(PYTHON) ETL/helper.py ingest-term-categories
	$(PYTHON) ETL/helper.py ingest-enrichr-terms --species human
	$(PYTHON) ETL/helper.py ingest-enrichr-terms --species mouse

