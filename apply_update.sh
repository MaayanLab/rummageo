#!/bin/bash

PYTHON=python3
VERSION=$1
source venv/bin/activate

if [ -z "$VERSION" ]; then
    echo "Usage: ./apply_update.sh VERSION"
    exit 1
fi

echo "Ingesting new GMT files..."
$PYTHON ETL/helper.py ingest -i ETL/out/gmts/human-geo-auto_$VERSION.gmt --species human
$PYTHON ETL/helper.py ingest -i ETL/out/gmts/mouse-geo-auto_$VERSION.gmt --species mouse
echo "Ingesting GSE metadata..."
$PYTHON ETL/helper.py ingest-gse-info --species human --path ETL/out/gse_info_human_$VERSION.json
$PYTHON ETL/helper.py ingest-gse-info --species mouse --path ETL/out/gse_info_mouse_$VERSION.json
echo "Ingesting GSE key terms..."
$PYTHON ETL/helper.py ingest-gse-attrs --species human --path ETL/out/keyterms/gse_key_terms_clean_human_$VERSION.json
$PYTHON ETL/helper.py ingest-gse-attrs --species mouse --path ETL/out/keyterms/gse_key_terms_clean_mouse_$VERSION.json
echo "Ingesting PubMed metadata..."
$PYTHON ETL/helper.py ingest-pb-info
echo "Ingesting term categories..."
$PYTHON ETL/helper.py ingest-term-categories --path ETL/out/keyterms/key_terms_categorized_human_$VERSION.json
$PYTHON ETL/helper.py ingest-term-categories --path ETL/out/keyterms/key_terms_categorized_mouse_$VERSION.json
echo "Ingesting Enrichr terms..."
$PYTHON ETL/helper.py ingest-enrichr-terms --species human --path ETL/out/enrichr/enrichr_terms_human_$VERSION.json
$PYTHON ETL/helper.py ingest-enrichr-terms --species mouse --path ETL/out/enrichr/enrichr_terms_mouse_$VERSION.json
echo "Updating backgrounds..."
$PYTHON ETL/helper.py update-background --species human
$PYTHON ETL/helper.py update-background --species mouse

