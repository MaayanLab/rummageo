PYTHON=python3

data/human-geo-auto.gmt: data/human-geo-auto.gmt.gz
	curl -L -o $< https://minio.dev.maayanlab.cloud/rummageo/human_12-04-2023.gmt.gz
	gunzip $<

data/mouse-geo-auto.gmt: data/mouse-geo-auto.gmt.gz
	curl -L -o $< https://minio.dev.maayanlab.cloud/rummageo/mouse_12-04-2023.gmt.gz
	gunzip $<

data/gse_processed_meta_mouse.json:
	curl -L -o $@ https://minio.dev.maayanlab.cloud/rummageo/gse_processed_meta_human.json

data/gse_processed_meta_mouse.json:
	curl -L -o $@ https://minio.dev.maayanlab.cloud/rummageo/gse_processed_meta_mouse.json

data/DisGeNET_terms.json:
	curl -L -o $@ https://minio.dev.maayanlab.cloud/rummageo/DisGeNET_terms.json

data/drugs_parsed.json:
	curl -L -o $@ https://minio.dev.maayanlab.cloud/rummageo/drugs_parsed.json

data/bto.owl:
	curl -L -o $@ http://purl.obolibrary.org/obo/bto.owl
	$(PYTHON) parseBTO.py

data/drugs_autoRIF.tsv.gz:
	curl -L -o $@ https://appyters.maayanlab.cloud/storage/DrugShot/DrugRIF.tsv.gz

data/tf-kinase.json:
	curl -L -o $@ https://minio.dev.maayanlab.cloud/table-mining/tf-kinases.json

data/tf-kinase.json:
	curl -L -o $@ https://minio.dev.maayanlab.cloud/table-mining/tf-kinases.json

data/drug_synonyms.json:
	curl -L -o $@ https://minio.dev.maayanlab.cloud/rummageo/compounds-to-filter.txt
	$(PYTHON) get_drug_synonyms.py

data/gse_attrs_clean_human.json:
	$(PYTHON) identify_funtional_terms.py

data/gse_attrs_clean_mouse.json:
	$(PYTHON) identify_funtional_terms.py

figures/fig1: data/human-geo-auto.gmt data/mouse-geo-auto.gmt
	$(PYTHON) figure_scripts/fig1.py

figures/fig2: data/gse_processed_meta_human.json data/gse_processed_meta_mouse.json data/bto.owl
	$(PYTHON) figure_scripts/fig2.py

.PHONY: export_figs
export_figs: figures/fig1


