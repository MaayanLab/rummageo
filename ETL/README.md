# RummaGEO ETL

This contains the code for computing new signatures and associated metadata elements for new releases of ARCHS4.
After downloading the new release to the process the pipeline can be run as follows:
`python3 new_release.py <species> <version> [base_path]`

## Details
- `process_ARCHS4.py`: is used to determine the valid Series to process and then to create sample partitions using metadata string embeddings.
- `compute_signatures.py` attempts to identify normal conditions and compute signatures pairwise for each study with partitioned samples.
- `create_meta_dict.py` creates a JSON file with automatically computed condition titles and determines if those titles are valid.
- `calc_confidence.py` compares the metadata clustering and normalized data clustering to compute silhouette scores for each processed Series.
- `enrichr_tags.py` precomputes significant Enrichr terms from selected libraries for each signature.
- `new_realse.py` runs all the above scripts in order to generate the files necessary to migrate the database to latest version of ARCHS4 as specified by the arguments described above.
- `helper.py` is used to ingest the output of the pipeline including the signatures and associated metadata into the RummaGEO database
- `plpy.py` is just a helper for accessing and querying the database in a way similar to how it is done from the database itself