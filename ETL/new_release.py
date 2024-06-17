import sys
from process_ARCHS4 import *
from compute_signatures import *
from calc_confidence import *
from create_meta_dict import *
from enrichr_tags import *
from create_gmt import *
from make_downloads import *

def new_release(species: str, version: str, base_path: str = ""):
    print("Partitioning new samples for", species, version)
    partition_samples(species, version, base_path)
    print("Computing signatures for", species, version)
    run_compute_sigs(species, version, base_path)
    print("Creating meta dict for", species, version)
    create_meta_dict(species, version, base_path)
    print("Computing confidence for", species, version)
    compute_confidence(species, version, base_path)
    print("Creating GMT for", species, version)
    create_gmt(species, version)
    print("Computing Enrichr tags for", species, version)
    compute_enrichr_labels(species, version)
    print("Creating updated download files")
    make_downloads(species, version)



if __name__ == '__main__':
    if len(sys.argv) != 3 and len(sys.argv) != 4:
        exit("Usage: python new_release.py <species> <version> [base_path]")
    species = sys.argv[1]
    version = sys.argv[2]
    base_path = sys.argv[3] if len(sys.argv) == 4 else ""
    new_release(species, version, base_path)
    

    