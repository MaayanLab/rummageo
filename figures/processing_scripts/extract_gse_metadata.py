import GEOparse
import json
from tqdm import tqdm

species = 'mouse'

def extract_gse_metadata(gse_id):
    gse = GEOparse.get_GEO(gse_id, destdir='data/geo', silent=True)
    print()
    gse_attrs = {
        'title': ' '.join(gse.metadata["title"]),
        'summary': ' '.join(gse.metadata["summary"]),
        'overall_design': ' '.join(gse.metadata["overall_design"]),
        'date': ' '.join(gse.metadata["submission_date"])
    }
    return gse_attrs

version = 2.4
species = 'mouse'

with open(f'data/gse_processed_meta_{species}_{version}_conf.json') as fr:
    gse_attrs = json.load(fr)

gse_meta = {}
for gse in tqdm(list(gse_attrs)):
    try:
        if ',' in gse:
            gse_split = gse.split(',')[0]
            gse_meta[gse] = extract_gse_metadata(gse_split)
        else:
            gse_meta[gse] = extract_gse_metadata(gse)
    except Exception as e:
        print(f"Error processing {gse}: {e}")



with open(f'GEO_gse_meta_{species}_{version}.json', 'w') as fw:
    json.dump(gse_meta, fw, indent=4)
