import os
import json
import shutil
import gzip
from tqdm import tqdm



def make_downloads():
    os.makedirs('out/downloads', exist_ok=True)
    gmts = os.listdir('out/gmts')
    metas = os.listdir('out/meta')
    enrichr_tags = os.listdir('out/enrichr')


    for species in ['human', 'mouse']:
        species_gmts = [f for f in gmts if species in f]
        sigs = set()
        with open(f'out/downloads/{species}-geo-auto.gmt','w') as wfd:
            for f in species_gmts:
                with open(f'out/gmts/{f}','r') as fd:
                    lines = fd.readlines()
                    for l in tqdm(lines):
                        l_split = l.split('\t')
                        term = l_split[0]
                        if term in sigs:
                            continue
                        sigs.add(term)
                        genes = '\t'.join(l_split[2:])
                        wfd.write(f"{term}\t\t{genes}")

                    

        species_metas = [f for f in metas if species in f]
        meta_dict_combined = {}
        for f in species_metas:
            with open(f'out/meta/{f}','rb') as fd:
                meta_dict_combined = meta_dict_combined | json.load(fd)
        with open(f'out/downloads/{species}-gse-processed-meta.json','w') as wfd:
            json.dump(meta_dict_combined, wfd)

        enrichr_tags_metas = [f for f in enrichr_tags if species in f]
        enrichr_tags_dict_combined = []
        for f in enrichr_tags_metas:
            with open(f'out/enrichr/{f}','rb') as fd:
                enrichr_tags_dict_combined.extend(json.load(fd))
        with open(f'out/downloads/enrichr-terms-{species}.json','w') as wfd:
            json.dump(enrichr_tags_dict_combined, wfd)