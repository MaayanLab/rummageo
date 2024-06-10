import os
import json
import shutil
import gzip



def make_downloads():
    os.makedirs('out/downloads', exist_ok=True)
    gmts = os.listdir('out/gmts')
    metas = os.listdir('out/meta')


    for species in ['human', 'mouse']:
        species_gmts = [f for f in gmts if species in f]
        with open(f'out/downloads/{species}-geo-auto.gmt','wb') as wfd:
            for f in species_gmts:
                with open(f'out/gmts/{f}','rb') as fd:
                    shutil.copyfileobj(fd, wfd)

        species_metas = [f for f in metas if species in f]
        meta_dict_combined = {}
        for f in species_metas:
            print(f)
            with open(f'out/meta/{f}','rb') as fd:
                meta_dict_combined = meta_dict_combined | json.load(fd)
        with open(f'out/downloads/{species}-gse-processed-meta.json','w') as wfd:
            json.dump(meta_dict_combined, wfd)


make_downloads() 


