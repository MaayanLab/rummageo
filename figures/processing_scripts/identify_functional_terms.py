#%%
import re
import os
import json
import pathlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
from common import data_dir, GMT

#%%
species = 'mouse'
version = '2.4'

# %%
with open(data_dir/f'gse_processed_meta_{species}_{version}_conf.json') as f:
    conditions_meta = json.load(f)
# %%
meta_df = pd.read_csv(data_dir/f'gse_gsm_meta_{species}.csv', index_col=0)
# %%
meta_df = meta_df[meta_df.index.isin(list(conditions_meta.keys()))]
# %%
with open(data_dir/'BTO_cell_lines.json') as f:
    cell_lines = json.load(f)

with open(data_dir/'BTO_tissue_cell_types.json') as f2:
    tissues_cell_types = json.load(f2)

# parsed from FDA download: https://download.open.fda.gov/drug/drugsfda/drug-drugsfda-0001-of-0001.json.zip
with open(data_dir/'drugs_parsed.json') as f3:
    drugs = json.load(f3)

with open(data_dir/'DisGeNET_terms.json') as f4:
    diseases = json.load(f4)
    
with open(data_dir/'drug_synonyms.json') as f5:
    drug_synonyms = json.load(f5)

drugs_to_search = list(drug_synonyms.values())
# %%
gse_attrs = {}

# %%
for gse, gsm_attrs in tqdm(meta_df.iterrows(), total=len(meta_df), desc='Identifying drugs tissues, and cell types, and cell lines from GSM'):
    if gse not in gse_attrs:
        gse_attrs[gse] = {
            'disease': [], 
            'tissue': [], 
            'cell_line': [],
            'drug': []
        }
    combined_gsm_desc = f"{gsm_attrs['title']} {gsm_attrs['characteristics_ch1']} {gsm_attrs['source_name_ch1']}"
    gse_attrs[gse]['tissue'].extend(re.findall(r'\b(?:(?=\S)' + '|'.join(map(re.escape, tissues_cell_types)) + r')\b', combined_gsm_desc))
    gse_attrs[gse]['cell_line'].extend(re.findall(r'\b(?:(?=\S)' + '|'.join(map(re.escape, cell_lines)) + r')\b', combined_gsm_desc))
    gse_attrs[gse]['drug'].extend(re.findall(r'\b(?:(?=\S)' + '|'.join(map(re.escape, drugs_to_search)) + r')\b', combined_gsm_desc))
    gse_attrs[gse]['disease'].extend(re.findall(r'\b(?:(?=\S)' + '|'.join(map(re.escape, diseases)) + r')\b', combined_gsm_desc))

#%%
#deduplicate
for gse in gse_attrs:
    gse_attrs[gse]['tissue'] = list(set(gse_attrs[gse]['tissue']))
    gse_attrs[gse]['cell_line'] = list(set(gse_attrs[gse]['cell_line']))
    gse_attrs[gse]['drug'] = list(set(gse_attrs[gse]['drug']))
    gse_attrs[gse]['disease'] = list(set(gse_attrs[gse]['disease']))

# %%
from maayanlab_bioinformatics.harmonization import ncbi_genes_lookup

def is_valid(g):
    """Remove terms that are just numbers -- these are likely 
    dates/ids in the table name and not gene ids"""
    if g == 'TF' or len(g) <= 2:
        return False
    try:
        return str(int(g)) != g
    except:
        return True

if species == 'human':
    lookup = ncbi_genes_lookup(organism='Mammalia/Homo_sapiens', filters=lambda ncbi: ncbi['type_of_gene'].isin(['protein-coding','ncRNA']))
else:
    lookup = ncbi_genes_lookup(organism='Mammalia/Mus_musculus', filters=lambda ncbi: ncbi['type_of_gene'].isin(['protein-coding','ncRNA']))
#%%
    
pattern = r'[-_.,:;]'
for gse, gsm_attrs in tqdm(meta_df.iterrows(), total=len(meta_df), desc='Identifying gene symbols/synonyms in sample info'):
    combined_gsm_desc = f"{gsm_attrs['title']} {gsm_attrs['characteristics_ch1']} {gsm_attrs['source_name_ch1']}"
    combined_gsm_desc_list = re.sub(pattern, ' ', combined_gsm_desc).split(' ')
    genes = [lookup(g) for g in combined_gsm_desc_list if lookup(g) and is_valid(g)]
    gse_attrs[gse]['genes'] = [g for g in genes if g]

#%%
for gse in gse_attrs:
    gse_attrs[gse]['genes'] = list(set(gse_attrs[gse]['genes']))

#%%
# load in series meta:
with open(f'data/GEO_gse_meta_{species}_{version}.json') as fr:
    gse_meta = json.load(fr)

# %%
pattern = r'[-_.,:;]'
# add gse (study) level metadata to attributes
for gse in tqdm(gse_meta, total=len(gse_meta), desc='Identifying drugs tissues, and cell types, and cell lines from GSE metadata'):
    if gse not in gse_attrs:
        gse_attrs[gse] = {
            'disease': [], 
            'tissue': [], 
            'cell_line': [],
            'drug': []
        }
    combined_gse_desc = f"{gse_meta[gse]['title']} {gse_meta[gse]['summary']} {gse_meta[gse]['overall_design']}"
    gse_attrs[gse]['tissue'].extend(re.findall(r'\b(?:(?=\S)' + '|'.join(map(re.escape, tissues_cell_types)) + r')\b', combined_gse_desc))
    gse_attrs[gse]['cell_line'].extend(re.findall(r'\b(?:(?=\S)' + '|'.join(map(re.escape, cell_lines)) + r')\b', combined_gse_desc))
    gse_attrs[gse]['drug'].extend(re.findall(r'\b(?:(?=\S)' + '|'.join(map(re.escape, drugs_to_search)) + r')\b', combined_gse_desc))
    gse_attrs[gse]['disease'].extend(re.findall(r'\b(?:(?=\S)' + '|'.join(map(re.escape, diseases)) + r')\b', combined_gse_desc))
    combined_gse_desc_list = re.sub(pattern, ' ', combined_gse_desc).split(' ')
    genes = [lookup(g) for g in combined_gse_desc_list if lookup(g) and is_valid(g)]
    gse_attrs[gse]['genes'].extend([g for g in genes if g])

#deduplicate
for gse in gse_attrs:
    gse_attrs[gse]['tissue'] = list(set(gse_attrs[gse]['tissue']))
    gse_attrs[gse]['cell_line'] = list(set(gse_attrs[gse]['cell_line']))
    gse_attrs[gse]['drug'] = list(set(gse_attrs[gse]['drug']))
    gse_attrs[gse]['disease'] = list(set(gse_attrs[gse]['disease']))
    gse_attrs[gse]['genes'] = list(set(gse_attrs[gse]['genes']))

#%%
## manually remove some of the BTO tissues/ cell lines that are not categorized correctly or just too general
# also fix some common mappings such as peripheral blood -> blood
tissues_to_remove = ['NA', 'stem', 'spike', 'gland', 'growth medium', 'blast', 'primary cell', 'v1', 'medium', 'arm', 'back', 'adult', 'trabecular meshwork', 'leg', 'needle', 'animal', 'all cell', 'egg', 'milk', 'forearm', 'yolk', 'scale', 'child', 'root', 'orbit', 'molar', 'jaw', 'incisor', 'finger', 'crypt', 'pregnant', 
                     'leg', 'tier', 'limb', 'pulp', 'bud', 'hip', 'knee', 'culture medium', 'pod', 'foot', 'wart', 'hand', 'ear', 'rosette', 'seed', 'bulb', 'disc', 'fat pad', 'back', 'chest', 'cap', 'head', 'hair', 'lh', 'mantle', '']
tissues_to_move = ['293 cell', '697 cell', 's2 cell', 'c6 cell', 'c2 cell', 'c10 cell', 'a9 cell', 'huvec']
tissue_mapping = {"peripheral blood": "blood", "marrow": "bone marrow", "whole blood": "blood", 'human induced pluripotent stem cell': 'hipsc cell'}

cell_lines_to_remove = ['NA', 'culture', 'hematopoietic stem', 'ALL cell line', 'neural crest', 'trophoblast stem', 'human aortic endothelial', 'progenitor cell', 'hematopoietic stem', 'trophoblast stem', 'neck squamous cell carcinoma cell line', 'primary effusion lymphoma', 'GBM cell', 'human hepatoma cell line', 'lymphoblastoid', 'myelogenous leukemia cell line', 'hepatocellular carcinoma cell line', 'primary culture']
cell_lines_to_move = ['lymphoma', 'myeloid leukemia', 'glioblastoma', 'glioma', 'osteosarcoma', 'lymphoblastoid', 'acute lymphoblastic leukemia', 'hepatoma', 'osteosarcoma', 'glioblastoma', 'liposarcoma', 'hepatoblastoma', 'hepatoma', "Ewing's sarcoma", 'oligodendroglioma', 'liposarcoma']
# %%
for gse in gse_attrs:
    gse_attrs[gse]['tissue'] = [tissue_mapping[t] if t in tissue_mapping else t for t in gse_attrs[gse]['tissue']]
    gse_attrs[gse]['tissue'] = [t for t in gse_attrs[gse]['tissue'] if t not in tissues_to_remove]
    if len(set(tissues_to_move).intersection(set(gse_attrs[gse]['tissue']))) > 0:
        for t in tissues_to_move:
            if t in gse_attrs[gse]['tissue']:
                gse_attrs[gse]['tissue'].remove(t)
                gse_attrs[gse]['cell_line'].append(t)
    gse_attrs[gse]['tissue'] = [cl.replace(' tissue', '') for cl in gse_attrs[gse]['tissue']]
    gse_attrs[gse]['tissue'] = list(set(gse_attrs[gse]['tissue']))

    gse_attrs[gse]['cell_line'] = [t for t in gse_attrs[gse]['cell_line'] if t not in cell_lines_to_remove]
    if len(set(tissues_to_move).intersection(set(gse_attrs[gse]['tissue']))) > 0:
        for cl in cell_lines_to_move:
            if cl in gse_attrs[gse]['cell_line']:
                gse_attrs[gse]['cell_line'].remove(cl)
                gse_attrs[gse]['disease'].append(cl)
    gse_attrs[gse]['cell_line'] = list(set(gse_attrs[gse]['cell_line']))
    gse_attrs[gse]['disease'] = list(set(gse_attrs[gse]['disease']))
    for cl in gse_attrs[gse]['cell_line']:
        if 'glimoa' in cl.lower() or 'sarcoma' in cl.lower() or 'blastoma' in cl.lower() or 'carcinoma' in cl.lower() or 'leukemia' in cl.lower() or 'lymphoma' in cl.lower() or 'melanoma' in cl.lower() or 'neuroblastoma' in cl.lower() or 'neuroendocrine' in cl.lower() or 'neuroepithelioma' in cl.lower() or 'neurofibroma' in cl.lower() or 'neuroglioma' in cl.lower() or 'neuroma' in cl.lower() or 'neuropathy' in cl.lower() or 'neurosis' in cl.lower() or 'neurotropic' in cl.lower() or 'neurotropism' in cl.lower() or 'neurotrophic' in cl.lower() or 'neurotoxic' in cl.lower() or 'neurotoxicity' in cl.lower() or 'neurotransmitter' in cl.lower() or 'cancer' in cl.lower() or 'line' in cl.lower():
            gse_attrs[gse]['cell_line'].remove(cl)
    gse_attrs[gse]['cell_line'] = [cl.replace(' cell', '') for cl in gse_attrs[gse]['cell_line']]
    gse_attrs[gse]['cell_line'] = list(set(gse_attrs[gse]['cell_line']))
        

# %%

with open(data_dir/f'gse_attrs_clean_{species}_v{version}.json', 'w') as fw:
    json.dump(gse_attrs, fw, indent=4)

# %%
import json
import pathlib
data_dir = pathlib.Path('data')
species = 'human'
version = '2.4'

with open(data_dir/f'gse_attrs_clean_{species}_v{version}.json') as fr:
    gse_attrs_2_4 = json.load(fr)

with open(data_dir/f'gse_attrs_clean_{species}_v4.json') as fr2:
    gse_attrs_4 = json.load(fr2)


# %%
gse_attrs_5 = gse_attrs_4 | gse_attrs_2_4

#%%
with open(f'data/gse_attrs_clean_{species}_v5.json', 'w') as fw:
    json.dump(gse_attrs_5, fw, indent=4)

# %%
