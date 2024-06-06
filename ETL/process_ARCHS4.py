# %%
import os
import json
import h5py as h5
import pandas as pd
import numpy as np
from tqdm import tqdm

import nltk
from nltk.corpus import stopwords
nltk.download('stopwords')
from sklearn.cluster import KMeans
from sentence_transformers import SentenceTransformer

sentence_bert_model = SentenceTransformer('all-mpnet-base-v2')

def get_embeddings(sentences):
    return sentence_bert_model.encode(sentences, batch_size=32, show_progress_bar=False)

pd.options.mode.chained_assignment = None

def partition_samples(species: str, version: str, base_path: str = ""):
    if os.path.exists(f'out/gse_groupings_{species}_{version}.json'):
        return
    os.makedirs('out', exist_ok=True)

    file = f'{base_path}{species}_gene_v{version}.h5'
    single_cell_prob_thresh = 0.5

    f = h5.File(file, "r")
    gse_scprob = np.array([
        f["meta"]["samples"]["series_id"], 
        f["meta"]["samples"]["geo_accession"],
        f["meta"]["samples"]["singlecellprobability"],
        f["meta"]["samples"]["title"],
        f["meta"]["samples"]["characteristics_ch1"],
        f["meta"]["samples"]["source_name_ch1"]
    ]).T
    f.close()

    # %%
    samps_df = pd.DataFrame(gse_scprob, columns =['gse', 'gsm', 'scprob', 'title','characteristics_ch1', 'source_name_ch1'])
    samps_df['gse'] = samps_df['gse'].apply(lambda s: s.decode("utf-8"))
    samps_df['gsm'] = samps_df['gsm'].apply(lambda s: s.decode("utf-8"))
    samps_df['title'] = samps_df['title'].apply(lambda s: s.decode("utf-8"))
    samps_df['characteristics_ch1'] = samps_df['characteristics_ch1'].apply(lambda s: s.decode("utf-8"))
    samps_df['source_name_ch1'] = samps_df['source_name_ch1'].apply(lambda s: s.decode("utf-8"))
    samps_df = samps_df[samps_df['scprob'] < single_cell_prob_thresh]
    samps_df.to_csv(f'out/gse_gsm_meta_{species}.csv', index=False)
    # %%
    with open('processed.json', 'r') as fr:
        processed = json.load(fr)

    gses_processed = []
    for v in processed[species]:
        gses_processed.extend(list(processed[species][v]))
    
    if version not in list(processed[species].keys()):
        to_process = list(set(samps_df['gse'].values).difference(gses_processed))
        processed[species][version] = to_process
    
    samps_df = samps_df[samps_df['gse'].isin(to_process)]

    
    valid = set()
    checked = set()
    for i, row in tqdm(samps_df.iterrows(), total=len(samps_df)):
        if row['gse'] in checked:
            continue
        samps = samps_df[samps_df['gse'] == row['gse']]
        n_samps = len(samps)
        if n_samps >= 6 and n_samps < 50:
            valid.add(row['gse'])
        checked.add(row['gse'])
    
    # %%
    words_to_remove = ['experiement', 'experiment', 'patient', 'batch', '1', '2', '3', '4', '5', '6', '7', '8', '9']
    stopwords_plus = set(stopwords.words('english') + (words_to_remove))

    gse_dict = {}

    for gse in tqdm(valid):
        gse_table = samps_df[samps_df['gse'] == gse]
        meta_names = gse_table['title'] + ' _ ' +  gse_table['characteristics_ch1'] + ' _ ' +  gse_table['source_name_ch1']
        data = list(map(lambda s: s.lower(), meta_names.values))
        data_clean = []
        for d in data:
            data_clean.append(' '.join(list(filter(lambda w: w not in stopwords_plus, d.replace(',', ' ').split(' ')))))

        e = get_embeddings(data_clean)
        embedding_df = pd.DataFrame(e, index=gse_table['gsm'].values)

        kmeans = KMeans(n_clusters= (len(embedding_df) // 3), n_init=10).fit(embedding_df.values)
        gse_table['label'] = kmeans.labels_

        gse_table = gse_table[gse_table['label'].map(gse_table['label'].value_counts()) >= 3]
        if len(gse_table) < 6:
            continue
        
        grouped_gse_table = gse_table.groupby('label')
        gse_dict[gse] = {}
        for label in set(gse_table['label'].values):
            gse_dict[gse][str(label)] = list(grouped_gse_table.get_group(label)['gsm'].values)

    with open(f'out/gse_groupings_{species}_{version}.json', 'w') as fw:
        json.dump(gse_dict, fw)
    
    with open('processed.json', 'w') as fw:
        json.dump(processed, fw)
