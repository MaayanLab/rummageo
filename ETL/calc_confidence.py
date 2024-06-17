import re
import os
import json
import h5py as h5
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import zscore
import plotly.graph_objects as go
from maayanlab_bioinformatics.normalization import quantile_normalize
from sklearn.metrics import silhouette_score
from tqdm import tqdm
import math
from itertools import combinations

def dist(p1, p2):
    (x1, y1), (x2, y2) = p1, p2
    return math.sqrt((x2 - x1)**2 + (y2 - y1)**2)

def calc_avg_distance(x, y):
    points = list(zip(x,y))
    distances = [dist(p1, p2) for p1, p2 in combinations(points, 2)]
    avg_distance = sum(distances) / len(distances)
    return avg_distance


def log2_normalize(x, offset=1.):
    return np.log2(x + offset)


def compute_confidence(species: str, version: str, base_path: str = ""):
    if os.path.exists(f'out/meta/gse_processed_meta_{species}_{version}_conf.json'):
        return
    
    with open(f'out/meta/gse_processed_meta_{species}_{version}.json') as f:
        gse_processed_meta = json.load(f)

    f = h5.File(f"{base_path}{species}_gene_v{version}.h5", "r")

    genes = [x.decode('UTF-8') for x in f['meta/genes/symbol']]
    samples = [x.decode('UTF-8') for x in f['meta/samples/geo_accession']] 

    for gse in tqdm(list(gse_processed_meta)):
        gsms = []
        [gsms.extend(gr) for gr in gse_processed_meta[gse]['samples'].values()]
        condition_dict = {}
        condition_names_valid = True

        # a few condition names are duplicated, so we need to make sure that the condition names are unique
        # this may indiciate that conditions were not grouped effectively, but we can observe this later

        if len(list(set(gse_processed_meta[gse]['titles'].values()))) !=  len(list(gse_processed_meta[gse]['titles'].values())):
            condition_names_valid = False

        for cond in gse_processed_meta[gse]['samples']:
            for gsm in gse_processed_meta[gse]['samples'][cond]:
                condition_dict[gsm] = cond
            

        samples_idx = sorted([i for i, x in enumerate(samples) if x in gsms])

        expression_data = f['data/expression'][:, samples_idx]
        #samples_readsaligned = f['meta']['samples']['readsaligned'][samples_idx]
        samples_readsaligned = np.sum(expression_data, axis=0)
        expression_data = np.divide(expression_data, samples_readsaligned)

        expr_df = pd.DataFrame(data=expression_data, index=genes, columns=gsms)
        reference = np.sqrt(expr_df.mean(axis=1))

        for col in expr_df.columns:
            expr_df[col] = np.true_divide(expr_df[col], reference)

        expr_df = expr_df.dropna()
        expr_df

        df_data_norm = log2_normalize(expr_df, offset=1)
        df_data_norm = quantile_normalize(expr_df, axis=0)
        #convert to zscores
        expr_df = pd.DataFrame(zscore(df_data_norm, axis=0), index=df_data_norm.index, columns=df_data_norm.columns)

        X = expr_df.T.values
        pca = PCA(n_components=2)
        result_expr = pca.fit_transform(X)
        pca_df = pd.DataFrame(result_expr, columns = ['PCA-1','PCA-2'])
        pca_df['sample'] = expr_df.columns
        pca_df['condition'] = pca_df['sample'].apply(lambda x: condition_dict[x])

        try:
            s_score = silhouette_score(pca_df[['PCA-1', 'PCA-2']], pca_df['condition'])
            gse_processed_meta[gse]['silhouette_score'] = s_score
        except Exception as e:
            print(gse, e)
            gse_processed_meta[gse]['silhouette_score'] = -2
        
    with open(f'out/meta/gse_processed_meta_{species}_{version}_conf.json', 'w') as f:
        json.dump(gse_processed_meta, f)