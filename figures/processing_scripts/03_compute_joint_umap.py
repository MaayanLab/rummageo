#%%
import json
import numpy as np
import pandas as pd
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.decomposition import TruncatedSVD
from umap import UMAP

from common import data_dir, GMT

#%%
random_state = 42

#%%
print('Loading Enrichr GMT...')
enrichr_gmt = GMT.from_file(data_dir/'enrichr-clean.gmt')

#%%
print('Loading RummaGEO GMT...')
rummagene_gmt = GMT.from_file(data_dir/'human-geo-auto.gmt')
gmt2 = GMT.from_file(data_dir/'mouse-geo-auto.gmt')

#%%
from maayanlab_bioinformatics.harmonization import ncbi_genes_lookup

lookup = ncbi_genes_lookup()

print('Merging RummaGEO GMTs... and mapping protein coding genes')
rummagene_gmt.gene_lists = rummagene_gmt.gene_lists + gmt2.gene_lists
rummagene_gmt.terms = rummagene_gmt.terms + gmt2.terms
for i in range(len(rummagene_gmt.gene_lists)):
  rummagene_gmt.gene_lists[i] = [g.upper() for g in rummagene_gmt.gene_lists[i] if lookup(g.upper())]


#%%
print('Collecting metadata...')
meta = pd.DataFrame(
  [
    { 'source': library, 'term': term }
    for library, term in enrichr_gmt.terms
  ] + [
    { 'source': 'rummagene', 'term': term }
    for term, _ in rummagene_gmt.terms
  ]
)

#%%
print('Categorizing terms...')
with (data_dir/'Enrichr'/'datasetStatistics.json').open('r') as fr:
  datasetStatistics = json.load(fr)

categories = {cat['categoryId']: cat['name'] for cat in datasetStatistics['categories']}
library_categories = {lib['libraryName']: categories[lib['categoryId']] for lib in datasetStatistics['statistics']}
library_categories['rummagene'] = 'Rummagene'
meta['category'] = meta['source'].apply(library_categories.get)

#%%
print('Computing IDF...')
vectorizer = TfidfVectorizer(analyzer=lambda gs: gs)
vectors = vectorizer.fit_transform(enrichr_gmt.gene_lists + rummagene_gmt.gene_lists)

#%%
print('Computing SVD...')
svd = TruncatedSVD(n_components=50, random_state=random_state)
svs = svd.fit_transform(vectors)

# %%
print('Computing UMAP...')
umap = UMAP(random_state=random_state)
embedding = umap.fit_transform(svs)

#%%
print('Computing outliers...')
x = embedding[:, 0]
y = embedding[:, 1]
x_min, x_mu, x_std, x_max = np.min(x), np.mean(x), np.std(x), np.max(x)
x_lo, x_hi = max(x_min, x_mu - x_std*1.68), min(x_max, x_mu + x_std*1.68)
y_min, y_mu, y_std, y_max = np.min(y), np.mean(y), np.std(y), np.max(y)
y_lo, y_hi = max(y_min, y_mu - y_std*1.68), min(y_max, y_mu + y_std*1.68)
outlier = (x>=x_lo)&(x<=x_hi)&(y>=y_lo)&(y<=y_hi)

#%%
print('Saving joint-umap...')
meta['UMAP-1'] = x
meta['UMAP-2'] = y
meta['outlier'] = (~outlier).astype(int)
meta.to_csv(data_dir / 'joint-umap_v2.tsv', sep='\t')
