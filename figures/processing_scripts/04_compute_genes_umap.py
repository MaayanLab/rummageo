#%%
import numpy as np
import pandas as pd
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.decomposition import TruncatedSVD
from umap import UMAP
import matplotlib.pyplot as plt

from common import data_dir, GMT

#%%
random_state = 42
merge=False
#%%
print('Loading Rummagene GMT...')
gmt = GMT.from_file(data_dir/'human-geo-auto.gmt')
if merge:
  gmt2 = GMT.from_file(data_dir/'mouse-geo-auto.gmt')

  from maayanlab_bioinformatics.harmonization import ncbi_genes_lookup
  lookup = ncbi_genes_lookup()

  print('Merging Rummagene GMTs... and mapping protein coding genes')
  gmt.gene_lists = gmt.gene_lists + gmt2.gene_lists
  gmt.terms = gmt.terms + gmt2.terms
  new_bg = set()
  for i in range(len(gmt.gene_lists)):
    gmt.gene_lists[i] = [g.upper() for g in gmt.gene_lists[i] if lookup(g.upper())]
    new_bg.update(gmt.gene_lists[i])
  gmt.background = list(new_bg)
#%%
print('Computing IDF...')
vectorizer = TfidfVectorizer(analyzer=lambda gs: gs)
vectors = vectorizer.fit_transform(gmt.gene_lists)
vectors = vectors.transpose()
umap_index = vectorizer.get_feature_names()
#%%
print('Computing SVD...')
svd = TruncatedSVD(n_components=50, random_state=random_state)
svs = svd.fit_transform(vectors)

pd.DataFrame(svs, index=umap_index).to_csv(data_dir / 'svd_genes_both.tsv', sep='\t')

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
print('Saving umap...')
df_umap = pd.DataFrame(embedding, columns=['UMAP-1', 'UMAP-2'], index=umap_index)
df_umap['outlier'] = (~outlier).astype(int)
df_umap

#%%
df_umap.to_csv(data_dir / 'umap_genes_both.tsv', sep='\t')

plt.scatter(df_umap['UMAP-1'].values, df_umap['UMAP-2'].values, s = .1, alpha=0.1, rasterized=True)
""" plt.legend(handles=[
  Line2D([0], [0], marker='o', color='w', label=label,
        markerfacecolor=cdict[label], markersize=10)
  for label in labels.unique()
]) """
plt.xlabel('UMAP-1')
plt.ylabel('UMAP-2')
plt.tight_layout()
plt.savefig('./figures/genes_umap_both.pdf', dpi=300)
plt.savefig('./figures/genes_umap_both.png', dpi=300)
