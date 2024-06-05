
#%%
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib.lines import Line2D
import matplotlib
import os
os.makedirs('figures/fig2', exist_ok=True)
matplotlib.rcParams['font.size'] = 14

umap_df = pd.read_csv('data/umap-human-mouse.tsv', sep='\t', index_col=0)

# %%
labels = umap_df.index.map(lambda t: t.split('-')[-1].split("'")[0].replace('.tsv', ''))
cdict = {'human up': '#e6194B', 'human dn': '#4363d8', 'mouse up': '#f58231', 'mouse dn': '#42d4f4'}

#cdict = {'human up': '#f20505', 'human dn': '#f20505', 'mouse up': '#0525f2', 'mouse dn': '#0525f2'}
# %%

fig = plt.figure(figsize=(10,8))
for g in np.unique(labels):
    ix = np.where(labels == g)
    plt.scatter(umap_df['UMAP-1'].values[ix], umap_df['UMAP-2'].values[ix], c = cdict[g], label =g, s = .01, alpha=.8, rasterized=True)
plt.legend(handles=[
  Line2D([0], [0], marker='o', color='w', label=label,
        markerfacecolor=cdict[label], markersize=10)
  for label in labels.unique()
])

plt.xlabel('UMAP-1')
plt.ylabel('UMAP-2')
plt.tight_layout()
plt.ylim(-3, 15)
plt.xlim(-1, 18)
plt.xticks([])
plt.yticks([])
plt.savefig('figures/fig2/2a.pdf', dpi=300)
plt.savefig('figures/fig2/2a.png', dpi=300)
plt.show()

# %%
