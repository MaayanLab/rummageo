#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import scanpy as sc
import glasbey
from adjustText import adjust_text
import os

os.makedirs('figures/fig4', exist_ok=True)
#%%
mouse_clusters = {15: 'Metabolism', 
                  33: 'Ion Channel Transport', 
                  2: 'Spliceosome', 
                  37: 'Phototransduction', 
                  7: 'Sex Determination',
                  14: 'Cilium Assembly',
                  17: 'Interferon Signaling Pathway',
                  25: 'Striated Muscle Contraction',
                  24: 'Cytokine Signaling Pathway', 
                  6: 'Humoral Immune Response',
                  9: 'Synaptic Transmission',
                  31: 'Cilliopathies',
                  30: 'Androgen Biosynthesis',
                  40: 'Cytoplasmic Translation',
                  27: 'Lipolysis Regulation',
                  19: 'Th17 Cell Differentiation',
                  11: 'Melanin Biosynthesis',
                  13: 'Innate Immune System',
                  21: 'Meosis',
                  20: 'Female Gonad Development',
                  22: 'Phototransduction',
                  8: 'Extracellular Matrix Organization',
                  32: 'Cell Cylce',
                  50: 'tRNA Aminoacylation',
                  41: 'Cellular Respiration',
                  65: 'Histone Demethylation',
}

human_clusters = {
   13: 'Extracellular Matrix Organization',
   16: 'Antigen Receptor Signaling',
   6: 'Autophagy Regulation',
   7: 'Carbon Dioxide Transport',
   5: 'Neuronal System',
   1: 'Organelle Assembly',
   23: 'Immune Cytokine Signaling',
   41: 'tRNA Aminoacylation',
   26: 'Triglyceride Biosynthesis',
   19: 'Inflammatory Response',
   9: 'Olfactory Transduction',
   40: 'Melanin Biosynthesis',
   12: 'Glutathione Metabolism',
   3: 'Olfactory transduction',
   24: 'Taste Transduction',
   42: 'Cholesterol Biosynthesis',
   38: 'Fanconi Anemia Pathway',
   35: 'p53 Signaling',
   51: 'Glycolysis',
   32: 'Ciliopathies',
   33: 'Sex Differentiation',
   27: 'Taste Transduction',
   46: 'Histone Demethylation',
   39: 'Cytoplasmic Translation',
   63: 'Protein Deubiquitination'
}

calc_clusters = False
for species in ['human', 'mouse']:
  if species == 'human':
      cluster_labels = human_clusters
      letter = 'a'
  else:
      letter = 'b'
      cluster_labels = mouse_clusters
  file = f'data/umap_genes_{species}_leiden.csv'
  vals = pd.read_csv(f'{file}', index_col=0)
  if calc_clusters:
    adata = sc.AnnData(vals.values)
    adata.obs['genes'] = vals.index.values
    # %%
    sc.pp.neighbors(adata, n_neighbors=5)
    sc.tl.umap(adata, spread=1, min_dist=0.1)
    sc.tl.leiden(adata)

    clusters = adata.obs['leiden'].astype(int).to_numpy()
    unique_clusters = np.unique(clusters)
    df = pd.DataFrame(adata.obsm['X_umap'], columns=['UMAP1', 'UMAP2'], index=adata.obs.genes)
    df['leiden'] = clusters

    # Write DataFrame to a CSV file
    df.to_csv(f'data/umap_genes_{species}_leiden.csv')
  else:
    unique_clusters = list(vals['leiden'].unique())
    clusters = list(vals['leiden'].values)
  cmap = dict(zip(unique_clusters, glasbey.create_palette(len(unique_clusters))))
  # Plot the UMAP with the cluster colors

  #%%
  fig = plt.figure(figsize=(16, 10))
  plt.rcParams.update({'font.size': 20})

  texts = []
  for cluster in list(sorted(unique_clusters)):
      idx = np.where(clusters == cluster)  # Get indices of points in this cluster
      if cluster in cluster_labels:
          x_mean = vals.iloc[idx]['UMAP1'].mean()
          y_mean = vals.iloc[idx]['UMAP2'].mean()
          texts.append(plt.text(x_mean, y_mean, cluster_labels[cluster], fontsize=10, fontweight='bold'))
      plt.scatter( vals.iloc[idx]['UMAP1'], vals.iloc[idx]['UMAP2'], s=.1, color=cmap[cluster], rasterized=True, label=f'Cluster {cluster}')
  plt.legend(handles=[
    Line2D([0], [0], marker='o', color='w', label=f"{label}", # ({cluster_labels[label]})" if label in cluster_labels else label,
          markerfacecolor=cmap[label], markersize=10)
    for label in list(sorted(unique_clusters))
    ], bbox_to_anchor=(1.05, 1), loc='upper left', ncol=3, fontsize=14)

  adjust_text(texts, force_text=(1, 3), min_arrow_len=0.1, arrowprops=dict(arrowstyle='->', color='black'))

  
  plt.xlabel('UMAP-1')
  plt.ylabel('UMAP-2')
  plt.yticks([])
  plt.xticks([])
  plt.tight_layout()
  plt.savefig(f'figures/fig4/4{letter}.png', dpi=300)
  plt.savefig(f'figures/fig4/4{letter}.pdf', dpi=300)
  plt.clf()
  

# %%
