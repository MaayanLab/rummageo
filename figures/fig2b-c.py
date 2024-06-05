#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib
import glasbey
import json
import ast

matplotlib.rcParams['font.size'] = 14
#%%
top_terms_num = 10
for category in ['tissue', 'disease']:
    if category == 'tissue':
        letter = 'b'
    elif category == 'disease':
        letter = 'c'

    with open(f'data/sig_{category}_mapper.json', 'r') as f:
        mapper = json.load(f)

    umap_df = pd.read_csv('data/umap-full.tsv', sep='\t', index_col=0)
    #umap_df = umap_df[umap_df['outlier'] == 0]

    labels = umap_df.index.map(lambda t: ast.literal_eval(t)[0].replace(' ', '-').replace('.tsv', ''))
    mapped_labels = labels.map(lambda t: mapper[t] if t in mapper else "other")

    top10 = mapped_labels.value_counts()[1:top_terms_num].index.values

    labels = mapped_labels.map(lambda t: t if t in top10 else 'other')

    unique_labels = np.unique(labels)
    # %%
    cdict = dict(zip(unique_labels, glasbey.create_palette(len(unique_labels))))
    cdict["other"] = "#b5b5b5"
    # %%
    plt.figure(figsize=(10, 8))
    for g in unique_labels:
        ix = np.where(labels == g)
        if g == "other":
            alpha = 0.1
            s = .1
        else:
            alpha = 1
            s = 1
        plt.scatter(umap_df['UMAP-1'].values[ix], umap_df['UMAP-2'].values[ix], c = cdict[g], label=g, s=s ,alpha=alpha, rasterized=True)
    plt.legend(handles=[
    Line2D([0], [0], marker='o', color='w', label=label,
            markerfacecolor=cdict[label], markersize=10)
    for label in unique_labels
    ], bbox_to_anchor =(0,0), loc='lower left', ncol=2)
    plt.xlabel('UMAP-1')
    plt.ylabel('UMAP-2')
    plt.tight_layout()
    plt.xticks([])
    plt.yticks([])
    plt.xlim(0, 18)
    plt.ylim(-3, 14)
    plt.savefig(f'figures/fig2/2{letter}.pdf', dpi=300)
    plt.savefig(f'figures/fig2/2{letter}.png', dpi=300)
    plt.show()
    plt.clf()
