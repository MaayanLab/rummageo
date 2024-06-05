#%%
import json
import pandas as pd
import matplotlib.pyplot as plt
from common import data_dir
import os
species = 'mouse'
#%%
with open(data_dir/f'gse_attrs_clean_{species}_v4.json') as fr:
    gse_attrs = json.load(fr)

# %%
attr_df = pd.DataFrame.from_records(gse_attrs).T
attr_df
# %%
num_terms_genes = len(attr_df[attr_df['genes'].str.len() != 0])
num_terms_drugs = len(attr_df[attr_df['drug'].str.len() != 0])
num_terms_tissues = len(attr_df[attr_df['tissue'].str.len() != 0])
num_terms_cell_lines = len(attr_df[attr_df['cell_line'].str.len() != 0])
num_terms_diseases = len(attr_df[attr_df['disease'].str.len() != 0])

unique_genes = set()
unique_tissues = set()
unique_cell_lines = set()
unique_drugs = set()
unique_diseases = set()

for gse in gse_attrs:
    unique_genes.update(gse_attrs[gse]['genes'])
    unique_tissues.update(gse_attrs[gse]['tissue'])
    unique_cell_lines.update(gse_attrs[gse]['cell_line'])
    unique_drugs.update(gse_attrs[gse]['drug'])
    unique_diseases.update(gse_attrs[gse]['disease'])
# %%
fig, ax = plt.subplots()
plot_data = [['GSEs containing tissue(s)/cell type(s)', num_terms_tissues], ['GSEs containing drug(s)', num_terms_drugs], ['GSEs mentioning gene(s)', num_terms_genes], ['GSEs containing disease(s)', num_terms_diseases], ['GSEs containing cell lines(s)', num_terms_cell_lines],
             ['Unique gene symbols mentioned', len(unique_genes)], ['Unique tissue(s)/cell type(s)', len(unique_tissues)], ['Unique cell line(s)', len(unique_cell_lines)], ['Unique drug(s)', len(unique_drugs)], ['Unique disease(s)', len(unique_diseases)]]
plot_data = list(sorted(plot_data, key=lambda x: x[1], reverse=True))
bars = ax.barh([p[0] for p in plot_data], [p[1] for p in plot_data], color='black')
ax.bar_label(bars, fmt='{:,.0f}')
ax.set_xlim(0, max([p[1] for p in plot_data]) + 2000)

os.makedirs('figures/fig2', exist_ok=True)
plt.tight_layout()
plt.savefig(f'figures/fig2/fig2{species}.pdf')
plt.savefig(f'figures/fig2/fig2{species}.png')
plt.clf()
# %%


try:
    with open(data_dir/'tf-kinases.json') as fr3:
        tf_kinases = json.load(fr3)
except FileNotFoundError:
    raise FileNotFoundError('Please run make create_gmt to download relatvent files')

tf_set = set(tf_kinases['tfs'])
tf_signatures = []
for gse in gse_attrs:
    gs_mapped = list(map(lambda x: x.upper(), gse_attrs[gse]['genes']))
    tf_intersect = list(tf_set.intersection(gs_mapped))
    if len(tf_intersect) > 0:
        tf_signatures.append(tf_intersect[0])
        continue


kinase_set = set(tf_kinases['kinases'])
kinase_signatures = []
for gse in gse_attrs:
    gs_mapped = list(map(lambda x: x.upper(), gse_attrs[gse]['genes']))
    kinase_intersect = list(kinase_set.intersection(gs_mapped))
    if len(kinase_intersect) > 0:
            
        kinase_signatures.append(kinase_intersect[0])
        continue


# %%
        
fig, ax = plt.subplots()
plot_data = [['GSEs containing TF(s)', len(tf_signatures)], ['GSEs containing kinases(s)', len(kinase_signatures)], ['GSEs containing unique kinases(s)', len(set(kinase_signatures))], ['GSEs containing unique TF(s)', len(set(tf_signatures))]]
plot_data = list(sorted(plot_data, key=lambda x: x[1], reverse=True))
bars = ax.barh([p[0] for p in plot_data], [p[1] for p in plot_data], color='black')
ax.bar_label(bars, fmt='{:,.0f}')
ax.set_xlim(0, max([p[1] for p in plot_data]) + 2000)

os.makedirs('figures/fig2', exist_ok=True)
plt.tight_layout()
plt.savefig(f'figures/fig2/fig2B{species}.pdf')
plt.savefig(f'figures/fig2/fig2B{species}.png')
plt.clf()

# %%
