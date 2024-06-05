#%%
import pathlib
import numpy as np
import json
import pandas as pd
import matplotlib.pyplot as plt
from common import data_dir, GMT

#%%
fig_dir = pathlib.Path('figures/fig1')
fig_dir.mkdir(parents=True, exist_ok=True)

#%%
species = 'human'
#%%
print(f'Loading Rummageo {species} GMT...')
gmt = GMT.from_file(data_dir/f'{species}-geo-auto.gmt')
df = gmt.to_df()

#%%
print('Computing set sizes...')
set_sizes = pd.Series(
  df.values.sum(axis=1, dtype=np.int64).astype(np.int64),
  index=df.index,
)
# %%
fig = plt.figure(figsize=(3,2))
plt.rcParams['font.size'] = 10
pd.Series(
  df.values.sum(axis=0, dtype=np.int64),
  index=df.columns,
).hist(bins=50, color='black')
plt.grid(False)
plt.ylabel('genes')
plt.xlabel('gene sets containing gene')
plt.tight_layout(pad=2.5)

plt.savefig(fig_dir/'1a.pdf', dpi=300)
plt.savefig(fig_dir/'1a.png', dpi=300)
# plt.show()
plt.clf()
#%%
fig = plt.figure(figsize=(3,2))
set_sizes.hist(bins=50, log='y', color='black')
plt.grid(False)
plt.ylabel('gene sets')
plt.xlabel('gene set length')
plt.tight_layout()
plt.savefig(fig_dir/'1b.pdf', dpi=300)
plt.savefig(fig_dir/'1b.png', dpi=300)
# plt.show()
plt.clf()

# %%
gses = df.index.map(lambda i: i[0].partition('-')[0])
fig = plt.figure(figsize=(3,2))
pd.Series(gses).value_counts().hist(bins=50, log='y', color='black')
plt.grid(False)
plt.ylabel('GSEs')
plt.xlabel('gene sets per GSE')
plt.tight_layout()
plt.savefig(fig_dir/'1c.pdf', dpi=300)
plt.savefig(fig_dir/'1c.png', dpi=300)
#plt.show()
plt.clf()


# %%

plt.rcParams['font.size'] = 10
species = 'mouse'
#%%
print(f'Loading Rummageo {species} GMT...')
gmt = GMT.from_file(data_dir/f'{species}-geo-auto.gmt')
df = gmt.to_df()

#%%
print('Computing set sizes...')
set_sizes = pd.Series(
  df.values.sum(axis=1, dtype=np.int64).astype(np.int64),
  index=df.index,
)
# %%
fig = plt.figure(figsize=(3,2))
pd.Series(
  df.values.sum(axis=0, dtype=np.int64),
  index=df.columns,
).hist(bins=50, color='black')
plt.grid(False)
plt.ylabel('genes')
plt.xlabel('gene sets containing gene')
plt.tight_layout(pad=2.5)
plt.savefig(fig_dir/'1d.pdf', dpi=300)
plt.savefig(fig_dir/'1d.png', dpi=300)
# plt.show()
plt.clf()
# %%
#%%
fig = plt.figure(figsize=(3,2))
set_sizes.hist(bins=50, log='y', color='black')
plt.grid(False)
plt.ylabel('gene sets')
plt.xlabel('gene set length')
plt.tight_layout()
plt.savefig(fig_dir/'1e.pdf', dpi=300)
plt.savefig(fig_dir/'1e.png', dpi=300)
# plt.show()
plt.clf()

# %%
gses = df.index.map(lambda i: i[0].partition('-')[0])
fig = plt.figure(figsize=(3,2))
pd.Series(gses).value_counts().hist(bins=50, log='y', color='black')
plt.grid(False)
plt.ylabel('GSEs')
plt.xlabel('gene sets per GSE')
plt.tight_layout()
plt.savefig(fig_dir/'1f.pdf', dpi=300)
plt.savefig(fig_dir/'1f.png', dpi=300)
plt.show()
plt.clf()
# %%

#%%
with open(data_dir/'gse_processed_meta_human.json') as fr:
    human_gse_attrs = json.load(fr)

num_groups = []
for gse in human_gse_attrs:
    if len(human_gse_attrs[gse]['samples']) > 1:
        num_groups.append(len(human_gse_attrs[gse]['samples']))
plt.rcParams['font.size'] = 10
fig = plt.figure(figsize=(6,4))
pd.Series(num_groups).value_counts().plot.bar(color='black')
plt.grid(False)
plt.ylabel('# GSEs')
plt.xlabel('Groups per Study')
plt.tight_layout()
plt.savefig(fig_dir/'1c-2.pdf', dpi=300)
plt.savefig(fig_dir/'1c-2.png', dpi=300)
plt.show()

# %%
with open(data_dir/'gse_processed_meta_mouse.json') as fr:
    mouse_gse_attrs = json.load(fr)

num_groups = []
for gse in mouse_gse_attrs:
    if len(mouse_gse_attrs[gse]['samples']) > 1:
        num_groups.append(len(mouse_gse_attrs[gse]['samples']))

plt.rcParams['font.size'] = 10
fig = plt.figure(figsize=(6,4))
pd.Series(num_groups).value_counts().plot.bar(color='black')
plt.grid(False)
plt.ylabel('# GSEs')
plt.xlabel('Groups per Study')
plt.tight_layout()
plt.savefig(fig_dir/'1f-2.pdf', dpi=300)
plt.savefig(fig_dir/'1f-2.png', dpi=300)
plt.show()
# %%
# %%

