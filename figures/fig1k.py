import json
import pandas as pd
from common import data_dir
import pathlib
import matplotlib.pyplot as plt


fig_dir = pathlib.Path('figures/fig1')
fig_dir.mkdir(parents=True, exist_ok=True)

with open(data_dir/'GEO_gse_meta_mouse.json') as fr:
    geo_meta_mouse = json.load(fr)

with open(data_dir/'GEO_gse_meta_human.json') as fr:
    geo_meta_mouse = json.load(fr)

geo_meta = {**geo_meta_mouse, **geo_meta_mouse}

years = []
for gse in geo_meta:
    years.append(geo_meta[gse]['date'].split(' ')[-1])
    
fig = plt.figure(figsize=(6,4))
plt.rcParams['font.size'] = 10
pd.Series(years).value_counts().sort_index().plot.bar(color='black')
plt.grid(False)
plt.ylabel('GSEs in Rummageo')
plt.xlabel('Year Published')
plt.tight_layout()
plt.savefig(fig_dir/'1h.pdf', dpi=300)
plt.savefig(fig_dir/'1h.png', dpi=300)
plt.show()