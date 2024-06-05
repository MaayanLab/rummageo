#%%
from common import data_dir
import json
import plotly.express as px
import pandas as pd
#%%
with open(data_dir/'gse_processed_meta_mouse_conf.json') as fr:
    mouse_conf = json.load(fr)

with open(data_dir/'gse_processed_meta_human_conf.json') as fr:
    human_conf = json.load(fr)


# %%
mouse_confs = []
for gse in mouse_conf:
    if  mouse_conf[gse]["silhouette_score"] != -2:
        mouse_confs.append(mouse_conf[gse]["silhouette_score"])

human_confs = []
for gse in human_conf:
    if human_conf[gse]["silhouette_score"] != -2:
        human_confs.append(human_conf[gse]["silhouette_score"])

#%%
conf_df = pd.DataFrame(data=[mouse_confs, human_confs]).T
conf_df.columns = ['mouse', 'human']
print(len(mouse_confs), len(mouse_conf))
print(len(human_confs), len(human_conf))
print(min(mouse_confs), min(human_confs))
#%%
import numpy as np
import plotly.graph_objects as go
fig = go.Figure()
for (species, vec) in [['human', human_confs], ['mouse', mouse_confs]]:
        fig.add_trace(
            go.Violin(
                y=vec,
                name=species,
                box_visible=True, 
                line_color='black',
                meanline_visible=True,
                points=False
            )
        )

fig.update_layout(yaxis_title='Silhouette Score', plot_bgcolor='white', yaxis_gridcolor='lightgrey',width=1000, height=1000, showlegend=False, font=dict(size=30), yaxis_range=[-1, 1.1], yaxis = dict(
        tickmode = 'linear',
        tick0 = -1,
        dtick = 0.25
    ))
fig.write_image('figures/silhoutette_score.png')
fig.write_image('figures/silhoutette_score.pdf')
# %%