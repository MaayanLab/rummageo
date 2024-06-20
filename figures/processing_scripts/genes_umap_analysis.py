#%%
import pandas as pd
import numpy as np
import numpy as np
import plotly.graph_objects as go
from scipy import stats
import numpy as np
import requests
import json
from time import sleep
from tqdm import tqdm
from common import GMT

def submit_to_enrichr(genes, gene_set_libraries, description=None):
    ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/addList'
    genes_str = '\n'.join(genes)
    payload = {
        'list': (None, genes_str),
        'description': (None, description)
    }

    response = requests.post(ENRICHR_URL, files=payload)
    if not response.ok:
        raise Exception('Error analyzing gene list')

    data = json.loads(response.text)
    user_list_id = data['userListId']
    link = f"https://maayanlab.cloud/Enrichr/enrich?dataset={data['shortId']}"

    ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/enrich'
    query_string = '?userListId=%s&backgroundType=%s'

    enrich_results = []
    for gene_set_library in gene_set_libraries:
        response = requests.get(
            ENRICHR_URL + query_string % (user_list_id, gene_set_library)
        )

        if not response.ok:
            print(response.text)
            print(response.status_code)
            raise Exception('Error fetching enrichment results')

        data = json.loads(response.text)
        if gene_set_library not in data:
            enrich_results.append([[1, 'NA']])
        else:
            top3_terms = data[gene_set_library][:3]
            if len(top3_terms) == 0:
                enrich_results.append([[1, 'NA']])
            else:
                enrich_results.append(top3_terms)
        sleep(.2)
    return link, enrich_results



def add_p_value_annotation(fig, array_columns, subplot=None, _format=dict(interline=0.07, text_height=1.07, color='black')):
    ''' Adds notations giving the p-value between two box plot data (t-test two-sided comparison)
    
    Parameters:
    ----------
    fig: figure
        plotly boxplot figure
    array_columns: np.array
        array of which columns to compare 
        e.g.: [[0,1], [1,2]] compares column 0 with 1 and 1 with 2
    subplot: None or int
        specifies if the figures has subplots and what subplot to add the notation to
    _format: dict
        format characteristics for the lines

    Returns:
    -------
    fig: figure
        figure with the added notation
    '''
    # Specify in what y_range to plot for each pair of columns
    y_range = np.zeros([len(array_columns), 2])
    for i in range(len(array_columns)):
        y_range[i] = [1.01+i*_format['interline'], 1.02+i*_format['interline']]

    # Get values from figure
    fig_dict = fig.to_dict()

    # Get indices if working with subplots
    if subplot:
        if subplot == 1:
            subplot_str = ''
        else:
            subplot_str =str(subplot)
        indices = [] #Change the box index to the indices of the data for that subplot
        for index, data in enumerate(fig_dict['data']):
            #print(index, data['xaxis'], 'x' + subplot_str)
            if data['xaxis'] == 'x' + subplot_str:
                indices = np.append(indices, index)
        indices = [int(i) for i in indices]
        print((indices))
    else:
        subplot_str = ''

    # Print the p-values
    for index, column_pair in enumerate(array_columns):
        if subplot:
            data_pair = [indices[column_pair[0]], indices[column_pair[1]]]
        else:
            data_pair = column_pair

        # Mare sure it is selecting the data and subplot you want
        #print('0:', fig_dict['data'][data_pair[0]]['name'], fig_dict['data'][data_pair[0]]['xaxis'])
        #print('1:', fig_dict['data'][data_pair[1]]['name'], fig_dict['data'][data_pair[1]]['xaxis'])

        # Get the p-value
        pvalue = stats.ttest_ind(
            fig_dict['data'][data_pair[0]]['y'],
            fig_dict['data'][data_pair[1]]['y'],
            equal_var=False,
        )[1]
        print(pvalue)
        if pvalue >= 0.05:
            symbol = 'ns'
        elif pvalue >= 0.01: 
            symbol = '*'
        elif pvalue >= 0.001:
            symbol = '**'
        else:
            symbol = '***'
        # Vertical line
        fig.add_shape(type="line",
            xref="x"+subplot_str, yref="y"+subplot_str+" domain",
            x0=column_pair[0], y0=y_range[index][0], 
            x1=column_pair[0], y1=y_range[index][1],
            line=dict(color=_format['color'], width=2,)
        )
        # Horizontal line
        fig.add_shape(type="line",
            xref="x"+subplot_str, yref="y"+subplot_str+" domain",
            x0=column_pair[0], y0=y_range[index][1], 
            x1=column_pair[1], y1=y_range[index][1],
            line=dict(color=_format['color'], width=2,)
        )
        # Vertical line
        fig.add_shape(type="line",
            xref="x"+subplot_str, yref="y"+subplot_str+" domain",
            x0=column_pair[1], y0=y_range[index][0], 
            x1=column_pair[1], y1=y_range[index][1],
            line=dict(color=_format['color'], width=2,)
        )
        ## add text at the correct x, y coordinates
        ## for bars, there is a direct mapping from the bar number to 0, 1, 2...
        fig.add_annotation(dict(font=dict(color=_format['color'],size=14),
            x=(column_pair[0] + column_pair[1])/2,
            y=y_range[index][1]*_format['text_height'],
            showarrow=False,
            text=symbol,
            textangle=0,
            xref="x"+subplot_str,
            yref="y"+subplot_str+" domain"
        ))
    return fig


def species_cluster_chr_perc(species, random=False):
    gene_info_df = pd.read_csv(f'Mammalia/{species}_gene_chr.txt', sep='\t')
    df = pd.read_csv(f'data/umap_genes_{species}_leiden.csv', index_col=0)      
    chromosome_info = {}
    gene_info = {}
    for i, r in gene_info_df.iterrows():
        gene_info[r['Gene stable ID']] = r['Chromosome/scaffold name']
        if r['Gene name']:
            gene_info[r['Gene name']] = r['Chromosome/scaffold name']
        if species == 'human':
            if r['Gene Synonym']:
                gene_info[r['Gene Synonym']] = r['Chromosome/scaffold name']

    percentages = []
    unique_clusters = df['leiden'].unique()
    for uc in unique_clusters:
        if random:
            n = len(df[df['leiden'] == uc])
            df_clus = df.sample(n, axis='index')
        else:
            df_clus = df[df['leiden'] == uc]
        chr_counts  = df_clus.index.map(lambda g: gene_info[g] if g in gene_info else 'NA').value_counts()
        chromosome_info[str(uc)] = (chr_counts.index.values[0], chr_counts.iloc[0] / chr_counts.sum())
        percentages.append(chr_counts.iloc[0] / chr_counts.sum())
    with open(f'data/cluster_chr_perc_{species}.json', 'w') as fw:
        json.dump(chromosome_info, fw, indent=4)
    print(species, 'clusters', 'max:', np.max(percentages), 'min:', np.min(percentages), 'mean:', np.mean(percentages), 'std:', np.std(percentages))
    return percentages


def create_violin_plot(name):
    percentages_human = species_cluster_chr_perc('human')
    percentages_human_random = species_cluster_chr_perc('human', True)
    percentages_mouse = species_cluster_chr_perc('mouse')
    percentages_mouse_random = species_cluster_chr_perc('mouse', True)

    fig = go.Figure()
    for (gs, vec) in [['Human', percentages_human], ['Human random', percentages_human_random], ['Mouse', percentages_mouse], ['Mouse random', percentages_mouse_random]]:
            fig.add_trace(
                go.Violin(
                    y=vec,
                    name=gs.split('_')[0].replace('single', 'GEO'),
                    box_visible=True, 
                    line_color='black',
                    meanline_visible=True,
                    points=False
                )
            )

    fig.update_layout(yaxis_title='Mean Cluster Percentage Same Chromosome', plot_bgcolor='white', yaxis_gridcolor='gray', width=1000, height=800, showlegend=False, font=dict(size=18), yaxis_range=[-.19, 1])
    fig = add_p_value_annotation(fig, [[0, 1], [2, 3]])
    fig.write_image(f'figures/chr_percentage_{name}.pdf')
    fig.write_image(f'figures/chr_percentage_{name}.png')


def enrich_clus_functions(species):
    df = pd.read_csv(f'umap_genes_{species}_leiden.csv', index_col=0)
    unique_clusters = df['leiden'].unique()
    data = []
    for uc in tqdm(unique_clusters):
        df_clus = df[df['leiden'] == uc]
        genes = df_clus.index.tolist()
        link, enrich_res = submit_to_enrichr(genes, ['GO_Biological_Process_2023', 'KEGG_2021_Human', 'Reactome_2022', 'WikiPathway_2023_Human'], f'{species} cluster {uc}')
        r = [f'{species} cluster {uc}', enrich_res, link]
        try:
            data.append([r[0], r[1][0][0][1], r[1][1][0][1], r[1][2][0][1], r[1][3][0][1], r[2]])
        except Exception as e:
            print(e)
            print(r)
    clus_annot_df = pd.DataFrame(data, columns=['Cluster', 'GO_Biological_Process_2023', 'KEGG_2021_Human', 'Reactome_2022', 'WikiPathway_2023_Human', 'Link'])
    clus_annot_df.to_csv(f'data/enrichr_{species}_clusters.tsv', sep='\t')

def graph_clus_yin_yang(species):
    df = pd.read_csv(f'umap_genes_{species}_leiden.csv', index_col=0)
    gmt = GMT.from_file(f'data/{species}-geo-auto.gmt')
    gmt
    


   
if __name__ == '__main__':
    species_cluster_chr_perc('human')
    species_cluster_chr_perc('mouse')


# %%
""" import pandas as pd
from common import GMT
species = 'human'
df = pd.read_csv(f'umap_genes_{species}_leiden.csv', index_col=0)
gmt = GMT.from_file(f'data/{species}-geo-auto.gmt') """


# %%
""" import random
sig_samps = random.sample(range(1, len(gmt.terms) - 1), 1000)
clus_totals = df['leiden'].value_counts().to_dict()
#%%
from tqdm import tqdm
up_down = []

for idx in range(len(gmt.terms)):
    term1 = gmt.terms[idx][0]
    term2 = gmt.terms[idx + 1][0]
    idx2 = idx + 1
    if term2.split(' ')[0] != term1.split(' ')[0]:
        term2 = gmt.terms[idx - 1][0]
        idx2 = idx - 1
        if term2.split(' ')[0] != term1.split(' ')[0]:
            continue
    genes1 = gmt.gene_lists[idx]
    genes2 = gmt.gene_lists[idx2]
    for clus, counts in df.loc[genes1]['leiden'].value_counts()[:10].to_dict().items():
        for clus2, counts2 in df.loc[genes2]['leiden'].value_counts()[:10].to_dict().items():
            if clus != clus2 and counts/clus_totals[clus] >= 0.25 and counts2/clus_totals[clus2] >= 0.25:
                if clus > clus2:
                    up_down.append(f"{clus2} {clus}")
                else:
                    up_down.append(f"{clus} {clus2}") """
# %%
