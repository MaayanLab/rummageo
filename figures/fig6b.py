import pathlib
import pandas as pd
import random
random.seed(123)
import json
from tqdm import tqdm
import os
import time
from collections import Counter
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import requests

class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        else:
            return super(NumpyEncoder, self).default(obj)
        
def read_gmt(path):
  with pathlib.Path(path).open('r') as fr:
    return {
      term: geneset
      for line in fr
      for term, _, *geneset in (line.strip().split('\t'),)
    }

from scipy.stats import mannwhitneyu

from transformers import AutoModel, AutoTokenizer
from sklearn.metrics.pairwise import cosine_similarity
import torch

os.makedirs('figures/fig6', exist_ok=True)

# Functions for BIOBERT
tokenizer = AutoTokenizer.from_pretrained("microsoft/BiomedNLP-BiomedBERT-base-uncased-abstract-fulltext")
model = AutoModel.from_pretrained("microsoft/BiomedNLP-BiomedBERT-base-uncased-abstract-fulltext")
# Move model to GPU if available
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model.to(device)

def embeddings(texts):
    inputs = tokenizer(texts, padding=True, truncation=False, return_tensors="pt")
    outputs = model(**inputs)
    embeddings = outputs.last_hidden_state.mean(dim=1)
    return embeddings

# Enrichment analysis function
url = "https://rummageo.com/graphql"
def geneset_enrichment(geneset: list):
    query = {
    "operationName": "EnrichmentQuery",
        "variables": {
            "filterTerm": "",
            "offset": 0,
            "first": 200000,
            "genes": geneset,
            "id": "15c56ba6-a293-4932-bcbc-27fc2e4327ab"
        },
        "query": """query EnrichmentQuery($genes: [String]!, $filterTerm: String = "", $offset: Int = 0, $first: Int = 10, $id: UUID!) {
            background(id: $id) {
                id
                species
                enrich(genes: $genes, filterTerm: $filterTerm, offset: $offset, first: $first) {
                nodes {
                    pvalue
                    adjPvalue
                    oddsRatio
                    nOverlap
                    geneSet {
                    id
                    term
                    nGeneIds
                    __typename
                    }
                    __typename
                }
                totalCount
                __typename
                }
                __typename
            }
        }
        """
    }

    headers = {
        "Accept": "application/json",
        "Content-Type": "application/json"
    }

    response = requests.post(url, data=json.dumps(query), headers=headers)

    if response.status_code == 200:
        res = response.json()
        return [(node['geneSet']['term'], node['adjPvalue']) for node in res['data']['background']['enrich']['nodes']]

def term_enrichment(category, matching):
    query = {
    "operationName": "TermEnrichment",
    "variables": {
      "sourceType": category,
      "species": "human",
      "enrichedTerms": matching
      },
    "query": """query TermEnrichment($enrichedTerms: [String], $sourceType: String = \"llm_attrs\", $species: String = \"human\") {
      enrichedFunctionalTerms(
        enrichedTerms: $enrichedTerms
          sourceType: $sourceType
          spec: $species
      ) {
          adjPvalue
          count
          oddsRatio
          notTermCount
          pvalue
          term
          totalNotTermCount
          totalTermCount
          __typename
        }
      }
    """
    }

    headers = {
        "Accept": "application/json",
        "Content-Type": "application/json"
    }

    response = requests.post(url, data=json.dumps(query), headers=headers)

    if response.status_code == 200:
        res = response.json()
        
    return res['data']['enrichedFunctionalTerms']

#### Process GWAS library ####
libraries = ["GWAS_Catalog_2023"]
term_categories = ["A"]
sample = 1000

for lib in libraries:
    print(f'Processing {lib}...')
    gmt = read_gmt(f"enrichr_libraries/{lib}.txt")
    genesets = {term: gmt[term] for term in random.sample(list(gmt), k=sample)}
    
    top_n_gse = 5000 # target number of matching GSEs to include in enrichment analysis
    min_n_gse = 100 # minimum number of matching GSEs to include in enrichment analysis

    matching_genesets = {}
    with tqdm(total=len(genesets)) as pbar:
        for id, gs in genesets.items():
            matching = geneset_enrichment(gs)
            # Extract the unique GSE IDs, preserving the order and keeping only the first occurrence
            unique_matching = []
            seen_gse_ids = set()
            i = 0
            while len(unique_matching) < top_n_gse and i < len(matching):
                tup = matching[i]
                name = tup[0]
                gse_id = name.split('-')[0]  # Extract GSE ID
                if gse_id not in seen_gse_ids:
                    unique_matching.append(tup)
                    seen_gse_ids.add(gse_id)
                i += 1
            if len(unique_matching) >= min_n_gse:
                matching_genesets[id] = unique_matching
            else:
                print(f"{len(unique_matching)} unique matches for {id}")
            pbar.update(1)

    for category_label in term_categories:
        results = {}
        print(f'Computing enrichment for Category {category_label} terms...')
        with tqdm(total=len(matching_genesets)) as pbar:
            for id, matching in matching_genesets.items():
                # Enrich
                results[id] = term_enrichment(category=category_label, matching=[gs[0] for gs in matching])
                pbar.update(1)

        # Save results
        with open(f'results/{lib}_terms_{category_label}.json', 'w') as file:
            json.dump(results, file, indent=4, cls=NumpyEncoder)

#### Process dexamethasone signatures ####
l1000_data_df = pd.read_csv("data/dex_files/Dexamethasone_L1000_ChemPert_data.tsv", sep='\t')
l1000_df_list = []
for row in l1000_data_df.itertuples(): 
    try:
        temp_df = pd.read_csv(row.persistent_id, sep='\t', index_col=0)
    except:
        print(f"Unable to access data from row {row.Index} at {row.persistent_id}")
        continue
    l1000_df_list.append(temp_df)
l1000_signature_df = pd.concat(l1000_df_list, axis=1)
l1000_signature_df.columns = [str(i) for i in range(0, len(l1000_signature_df.columns))]

# Get up- and down- genesets of specified size
size = 250
genesets = {f"{i}_dn": l1000_signature_df[i].sort_values()[:size].index.to_list() for i in l1000_signature_df.columns}

for i in l1000_signature_df.columns:
    genesets[f"{i}_up"] = l1000_signature_df[i].sort_values(ascending=False)[:size].index.to_list()

top_n_gse = 5000 # target number of matching GSEs to include in enrichment analysis
min_n_gse = 100 # minimum number of matching GSEs to include in enrichment analysis

matching_genesets = {}
with tqdm(total=len(genesets)) as pbar:
    for id, gs in genesets.items():
        if id in matching_genesets:
            pbar.update(1)
            continue
        try:
            matching = geneset_enrichment(gs)
        except TypeError as e:
            time.sleep(2)
            matching = geneset_enrichment(gs)
        # Extract the unique GSE IDs, preserving the order and keeping only the first occurrence
        unique_matching = []
        seen_gse_ids = set()
        i = 0
        while len(unique_matching) < top_n_gse and i < len(matching):
            tup = matching[i]
            name = tup[0]
            gse_id = name.split('-')[0]  # Extract GSE ID
            if gse_id not in seen_gse_ids:
                unique_matching.append(tup)
                seen_gse_ids.add(gse_id)
            i += 1
        if len(unique_matching) >= min_n_gse:
            matching_genesets[id] = unique_matching
        else:
            print(f"{len(unique_matching)} unique matches for {id}")
        pbar.update(1)

for category_label in ["B"]:
    results = {}
    print(f'Computing enrichment for Category {category_label} terms...')
    with tqdm(total=len(matching_genesets)) as pbar:
        for id, matching in matching_genesets.items():
            # Enrich
            results[id] = term_enrichment(category=category_label, matching=[gs[0] for gs in matching])
            pbar.update(1)

    # Save results
    with open(f'fig6/dexL1000_terms_{category_label}.json', 'w') as file:
        json.dump(results, file, indent=4, cls=NumpyEncoder)

#### Process random L1000 signatures ####
all_chem_pert_df = read_gmt('data/dex_files/l1000_cp.gmt')
## Remove dexamethasone signatures
for k in list(all_chem_pert_df.keys()):
    if "dexamethasone" in k:
        del all_chem_pert_df[k]

n = len(genesets) # number of random signatures
random_signatures = {k: all_chem_pert_df[k] for k in random.sample(list(all_chem_pert_df.keys()), n)}

matching_genesets = {}
with tqdm(total=len(random_signatures)) as pbar:
    for id, gs in random_signatures.items():
        if id in matching_genesets:
            pbar.update(1)
            continue
        try:
            matching = geneset_enrichment(gs)
        except TypeError as e:
            time.sleep(2)
            matching = geneset_enrichment(gs)
        # Extract the unique GSE IDs, preserving the order and keeping only the first occurrence
        unique_matching = []
        seen_gse_ids = set()
        i = 0
        while len(unique_matching) < top_n_gse and i < len(matching):
            tup = matching[i]
            name = tup[0]
            gse_id = name.split('-')[0]  # Extract GSE ID
            if gse_id not in seen_gse_ids:
                unique_matching.append(tup)
                seen_gse_ids.add(gse_id)
            i += 1
        if len(unique_matching) >= min_n_gse:
            matching_genesets[id] = unique_matching
        else:
            print(f"{len(unique_matching)} unique matches for {id}")
        pbar.update(1)

for category_label in ["B"]:
    results = {}
    print(f'Computing enrichment for Category {category_label} terms...')
    with tqdm(total=len(matching_genesets)) as pbar:
        for id, matching in matching_genesets.items():
            # Enrich
            results[id] = term_enrichment(category=category_label, matching=[gs[0] for gs in matching])
            pbar.update(1)

    # Save results
    with open(f'fig6/randomL1000_terms_{category_label}.json', 'w') as file:
        json.dump(results, file, indent=4, cls=NumpyEncoder)

#### GWAS results ####
lib = "GWAS_Catalog_2023"

with open(f'fig6/{lib}_terms_A.json', 'r') as file:
    gwas_results = json.load(file)

adj_pvals = []
for key, terms in gwas_results.items():
    gwas_enriched_terms[key] = [r['term'] for r in terms]
    adj_pvals.extend([r['adjPvalue'] for r in terms])

# Filter noisy terms
#counts = Counter([term for termlist in gwas_enriched_terms.values() for term in termlist])
#df_counts = pd.DataFrame(counts.items(), columns=['Term', 'Count'])
#thresh = np.percentile(df_counts['Count'], 90)
#common_terms = [term for term, count in counts.items() if count >= thresh]
with open('data/filters/common_terms.txt', 'r') as file: common_terms = file.readlines()
gwas_results = {k: list(filter(lambda r: r['term'] not in common_terms, terms)) for k, terms in gwas_results.items()}

keys = list(gwas_results.keys())
enriched_terms = list(df_counts['Term'])
key_embeddings = embeddings(keys).detach().numpy()
term_embeddings = embeddings(enriched_terms).detach().numpy()
similarities = cosine_similarity(key_embeddings, term_embeddings)

pval_thresh = np.percentile(adj_pvals, 25)

gwas_similarities = {}
gwas_counts = []
shuffled_similarities = [] # baseline comparison
shuffled_counts = []
for key, res in gwas_results.items():
    if len(res) > 0:
        i = keys.index(key)
        terms = [r['term'] for r in res if r['adjPvalue'] < pval_thresh]
        gwas_similarities[key] = []
        count = 0
        rand_count = 0
        random_i = random.randint(0, len(keys)-1)
        for term in terms:
            j = enriched_terms.index(term)
            gwas_similarities[key].append(similarities[i,j])
            shuffled_similarities.append(similarities[random_i, j])
            if similarities[i,j] > 0.95:
                count += 1
            if similarities[random_i, j] > 0.95:
                rand_count += 1
        gwas_counts.append(count)
        shuffled_counts.append(rand_count)

#### Add L1000 dexamethasone results ####
dex_roots = ["dexamethasone", "betamethasone", "prednisolone", "cortiso", "steroid", "nr3c1", "nr0b1", "nr1i2", "annexin", "anxa1", "interleukin", "il10", "il 10", "il1\u03bab", "il6", "il 6", "notch", "nitric oxide synthase", "nos2", "glucocortico", "corticotrop", "anti inflamm", "arthritis", "immunosuppress", "tumor necrosis factor alpha", "tumor necrosis factor \u03b1", "tnf", "tnf \u03b1", "mip2", "mip", "glutathione"]

with open(f'fig6/dexL1000_terms_B.json', 'r') as file:
    dex = json.load(file)

with open(f'fig6/randomL1000_terms_B.json', 'r') as file:
    random_sig = json.load(file)

adj_pvals = []
for key, terms in dex.items():
    dex_enriched_terms[key] = [r['term'] for r in terms]
    adj_pvals.extend([r['adjPvalue'] for r in terms])

# Filter common terms
#counts = Counter([term for termlist in dex_enriched_terms.values() for term in termlist])
#df_counts = pd.DataFrame(counts.items(), columns=['Term', 'Count'])
#thresh = np.percentile(df_counts['Count'], 75)
#common_terms = [term for term, count in counts.items() if count >= thresh]
dex = {k: list(filter(lambda r: r['term'] not in common_terms, terms)) for k, terms in dex.items()}

pval_thresh = np.percentile(adj_pvals, 25)

dex_counts = []
random_sig_counts = []
for key, res in dex.items():
    if len(res) > 0:
        terms = [r['term'] for r in res if r['adjPvalue'] < pval_thresh]
        random_key = random.choice(list(random_sig.keys()))
        random_terms = [r['term'] for r in random_sig[random_key] if r['adjPvalue'] < pval_thresh]
        dex_terms = [s for s in terms if any(root in s for root in dex_roots)]
        random_terms = [s for s in random_terms if any(root in s for root in dex_roots)]
        count = len(dex_terms)
        dex_counts.append(count)
        random_sig_counts.append(len(random_terms))


#### Plot ####
results = [gwas_counts, dex_counts
]

benchmarks = [shuffled_counts, random_sig_counts
]

# Prepare data for seaborn
data = []
labels = []
categories = []
for i, (res, bench) in enumerate(zip(results, benchmarks)):
    data.extend(res)
    data.extend(bench)
    labels.extend([f'Pair {i+1}'] * (len(res) + len(bench)))
    categories.extend(['Results'] * len(res) + ['Benchmark'] * len(bench))

# Create a DataFrame
df = pd.DataFrame({
    'Recovered Terms': data,
    'Benchmark': labels,
    'Category': categories
})

plt.figure(figsize=(6, 6))
palette = {'Results': sns.color_palette()[0], 'Benchmark': 'grey'}
sns.violinplot(x='Benchmark', y='Recovered Terms', hue='Category', data=df, split=True, palette=palette)
ax = plt.gca()

# Set x-tick labels
xtick_labels = ['GWAS Catalog\n(Category A)', 'L1000 Dexamethasone\n(Category B)']
ax.set_xticks([i for i in range(len(results))])
ax.set_xticklabels(xtick_labels)

ax.set_ylim(-2, None)

# Remove the first y-tick label
yticks = ax.get_yticks()
ax.set_yticks(yticks)
ytick_labels = [str(int(tick)) if i > 0 else '' for i, tick in enumerate(yticks)]
ax.set_yticklabels(ytick_labels)

# Legend
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[:2], ['RummaGEO', 'Random'], loc='upper left')

# Perform Mann-Whitney U test and annotate each violin
for i, (res, bench) in enumerate(zip(results, benchmarks)):
    stat, p_value = mannwhitneyu(res, bench, alternative='two-sided')
    x_pos = i+0.25
    y_pos = max(max(res), max(bench))
    if p_value < 0.0001:
        label = 'p < 0.0001'
    else:
        label = f'p={p_value:.4f}'
    ax.text(x_pos, y_pos, label, ha='center', va='bottom', fontsize=12, color='black')

plt.savefig('figures/fig6B.pdf', dpi=400)
