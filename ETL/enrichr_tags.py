#%%
import os
import json
from tqdm import tqdm
from maayanlab_bioinformatics.enrichment import enrich_crisp
import multiprocessing as mp
from statsmodels.stats.multitest import multipletests
import pyenrichr as pye

libraries = [
     'ChEA_2022', 'KEGG_2021_Human', 'WikiPathway_2023_Human', 'GO_Biological_Process_2023', 'MGI_Mammalian_Phenotype_Level_4_2021', 'Human_Phenotype_Ontology', 'GWAS_Catalog_2023']
loadedLibs = {
}

#%%
for l in libraries:
    lib = {}
    with open(f'enrichr_libs/{l}.txt') as f:
        lines = f.readlines()
    for line in lines:
        split_line = line.replace('\n', '').split('\t')
        lib[split_line[0]] = set(split_line[2:])
    loadedLibs[l] = lib

fisher = pye.enrichment.FastFisher(34000)

def get_enrichr_labels(term, gene_list, loadedLibs = loadedLibs):
    enrichr_results = {}   
    term = term[0]
    enrichr_results[term] = {}
    gene_list = [g.upper() for g in gene_list]
    for lib in libraries:
        enrich_res = pye.enrichment.fisher(gene_list, loadedLibs[lib], fisher=fisher, background_size=21000)
        sig_enrich_res = enrich_res[enrich_res['fdr'] < 0.05]

        if len(sig_enrich_res) == 0:
            enrichr_results[term][lib] = []
            continue
        sig_enrich_res.reset_index(drop=True, inplace=True)
        enriched_term = []
        for i, r in sig_enrich_res.iterrows():
            enriched_term.append([r['term'], r['p-value'], r['fdr'], r['odds'], r['overlap']])
            if i == 2:
                break
        enrichr_results[term][lib] = enriched_term
    return enrichr_results


def compute_enrichr_labels(species: str, version: str):
        gene_lists = []
        terms = []
        with open(f'out/gmts/{species}-geo-auto_{version}.gmt') as fr:
            lines = fr.readlines()
        for l in tqdm(lines):
            l = l.replace('\n', '')
            split_line = l.split('\t')
            terms.append((split_line[0], split_line[1]))
            gene_lists.append(split_line[2:])
        results = []
        for term, gene_list in tqdm(zip(terms, gene_lists), total=len(terms)):
            results.append(get_enrichr_labels(term, gene_list))
        with open(f'out/enrichr_terms_{species}_{version}.json', 'w') as f:
            json.dump(results, f)