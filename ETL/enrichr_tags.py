#%%
import os
import json
from tqdm import tqdm
from maayanlab_bioinformatics.enrichment import enrich_crisp
import multiprocessing as mp
from statsmodels.stats.multitest import multipletests


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
        lib[split_line[0]] = split_line[2:]
    loadedLibs[l] = lib

def get_enrichr_labels(term, gene_list, loadedLibs = loadedLibs):
    enrichr_results = {}   
    term = term[0]
    enrichr_results[term] = {}
    gene_list = [g.upper() for g in gene_list]
    for lib in libraries:
        enrich_res = enrich_crisp(gene_list, loadedLibs[lib], n_background_entities=21000)
        res_sorted = list(sorted(enrich_res, key=lambda x: x[1].pvalue))
        p_values = list(map(lambda x: x[1].pvalue, res_sorted))
        if res_sorted == []:
            enrichr_results[term][lib] = []
            continue
        try:
            adjusted_p_values = multipletests(p_values, method='fdr_bh', is_sorted=True)[1]
        except Exception as e:
            print(e, p_values, res_sorted, term, gene_list, lib)
            adjusted_p_values = [1, 1, 1]
        enriched_term = []
        for i, r in enumerate(res_sorted[:3]):
            enriched_term.append([r[0], r[1].pvalue, adjusted_p_values[i], r[1].odds_ratio, r[1].n_overlap])
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
        
        if os.path.exists(f'out/enrichr_terms_{species}_{version}.json'):
            with open(f'out/enrichr_terms_{species}_{version}.json') as f:
                enrichr_terms = json.load(f)
            terms_to_enrich = []
            genesets_to_enrich = []
            enriched_terms_done = set([item.keys()[0] for item in enrichr_terms])
            for i in range(len(terms)):
                if terms[i][0] not in enriched_terms_done:
                    terms_to_enrich.append(terms[i])
                    genesets_to_enrich.append(gene_lists = gene_lists[i])
            inputs = zip(terms_to_enrich, genesets_to_enrich)
            with mp.Pool(os.cpu_count()) as pool:
                results = pool.starmap(get_enrichr_labels, tqdm(inputs, total=len(terms)))
            with open(f'out/enrichr_terms_{species}_{version}.json', 'w') as f:
                json.dump(results + enrichr_terms, f)
        else:
            inputs = zip(terms, gene_lists)
            with mp.Pool(os.cpu_count() // 2) as pool:
                results = pool.starmap(get_enrichr_labels, tqdm(inputs, total=len(terms)))
            with open(f'out/enrichr_terms_{species}_{version}.json', 'w') as f:
                json.dump(results, f)