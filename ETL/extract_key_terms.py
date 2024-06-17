import json
import pandas as pd
import requests
from transformers import AutoTokenizer
import pandas as pd
from tqdm import tqdm
from datetime import datetime
import time
import codecs
import os
import re
from itertools import islice
import json
import requests
import pandas as pd
from tqdm import tqdm
import time
import xml.etree.ElementTree as ET
from urllib.error import HTTPError
from socket import error as SocketError
from Bio import Entrez
import jellyfish
import GEOparse
from dotenv import load_dotenv
load_dotenv()

os.makedirs('out/keyterms', exist_ok=True)


API_URL = "https://api-inference.huggingface.co/models/mistralai/Mistral-7B-Instruct-v0.2"
api_token = os.getenv('HF_API_KEY') # Get yours at hf.co/settings/tokens
# Log into huggingface using huggingface-cli login
headers = {"Authorization": f"Bearer {api_token}"} 

tokenizer = AutoTokenizer.from_pretrained("mistralai/Mistral-7B-Instruct-v0.2")

def query(payload):
	response = requests.post(API_URL, headers=headers, json=payload)
	return response.json()

# function to format any weird responses
def format_string(s):
    # Strip leading and trailing whitespace
    s = s.strip()
    # Get only the first sentence
    s = str.split(s, ".")[0] 
    # Replace any commas with semicolons
    s = s.replace(',', ';')
    # Replace the last period with ']' and add '[' at the start
    if not s.endswith('.'):
        s += '.'
    s = '[' + s[:-1] + ']'
    return s

def string_process(s):
    # parse and extract text using BeautifulSoup
    #processed_s = BeautifulSoup(s, "html.parser").get_text()
    # remove newline characters
    processed_s = s.replace("\n", "").replace("\r", "")
    # remove extra spaces
    processed_s = ' '.join(processed_s.split())
    # to lower case
    processed_s = processed_s.lower()
    # dashes
    processed_s = re.sub(r'[-–—]', ' ', processed_s)
    # remove punctuation
    processed_s = re.sub(r'[^\w\s()/+]', '', processed_s)
    return processed_s

def get_gse_info(gse_list, species):
    gse_info = {}
    for gse in tqdm(gse_list, desc='Fetching GSE info...'):
        if ',' in gse:
            gse_split = gse.split(',')[0]
        else:
            gse_split = gse
        for i in range(10):
            try:
                geo_meta = GEOparse.GEOparse.get_GEO(geo=gse_split, silent=True, destdir='../data/geo')
            except:
                print(f'Failed to fetch {gse}')
                continue
            gse_info[gse] = {}
            gse_info[gse]['title'] = geo_meta.get_metadata_attribute('title')
            gse_info[gse]['summary'] = geo_meta.get_metadata_attribute('summary')
            try:
                gse_info[gse]['pmid'] = geo_meta.get_metadata_attribute('pubmed_id')
            except:
                gse_info[gse]['pmid'] = None

            date_object = datetime.strptime(geo_meta.get_metadata_attribute('status'), "Public on %b %d %Y")
            formatted_date = date_object.strftime("%Y-%m-%d")
            gse_info[gse]['publication_date'] = formatted_date  # yyyy-mm-dd
            gse_info[gse]['platform'] = geo_meta.get_metadata_attribute('platform_id')
            gse_info[gse]['species'] = species
    return gse_info


def fetch_pubmed_abstract(pmid: str):
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        'db': 'pubmed',
        'id': pmid,
        'retmode': 'json',
        'rettype': 'abstract',
        'api_key': os.getenv('PB_API_KEY')
    }

    response = requests.get(base_url, params=params)
    try:
        res = ''.join(response.text.strip().split('\n\n')[4]).replace('\n', ' ')
    except IndexError:
        print(pmid, response.text)
        res = ''
    return res


def mistral_keywords(gse_id: str, abstract: str):
    example = """Example abstract: \nAge-related macular degeneration (AMD) is a leading cause of blindness, affecting 200 million people worldwide. To identify genes that could be targeted for treatment, we created a molecular atlas at different stages of AMD. Our resource is comprised of RNA sequencing (RNA-seq) and DNA methylation microarrays from bulk macular retinal pigment epithelium (RPE)/choroid of clinically phenotyped normal and AMD donor eyes (n = 85), single-nucleus RNA-seq (164,399 cells), and single-nucleus assay for transposase-accessible chromatin (ATAC)-seq (125,822 cells) from the retina, RPE, and choroid of 6 AMD and 7 control donors. We identified 23 genome-wide significant loci differentially methylated in AMD, over 1,000 differentially expressed genes across different disease stages, and an AMD Müller state distinct from normal or gliosis. Chromatin accessibility peaks in genome-wide association study (GWAS) loci revealed putative causal genes for AMD, including HTRA1 and C6orf223. Our systems biology approach uncovered molecular mechanisms underlying AMD, including regulators of WNT signaling, FRZB and TLE2, as mechanistic players in disease.\n
    Example Output: \n[age-related macular degeneration; geographic atrophy; single-cell RNA-seq; single-cell ATAC-seq; rare variant genetics; Muller glia; retinal pigment epithelium]\n"""
    prompt = f"Your role is to extract biomedical keyterms from research abstracts. Return ONLY the list of terms formatted as [term1; term2; term3], with at most 10 terms.\n\nHere's an example:\n\n{example}\n\nNow, do this for the following abstract. {abstract}"
    conv = [
            {"role": "user", "content": prompt}]
    inputs = tokenizer.apply_chat_template(conv, tokenize=False, add_generation_prompt=True)

    retry_count = 0
    terms = None
    while retry_count < 2:
        try:
            output = query({
                "inputs": inputs,
                "parameters": {
                    "max_new_tokens": 1000,
                    "return_full_text": False,
                }
            })
            print(output)
            if 'error' in output:
                if 'rate limit' in output['error'].lower():
                    print('Rate limit exceeded, sleeping for 1hr...')
                    time.sleep(3600)
                    continue
            generated_text = output[0]['generated_text']
            generated_text = codecs.escape_decode(bytes(generated_text,"utf-8"))[0].decode("utf-8") 
            terms = re.search(r'\[(.*?)\]', generated_text).group(1)
        except AttributeError as e: 
            # list is missing close bracket
            if generated_text[-1] == ';':
                generated_text = generated_text[:-1] + ']'
                terms = re.search(r'\[(.*?)\]', generated_text).group(1)
                retry_count += 2 # finish
            # list is formatted in some other way
            else:
                terms = format_string(generated_text)
                retry_count += 2 # finish
        except (SyntaxError, TypeError, ValueError) as e:
            retry_count += 1
            print(f'Error {e} for {gse_id}, retrying in 2s...')
            # log the generated text
            time.sleep(2)
        except (MemoryError, OverflowError) as e:
            # break the loop
            print(f'Memory error {e} for {gse_id}')
            break
        except Exception as e:
            print(f'Error {e} for {gse_id}, retrying in 2s...')
            retry_count += 1
            time.sleep(2)
        else:
            break
    return terms


term_filter = set(['transcription factor', 'single cell rna seq', 
             'mirna', 'enhancer', 'long non coding rna', 'cancer', 
             'differentiation', 'metabolism', 'cell proliferation',
             'rna polymerase ii', 'histone modification', 'histone acetylation'
             'translation', 'mouse model', 'pluripotent stem cell', 'protein synthesis', 
             'human embryonic stem cell', 'transcriptional regulation', 'chromatin remodeling', 
             'differentially expressed gene', 'chip seq', 'atac seq', 'self renewal', 'mouse embryonic stem cell'
             'epigenome', 'gene regulation', 'interferon Œ≥', 'chromatin', 'gene ontology', 'stem cell', 'super enhancer', 
             'transcriptional profiling', 'protein protein interaction', 'proteomic', 'gene expression profile',
             'metabolomic', 'histone deacetylase', 'signaling pathway', 'transcriptome analysis', 'mrna',
             'drug resistance', 'chromatin structure', 'next generation sequencing', 'chromatin immunoprecipitation', 
             'mrna stability', 'chemotherapy', 'transcriptional program', 'chromatin immunoprecipitation sequencing', 
             'metabolic reprogramming', 'epigenetic reprogramming', 'single cell transcriptomic', 'proteome',
             'chromatin binding', 'reprogramming', 'gene regulatory network', 'transcriptomic analysis', 'chromatin organization', 
             'tumor progression', 'bulk rna seq', 'promoter', 'in vitro', 'small molecule', 'genome wide association studies',
             'mutation', 'gene expression profiling', 'organoid', 'transcriptional response', 'clonal expansion',
             'cancer cell', 'amino acid metabolism', 'resistance', 'in vivo', 'chromatin landscape', 'histone h3',
             'machine learning', 'transcriptional profile', 'transcriptional change', 'gene expression program',
             'transcriptome profiling', 'post translational modification', 'hub gene', 'tissue repair', 'single cell analysis', 
             'transcriptional network', 'regulatory elements', 'protein translation', 'lineage plasticity', 'histone methylation',
             'pluripotent cells', 'development', 'histone mark', 'mrna translation', 'transcriptional repression', 'transcriptional signature',
             'cancer cells', 'metabolic pathway', 'rna interference', 'gene set enrichment analysis', 'enhancer promoter interaction', 
             'non coding rna', 'transcription regulation', 'differential expression', 'single cell transcriptomics', 
             'single cell transcriptome', 'protein expression', 'cellular response', 'genetic variation', 'gene network',
             'transcriptional activity', 'differential gene expression', 'mass spectrometry', 'germ free mice', 'gene expression analysis',
             'gene expression regulation', 'cellular proliferation', 'copy number variation', 'cancer cell lines', 
             'genome wide association study', 'gene transcription', 'cell signaling', 'rna binding protein', 'gene expression data',
             'cancer immunotherapy', 'smart seq', 'age', 'mrna sequencing', 'inhibitory receptors', 'mice', 'drug response', 'selectivity',
             'gene expression pattern', 'gain of function', 'in vitro assay', 'protein coding genes', 'gene expression changes', 
             'transcription factor network', 'transcriptomic changes', 'biological pathways', 'mouse development', 
             'regulatory regions', 'signaling', 'gene editing', 'homeobox gene', 'gene expression comparisons', 'rna seq', 'gene expression', 'transcriptome', 'human', 'transcription', 'mouse', 'gene expression signature', 'genetic', 'genomic', 'biological sciences'])


def clean_terms(terms, plurals_dict, LLM_ss_dict, manual_dict):
    if terms == []:
        return []
    term_list = [string_process(s) for s in terms.split('; ')]
    term_list = [s for s in term_list if len(s) > 2]
    term_list = [plurals_dict.get(term, term) for term in term_list]
    prev_state = []
    while term_list != prev_state: 
        prev_state = term_list
        term_list = []
        for term in prev_state:
            # Replace terms with their synonym
            mapped = LLM_ss_dict.get(term, term)
            if isinstance(mapped, list):
                for x in LLM_ss_dict.get(term, term):
                    term_list.append(x)
            elif isinstance(mapped, str):
                term_list.append(mapped)
    replaced_list = [manual_dict.get(term, term) for term in term_list] # Replace terms with manual fix
    replaced_list_cleaned = [term for term in replaced_list if term not in term_filter] # Remove general terms
    return list(set(replaced_list_cleaned))


def categorize_terms(species, version):
    """
    Categorize a list of terms into a list of categories based on the Levenshtein distance.
    :param terms: list of terms to categorize
    :param categories: list of categories
    :return: dictionary with categories as keys and lists of terms as values
    """
    new_categorizations = {}

    with open(f'out/keyterms/gse_key_terms_clean_{species}_{version}.json') as f:
        new_keyterms = json.load(f)

    cat_df = pd.read_csv('../data/LLM_keyterm_categories.csv')

    terms = set(cat_df['term'].values)

    for gse in tqdm(new_keyterms, 'Categorizing terms...', total=len(new_keyterms)):
        for t in new_keyterms[gse]:
            if t in terms:
                continue
            else:
                cat_df['sim'] = cat_df['term'].apply(lambda x: jellyfish.jaro_winkler_similarity(x, t))
                term_sim = cat_df.sort_values('sim', ascending=False)
                new_categorizations[t] = term_sim['manual_category'][:100].value_counts().idxmax()
    return new_categorizations

def generate_key_terms(species, version):
    with open(f'out/meta/gse_processed_meta_{species}_{version}_conf.json', 'r') as f:
        gse_info_conf = json.load(f)
    gse_list = list(gse_info_conf.keys())
    if not os.path.exists(f'out/meta/gse_info_{species}_{version}.json'):
        gse_info = get_gse_info(gse_list, species)
        with open(f'out/meta/gse_info_{species}_{version}.json', 'w') as f:
            json.dump(gse_info, f)
    else:
        with open(f'out/meta/gse_info_{species}_{version}.json', 'r') as f:
            gse_info = json.load(f)

    if os.path.exists(f'out/keyterms/gse_key_terms_{species}_{version}.json'):
        with open(f'out/keyterms/gse_key_terms_{species}_{version}.json', 'r') as f:
            key_terms = json.load(f)
    else:
        key_terms = {}
    i = 0
    for gse in tqdm(gse_list, desc='Extracting key terms...'):
        if gse in key_terms and key_terms[gse] and key_terms[gse] != []:
            continue
        if gse not in gse_info:
            key_terms[gse] = []
            print(gse)
            continue
        pmid = gse_info[gse]['pmid']
        if pmid:
            abstract = fetch_pubmed_abstract(pmid)
            res = mistral_keywords(gse, abstract)   
        else:
            res = mistral_keywords(gse, gse_info[gse]['summary'])
        key_terms[gse] = res
        i += 1
        if i % 20 == 0:
            with open(f'out/keyterms/gse_key_terms_{species}_{version}.json', 'w') as f:
                json.dump(key_terms, f)

    with open(f'out/keyterms/gse_key_terms_{species}_{version}.json', 'w') as f:
        json.dump(key_terms, f)

    with open(f'keyterm-maps/{species}/plurals.json', 'r') as f:
        plurals_dict = json.load(f)
    with open(f'keyterm-maps/{species}/LLM_substrings.json', 'r') as file:
        LLM_ss_dict = json.load(file) 
    with open(f'keyterm-maps/manual_map.json', 'r') as file:
        manual_dict = json.load(file)

    key_terms_clean = {}
    for gse in key_terms:
        cleaned_terms = clean_terms(key_terms[gse], plurals_dict, LLM_ss_dict, manual_dict)
        key_terms_clean[gse] = cleaned_terms

    with open(f'out/keyterms/gse_key_terms_clean_{species}_{version}.json', 'w') as f:
        json.dump(key_terms_clean, f)

    new_categorizations = categorize_terms(species, version)

    with open(f'out/keyterms/key_terms_categorized_{species}_{version}.json', 'w') as f:
        json.dump(new_categorizations, f, indent=4)


generate_key_terms('human', '2.4')

generate_key_terms('mouse', '2.4')