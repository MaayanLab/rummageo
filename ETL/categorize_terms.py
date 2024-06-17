import jellyfish
import pandas as pd
import json


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

    cat_df = pd.read_csv('data/LLM_keyterm_categories.csv')

    terms = set(cat_df['terms'].values)

    for t in new_keyterms:
        if t in terms:
            continue
        else:
            cat_df['sim'] = cat_df['terms'].apply(lambda x: jellyfish.levenshtein_distance(x, t))
            term_sim = cat_df.sort_values('sim', ascending=False)
            new_categorizations[t] = term_sim['manual_category'][:100].value_counts().idxmax()
    
    with open(f'out/keyterms/key_terms_categorized_{species}_{version}.json', 'w') as f:
        json.dump(new_categorizations, f, indent=4)
    
    return new_categorizations


            
    

# %%
