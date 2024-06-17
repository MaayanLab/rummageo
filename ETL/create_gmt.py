import pandas as pd
import os
from tqdm import tqdm


def create_gmt(species: str, version: str):
    num_gs =0
    os.makedirs('out/gmts', exist_ok=True)
    sig_files = os.listdir(f'data_{species}_{version}')
    with open(f'out/gmts/{species}-geo-auto_{version}.gmt', 'a') as f:
        for signame in tqdm(list(sig_files)):
            try:
                sig = pd.read_csv(f'data_{species}_{version}/{signame}', index_col=0, sep='\t', compression='gzip')
                sig_signif = sig[sig['adj.P.Val'] < .05]

                genes_up = sig_signif[sig_signif['t'] > 0].index.values

                if len(genes_up) > 2000:
                    sig_signif_up = sig_signif[sig_signif['t'] > 0]
                    cutoff = 0.01
                    i = 1
                    while len(genes_up) > 2000:
                        genes_up = sig_signif_up[sig_signif_up['adj.P.Val'] < cutoff].index.values
                        cutoff = 0.05 * 10**(-i)
                        i += 1

                genes_down = sig_signif[sig_signif['t'] < 0].index.values

                if len(genes_down) > 2000:
                    sig_signif_dn = sig_signif[sig_signif['t'] < 0]
                    cutoff = 0.01
                    i = 1
                    while len(genes_down) > 2000:
                        genes_down = sig_signif_dn[sig_signif_dn['adj.P.Val'] < cutoff].index.values
                        cutoff = 0.05 * 10**(-i)
                        i += 1
                if len(genes_up) >= 5:
                    genes_up_str = '\t'.join(genes_up)
                    f.write(f"{signame.replace('.tsv.gz', '')} up\t\t{genes_up_str}\n")
                    num_gs+=1
                #else: print(sig_signif, len(genes_up), genes_up)
                if len(genes_down) >= 5:
                    genes_down_str = '\t'.join(genes_down)
                    f.write(f"{signame.replace('.tsv.gz', '')} dn\t\t{genes_down_str}\n")
                    num_gs+=1
                #else: print(sig_signif, len(genes_down), genes_down)
            except Exception as e:
                print(sig, e)

    print("Exported", num_gs, "gene sets.")
