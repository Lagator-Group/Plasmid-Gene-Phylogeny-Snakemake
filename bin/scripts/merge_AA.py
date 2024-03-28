# %%
import pandas as pd
import os

# %%
plasmid_df = pd.read_csv(snakemake.input[0])
gene_of_interest = snakemake.wildcards.gene

print(gene_of_interest)

# %%
plasmid_list = []

n=0
for gene_list in plasmid_df['TF Gene Simplified']:
    gene_list = gene_list.replace('[', '').replace(']', '').replace("'", '')
    gene_list = gene_list.split(', ')
    if gene_of_interest in gene_list:
        plasmid_list.append(plasmid_df['Plasmid'][n])
    else:
        pass
    n+=1

# %%
try:
    os.mkdir('merged_AA')
except:
    pass

for plasmid in plasmid_list:

    sprot_df = pd.read_csv(f"sprot_TF/{plasmid}.tsv", sep="\t")
    prokka = f'prokka/{plasmid}/{plasmid}.faa'

    locus_list = []

    n=0
    for gene_list in sprot_df['Gene Names']:
        gene_list = str(gene_list)
        gene_list = gene_list.split(' ')
        if gene_of_interest in gene_list:
            print(prokka)
            locus = sprot_df['Locus Tag'][n]
            locus_list.append(locus)
            with open(prokka, 'r') as f:
                data=f.read()
                for locus in locus_list:
                    start=data.find(locus)
                    end=data.find('>', start+1)
                    result = data[start:end]
                    result = result.replace(' ', '_')

                    print(result)
                    with open(f'merged_AA/{gene_of_interest}.fasta', 'a') as out:
                        out.write(f'>{plasmid}_{result}')
                        out.close()
        else:
            pass
        n+=1    



