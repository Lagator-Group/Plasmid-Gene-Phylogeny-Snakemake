# %%
import pandas as pd
import os
import shutil

# %%
plasmid_df = pd.read_csv(snakemake.input[0])
gene_of_interest = snakemake.wildcards.gene

'''#debugging
plasmid_df = pd.read_csv('plasmid_summary.csv')
gene_of_interest = 'finO'
'''

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

try:
    os.mkdir('prokka_gff')
except:
    pass

try:
    os.mkdir(f'prokka_gff/{gene_of_interest}')
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

            gff_source = f'prokka/{plasmid}/{plasmid}.gff'
            gff_destination = f'prokka_gff/{gene_of_interest}/{plasmid}.gff'

            shutil.copy(gff_source, gff_destination)

            print(prokka)
            locus = sprot_df['Locus Tag'][n]
            locus_list.append(locus)
        else:
            pass
        n+=1    

    with open(prokka, 'r') as f:
        data=f.read()
        
        for locus in locus_list:
            start=data.find(locus)
            end=data.find('>', start+1)
            result = data[start:end]
            result = result.replace(' ', '_')
            print(result)

            with open(f'merged_AA/{gene_of_interest}_temp.fasta', 'a') as out:
                out.write(f'>{plasmid}_{result}\n')
                out.close()

# Open the input and output files
with open(f'merged_AA/{gene_of_interest}_temp.fasta', 'r') as input_file, open(f'merged_AA/{gene_of_interest}.fasta', 'w') as output_file:
    # Read each line from the input file
    for line in input_file:
        # Check if the line contains only '\n'
        if line.strip() != '':
            # Write non-empty lines to the output file
            output_file.write(line)

os.remove(f'merged_AA/{gene_of_interest}_temp.fasta')



