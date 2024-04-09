# This script merges GFF files for plasmids that contain a gene of interest
# from the plasmid_summary.csv file.
# The GFF files are located in prokka/<plasmid>/<plasmid>.gff
# The output is a FASTA file containing the protein sequences for the gene of interest
# for each plasmid that contains it.

import pandas as pd
import os
import shutil

gene_of_interest = snakemake.wildcards.gene
print(gene_of_interest)

# Create the output directories if they don't already exist
try:
    os.mkdir('merged_AA_temp')
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

# For each plasmid in the list of plasmids to process,
# retrieve the sprot_df for that plasmid and extract the gene names and locus tags
for tsv in os.listdir('sprot'):
    plasmid = tsv.replace('.tsv', '')
    print(f"Processing plasmid {plasmid}")
    sprot_df = pd.read_csv(f"sprot/{plasmid}.tsv", sep="\t")

    # Extract the list of gene names and locus tags for this plasmid
    locus_list = []

    n = 0  # index of row in sprot_df
    for gene_list in sprot_df['Gene Names']:
        gene_list = str(gene_list)
        gene_list = gene_list.split(' ')
        print(f"Gene list: {gene_list}")
        if gene_of_interest in gene_list:
            print(f"{gene_of_interest} is in gene list for plasmid {plasmid}")
            # If the gene of interest is in the list of genes for this plasmid,
            # copy the GFF file to the output directory and extract the locus tag
            gff_source = f'prokka/{plasmid}/{plasmid}.gff'
            gff_destination = f'prokka_gff/{gene_of_interest}/{plasmid}.gff'

            print(f"Copying {gff_source} to {gff_destination}")
            shutil.copy(gff_source, gff_destination)

            locus = sprot_df['Locus Tag'][n]
            locus_list.append(locus)
            print(f"Locus tag for {plasmid}: {locus}")
        else:
            pass  # Do nothing if the gene of interest is not in this plasmid's list

        n += 1

    # Read the FASTA file for this plasmid and extract the sequence for each locus tag
    # that was found in the sprot_df
    prokka = f'prokka/{plasmid}/{plasmid}.faa'

    with open(prokka, 'r') as f:
        data = f.read()

        for locus in locus_list:
            start = data.find(locus)
            end = data.find('>', start + 1)
            result = data[start:end]
            result = result.replace(' ', '_')

            # Write the extracted sequence to the output file, along with the plasmid name and locus tag
            with open(f'merged_AA_temp/{gene_of_interest}_temp.fasta', 'a') as out:
                out.write(f'>{plasmid}_{result}\n')
                out.close()
                print(f"Wrote {plasmid}_{result} to output file")

# Open the input and output files
with open(f'merged_AA_temp/{gene_of_interest}_temp.fasta', 'r') as input_file, open(f'merged_AA_temp/{gene_of_interest}.fasta', 'w') as output_file:
    print(f"Reading from {f'merged_AA_temp/{gene_of_interest}_temp.fasta'}...")
    # Read each line from the input file
    for line in input_file:
        print(f"Read line: {line.strip()}")
        # Check if the line contains only '\n'
        if line.strip() != '':
            # Write non-empty lines to the output file
            output_file.write(line)
            print(f"Wrote line to {f'merged_AA_temp/{gene_of_interest}.fasta'}")

print(f"Deleting {f'merged_AA_temp/{gene_of_interest}_temp.fasta'}...")
os.remove(f'merged_AA_temp/{gene_of_interest}_temp.fasta')

