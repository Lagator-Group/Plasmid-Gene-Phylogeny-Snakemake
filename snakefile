configfile: 'config.yml'

rule all:
    input:
        expand('muscle/protein_{gene}.afa',gene=config['genes'])

rule get_protein_faa:
    input: 
        'plasmid_summary.csv',
    output:
        'merged_AA/{gene}.fasta'
    threads: 2
    script:
        'bin/scripts/merge_AA.py'

rule muscle:
    input:
        'merged_AA/{gene}.fasta'
    output:
        'muscle/protein_{gene}.afa'
    threads: 8
    conda:
        'muscle'
    shell:
        'muscle -align {input} -output {output}'

