configfile: 'config.yml'

rule all:
    input:
        expand('merged_AA/{gene}.fasta',gene=config['genes'])

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
        'bin/env/muscle.yml'
    shell:
        'muscle -align {input} -output {output}'

rule mkdir_iqtree:
    input:
    output:
        'iqtree/protein_{gene}/{gene}.log'
    shell:
        'touch {output}'
    
rule iqtree:
    input:
        'muscle/protein_{gene}.afa',
        'iqtree/protein_{gene}/{gene}.log'
    output:
        'iqtree/protein_{gene}/{gene}.iqtree'
    params:
        prefix = 'iqtree/protein_{gene}/{gene}'
    threads: 8
    conda:
        'bin/env/iqtree.yml'
    shell:
        'iqtree -s {input[0]} -st AA -nt {threads} --prefix {params.prefix}'

