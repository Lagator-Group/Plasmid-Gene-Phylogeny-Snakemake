configfile: 'config.yml'

rule all:
    input:
        expand('data/mobsuite/{plasmid}.tsv',plasmid=config['plasmids']),
        expand('data/panaroo/{gene}',gene=config['genes']),
        expand('data/muscle/protein_{gene}.afa',gene=config['genes']),
        'data/plasmid_grouping/metadata.tsv'

rule merge_AA_gff:
    input: 
        'data/sprot',
    output:
        temp('data/merged_AA_temp/{gene}.fasta'),
        directory('data/prokka_gff/{gene}')
    threads: 2
    script:
        'bin/scripts/merge_AA_gff.py'

rule remove_protein:
    input:
        'data/merged_AA_temp/{gene}.fasta'
    output:
        'data/merged_AA/{gene}.fasta'
    threads: 2
    script:
        'bin/scripts/remove_protein.py'

rule muscle:
    input:
        'data/merged_AA/{gene}.fasta'
    output:
        'data/muscle/protein_{gene}.afa'
    threads: config['threads']
    conda:
        'bin/env/muscle.yml'
    shell:
        'muscle -super5 {input} -output {output}'

rule mobtyper:
    input:
        'data/fasta_plasmid/{plasmid}.fasta'
    output:
        'data/mobsuite/{plasmid}.tsv'
    threads: 4
    conda:
        'bin/env/mobsuite.yml'
    shell:
        'mob_typer --infile {input} --out_file {output}'

rule panaroo:
    input:
        'data/prokka_gff/{gene}'
    output:
        directory('data/panaroo/{gene}')
    threads: 8
    conda:
        'bin/env/panaroo.yml'
    shell:
        'panaroo -i {input}/*.gff -o {output} -t {threads} --clean-mode moderate'

rule metadata_touch:
    input:
        expand('data/mobsuite/{plasmid}.tsv',plasmid=config['plasmids'])
    output:
        temp('data/temp/metadata_touch.txt')
    threads: 1
    shell:
        'touch {output}'

rule metadata:
    input:
        'data/temp/metadata_touch.txt'
    output:
        'data/plasmid_grouping/metadata.tsv'
    threads: 2
    script:
        'bin/scripts/metadata.py'