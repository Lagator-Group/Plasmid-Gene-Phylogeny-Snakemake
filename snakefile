configfile: 'config.yml'

rule all:
    input:
        expand('iqtree/protein_{gene}/{gene}.tree',gene=config['genes']),
        expand('mobsuite/{plasmid}.tsv',plasmid=config['plasmids']),
        expand('panaroo/{gene}',gene=config['genes']),
        'plasmid_grouping/metadata.tsv'

rule merge_AA_gff:
    input: 
        'plasmid_summary.csv',
    output:
        temp('merged_AA_temp/{gene}.fasta'),
        directory('prokka_gff/{gene}')
    threads: 2
    script:
        'bin/scripts/merge_AA_gff.py'

rule remove_protein:
    input:
        'merged_AA_temp/{gene}.fasta'
    output:
        'merged_AA/{gene}.fasta'
    threads: 2
    script:
        'bin/scripts/remove_protein.py'

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
        'iqtree/protein_{gene}/{gene}.tree'
    params:
        prefix = 'iqtree/protein_{gene}/{gene}',
        treefile = 'iqtree/protein_{gene}/{gene}.treefile'
    threads: 8
    conda:
        'bin/env/iqtree.yml'
    shell:
        'iqtree -s {input[0]} -st AA -nt {threads} --prefix {params.prefix} &&'
        'mv {params.treefile} {output}'

rule mobtyper:
    input:
        'fasta_plasmid/{plasmid}.fasta'
    output:
        'mobsuite/{plasmid}.tsv'
    threads: 4
    conda:
        'bin/env/mobsuite.yml'
    shell:
        'mob_typer --infile {input} --out_file {output}'

rule panaroo:
    input:
        'prokka_gff/{gene}'
    output:
        directory('panaroo/{gene}')
    threads: 8
    conda:
        'bin/env/panaroo.yml'
    shell:
        'panaroo -i {input}/*.gff -o {output} -t {threads} --clean-mode moderate'

rule metadata:
    input:
    output:
        'plasmid_grouping/metadata.tsv'
    threads: 2
    script:
        'bin/scripts/metadata.py'