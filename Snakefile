rule all:
    input:
        "phylodata.txt",
        "FASTA"


rule fetch_fasta:
    input:
        r"ena_genome_assembly.xml"
    output:
        directory("FASTA")
    shell:
        "python fasta_from_ena.py -i {input} -o {output}"


rule fetch_metadata:
    input:
        "FASTA"
    output:
        "metadata.txt"
    shell:
        "python fetch_Entrez_metadata.py {input}"


rule phylotype:
    input:
        wd = "FASTA",
        _ = "metadata.txt"
    output:
        "phylodata.txt"
    shell:
        "python utils.py -i {input.wd} -p"


rule serotype:
    input:
        fasta = "FASTA",
        serotype = "serotype"
    output:
        "serodata.txt"
    shell:
        "python utils.py -i {input.fasta} -s {input.serotype}"
    
