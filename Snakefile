rule run_ChainX:
    input:
        "../algbio-ChainX/data/chromosome.fasta",
        "../algbio-ChainX/data/mutated.fasta",
        "../algbio-ChainX/data/anchors.mems"
    output:
        "../algbio-ChainX/result/result.txt"
    shell:
        "./../algbio-ChainX/chainX -m g -q {input[1]} -t {input[0]} --anchors {input[2]} > {output}"


rule run_bdbwtmemchaining:
    input:
        "../algbio-ChainX/data/chromosome.fasta",
        "../algbio-ChainX/data/mutated.fasta"
    output:
        "../bdbwt-mem-chaining/result/result.txt"
    shell:
        "./../bdbwt-mem-chaining/main ../bdbwt-mem-chaining/config > {output}"