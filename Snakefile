#outline of pipeline
#   compute MEM
#   compute MUM
#   compute minimizers
#   give anchors to (cahinX) aligner
#       test also the bdwt-mem-chaining

rule compute_MEM:
    input:
        "../br-index-mems/inputs/text.txt.bri",
        "../br-index-mems/inputs/patterns.txt.bri"
    output:
        "../br-index-mems/outputs/MEMs.txt"
    shell:
        "./../br-index-mems/build/bri-mem -k 4 -o {output} {input[0]} {input[1]}"

rule run_ChainX_with_MEMs_as_anchors:
    input:
        "../br-index-mems/inputs/text.txt",
        "../br-index-mems/inputs/patterns.txt",
        "../br-index-mems/outputs/MEMs.txt"
    output:
        "../algbio-ChainX/result/result_new.txt"
    shell:
        "./../algbio-ChainX/chainX -m g -q {input[1]} -t {input[0]} --anchors {input[2]} > {output}"

'''
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
'''