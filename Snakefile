#outline of pipeline
#   compute MEM
#   compute MUM
#   compute minimizers
#   give anchors to (chainX) aligner

target = "data/test/ecoli.fasta"
query = "data/test/reads/read{r}.fasta"
n = 50 #number of reads
k=22 #minumum mem length, k-mer length, hox:minimap2 only allows max 28 as a k-mer size
rule all:
    input:
        "results/mummer_mem_anchors.txt",
        "results/mummer_mum_anchors.txt",
        "results/bdbwt_ext_mini_anchors.txt",
        "results/bdbwt_mem.txt",
        "results/minimap_mm.txt"

#MEM mummer
rule benchmark_compute_mem_mummer:
    input: query, target
    output: "data/anchors/mummer_mem/read{r}.anchor"
    shell: "./mummer/mummer -maxmatch -l {k} {input[0]} {input[1]} >> {output}"
#MUM mummer
rule benchmark_compute_mum_mummer:
    input: query, target
    output: "data/anchors/mummer_mum/read{r}.anchor"
    shell: "./mummer/mummer -mum -l {k} {input[0]} {input[1]} >> {output}"
#Ext mini from bdbwt
rule generate_config_file_for_bdbwt_extended_minimizers:
    input: query, target
    output: "bdbwt-mem/configs/ext_min/config{r}"
    run:
        f2 = open(output[0], 'w')
        f2.write(f"Verbosity > 4\n")
        f2.write(f"Text1	> {input[0]}\n")
        f2.write(f"Index1	> 1\n")
        f2.write(f"Text2	> {input[1]}\n")
        f2.write(f"Index2	> 1\n")
        f2.write(f"Mode 	> 1\n")
        f2.write(f"Depth	> {k}\n")
        f2.write(f"Window 	> {k}\n")
        f2.write("Mergers > 0\n")
        f2.write("BWTThrd > 1\n")
        f2.write("VerbCA  > 0\n")
        f2.write("VerbED > 0\n")
        f2.write("RawChain > 0\n")
        f2.write("strChain > 0\n")
        f2.write("recombAbs > 0\n")
        f2.write("linearRMQ > 0")
        f2.close()

rule benchmark_compute_extended_minimizers_bdbwt:
    input: "bdbwt-mem/configs/ext_min/config{r}"
    output: "data/anchors/bdbwt-ext-mini/read{r}.anchor"
    shell: "./bdbwt-mem/main {input} >> {output}"

#MEM bdbwt
rule generate_config_file_for_bdbwt_MEM:
    input: query, target
    output: "bdbwt-mem/configs/mem/config{r}"
    run:
        f2 = open(output[0], 'w')
        f2.write(f"Verbosity > 4\n")
        f2.write(f"Text1	> {input[0]}\n")
        f2.write(f"Index1	> 1\n")
        f2.write(f"Text2	> {input[1]}\n")
        f2.write(f"Index2	> 1\n")
        f2.write(f"Mode 	> 0\n")
        f2.write(f"Depth	> {k}\n")
        f2.write(f"Window 	> {k}\n")
        f2.write("Mergers > 0\n")
        f2.write("BWTThrd > 1\n")
        f2.write("VerbCA  > 0\n")
        f2.write("VerbED > 0\n")
        f2.write("RawChain > 0\n")
        f2.write("strChain > 0\n")
        f2.write("recombAbs > 0\n")
        f2.write("linearRMQ > 0")
        f2.close()

rule benchmark_compute_MEM_bdbwt:
    input: "bdbwt-mem/configs/mem/config{r}"
    output: "data/anchors/bdbwt-mem/read{r}.anchor"
    shell: "./bdbwt-mem/main {input} >> {output}"
#minimap minimizers
rule benchmark_compute_minimap2_minimizers:
    input: query, target
    output: "data/anchors/minimap/read{r}.min"
    shell: "./minimap2/minimap2 -k {k} {input[0]} {input[1]} --print-seeds &>> {output}"
#k-mer
#stats
rule benchmark_anchor_stats_mummer_mem:
    input: expand("data/anchors/mummer_mem/read{r}.anchor", r = [_ for _ in range(0,n)])
    output: "results/mummer_mem_anchors.txt"
    benchmark: "benchmarks/anchors/mummer_mem.tsv"
    run:
        sum_of_bases = 0
        number_of_anchors = 0
        for file in input:
            with open(file) as f:
                for line in f:
                    if line[0] != '>':
                        parts = [int(x) for x in line.split()]
                        sum_of_bases += parts[2]
                        number_of_anchors += 1
        f2 = open(output[0], 'w')
        f2.write(f'number of bases: {sum_of_bases}\n')
        f2.write(f'number of anchors: {number_of_anchors}')
        f2.close()

rule benchmark_anchor_stats_mummer_mum:
    input: expand("data/anchors/mummer_mum/read{r}.anchor", r = [_ for _ in range(0,n)])
    output: "results/mummer_mum_anchors.txt"
    benchmark: "benchmarks/anchors/mummer_mum.tsv"
    run:
        sum_of_bases = 0
        number_of_anchors = 0
        for file in input:
            with open(file) as f:
                for line in f:
                    if line[0] != '>':
                        parts = [int(x) for x in line.split()]
                        sum_of_bases += parts[2]
                        number_of_anchors += 1
        f2 = open(output[0], 'w')
        f2.write(f'number of bases: {sum_of_bases}\n')
        f2.write(f'number of anchors: {number_of_anchors}')
        f2.close()

rule benchmark_anchor_stats_extended_minimizers:
    input: expand("data/anchors/bdbwt-ext-mini/read{r}.anchor", r = [_ for _ in range(0,n)])
    output: "results/bdbwt_ext_mini_anchors.txt"
    benchmark: "benchmarks/anchors/extended_minimizers.tsv"
    run:
        sum_of_bases = 0
        number_of_anchors = 0
        for file in input:
            MEMs = False
            with open(file) as f:
                for line in f:
                    if MEMs:
                        parts = [int(x) for x in line.split(',')]
                        sum_of_bases += parts[2]
                        number_of_anchors += 1
                    if line.strip() == "MEMs:":
                        MEMs = True

        f2 = open(output[0], 'w')
        f2.write(f'number of bases: {sum_of_bases}\n')
        f2.write(f'number of anchors: {number_of_anchors}')
        f2.close()

rule benchmark_anchor_stats_bdbwt_mems:
    input: expand("data/anchors/bdbwt-mem/read{r}.anchor", r = [_ for _ in range(0,n)])
    output: "results/bdbwt_mem.txt"
    benchmark: "benchmarks/anchors/bdbwt_mem.tsv"
    run:
        sum_of_bases = 0
        number_of_anchors = 0
        for file in input:
            MEMs = False
            with open(file) as f:
                for line in f:
                    if MEMs:
                        parts = [int(x) for x in line.split(',')]
                        sum_of_bases += parts[2]
                        number_of_anchors += 1
                    if line.strip() == "MEMs:":
                        MEMs = True

        f2 = open(output[0], 'w')
        f2.write(f'number of bases: {sum_of_bases}\n')
        f2.write(f'number of anchors: {number_of_anchors}')
        f2.close()

rule benchmark_anchor_stats_minimap_mm:
    input: expand("data/anchors/minimap/read{r}.min", r = [_ for _ in range(0,n)])
    output: "results/minimap_mm.txt"
    benchmark: "benchmarks/anchors/minimap_mm.tsv"
    run:
        sum_of_bases = 0
        number_of_anchors = 0
        for file in input:
            with open(file) as f:
                for line in f:
                    if line[0:2] == "SD":
                        parts = line.split()
                        sum_of_bases += int(parts[5])
                        number_of_anchors += 1

        f2 = open(output[0], 'w')
        f2.write(f'number of bases: {sum_of_bases}\n')
        f2.write(f'number of anchors: {number_of_anchors}')
        f2.close()

rule process_min_anchors:
    input: "data/anchors/minimap/read{r}.min"
    output: "data/anchors/minimap-tidy/read{r}.min"
    run:
        f2 = open(output[0], 'w')
        with open(input[0]) as f:
            for j,line in enumerate(f):
                if line[0:2] == 'SD':
                    parts = line.split()
                    x = int(parts[2])
                    y = int(parts[4])
                    k = int(parts[5])

                    f2.write(f"{x-k+1},{y-k+1},{k}\n")
        f2.close()

rule run_chainX_with_minimizers:
    input: target, query, "data/anchors/minimap-tidy/read{r}.min"
    output: "results/minimap/read{r}_result.txt"
    shell: "./ChainX/chainX -m g -q {input[1]} -t {input[0]} --anchors {input[2]} &>> {output}"



rule process_mem_anchors:
    input: "data/anchors/bdbwt-mini/{r}.minie"
    output: "data/anchors/bdbwt-mini-tidy/{r}.minie"
    run:
        f2 = open(output[0], 'w')
        with open(input[0]) as f:
            write2 = False
            for line in f:
                print(line)
                if write2:
                    f2.write(line)
                if line.strip() == "MEMs:":
                    write2 = True
        f2.close()

rule run_chainX_with_extended_minim:
    input: target, query, "data/anchors/bdbwt-mini-tidy/{r}.minie"
    output: "results/bdbwt-minie/result{r}.txt"
    shell: "./ChainX/chainX -m g -q {input[1]} -t {input[0]} --anchors {input[2]} &>> {output}"

#draw some figures


