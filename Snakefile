#outline of pipeline
#   compute MEM
#   compute MUM
#   compute minimizers
#   give anchors to (chainX) aligner

rule compute_mem_bdwt:
    input: "data/test/MT-human.fa", "data/test/MT-orang.fa"
    output: "data/anchors/mems.mem"
    shell: "./bdbwt-mem/main bdbwt-mem/config >> {output}"

rule compute_minimizers:
    input: "data/test/MT-human.fa", "data/test/MT-orang.fa"
    output: "data/anchors/minimizers.min"
    shell: "./minimap2/minimap2 data/test/MT-human.fa data/test/MT-orang.fa --print-seeds &>> {output}"

rule process_mem_anchors:
    input: "data/anchors/mems.mem"
    output: "data/anchors/tidy-mems.mem"
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

rule process_min_anchors:
    input: "data/anchors/minimizers.min"
    output: "data/anchors/tidy-minimizers.min"
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
    input: "data/test/MT-human.fa", "data/test/MT-orang.fa", "data/anchors/tidy-minimizers.min"
    output: "results/result_minimizer.txt"
    shell: "./ChainX/chainX -m g -q {input[1]} -t {input[0]} --anchors {input[2]} >> {output}"

rule run_chainX_with_mems:
    input: "data/test/MT-human.fa", "data/test/MT-orang.fa", "data/anchors/tidy-mems.mem"
    output: "results/result_mems.txt"
    shell: "./ChainX/chainX -m g -q {input[1]} -t {input[0]} --anchors {input[2]} >> {output}"
