import re
import math
import csv
import pandas as pd
import numpy as np

#TODO generate different reads? congigurable k value

target = "data/test/ecoli.fasta"
query = "data/test/reads/read{r}.fasta"
anchor_types=['mummer-mum', 'mummer-mem', "bdbwt-ext-mini", "bdbwt-mem", 'minimap']
n = 50 #number of reads
k=22 #minumum mem length, k-mer length, hox:minimap2 only allows max 28 as a k-mer size
minimap_k = 28
rule all:
    input:
        expand("results/{anchor_type}-summary.txt", anchor_type = anchor_types),
        expand("results/anchors/{anchor_type}.csv", anchor_type = anchor_types)

#generate reads with simlord
rule generate_reads:
    input: target
    output: "data/test/ecoli_reads.fastq"
    shell: "simlord  --read-reference {input} -n {n} --no-sam  data/test/ecoli_reads"

rule parse_reads:
    input: "data/test/ecoli_reads.fastq"
    output: expand("data/test/reads/read{r}.fasta", r = [_ for _ in range(0,n)])
    run:
        j=0
        with open('data/test/ecoli_reads.fastq') as f:
            for line in f:
                if line[0] == "@":
                    f2 = open(f"data/test/reads/read{j}.fasta", 'w')
                    f2.write(f">{line[1:]}")
                    f2.write(next(f))
                    j += 1
                    f2.close()

def k_length(wildcards):
    if config['constant_k']:
        return k
    read_properties = get_read_properties(f'data/test/reads/read{wildcards.r}.fasta')
    total_error_probability = read_properties[4]
    alpha = -math.log(1-total_error_probability)
    C = (2+total_error_probability)/(1-2*alpha)
    target_length = 4641652
    return f'{int(C*math.log(target_length,4))}'

def minimap_vairable_k_length(wildcards):
    if config['constant_k']:
        return k
    read_properties = get_read_properties(f'data/test/reads/read{wildcards.r}.fasta')
    total_error_probability = read_properties[4]
    alpha = -math.log(1-total_error_probability)
    C = (2+total_error_probability)/(1-2*alpha)
    target_length = 4641652
    if int(C*math.log(target_length,4)) < 28:
        return f'{int(C*math.log(target_length,4))}'
    else:
        return f'28'

#anchor computing
rule benchmark_compute_mem_mummer_variable_k:
    input: query, target
    output: "data/anchors/mummer_mem/read{r}.anchor"
    params: k = k_length
    shell: "./mummer/mummer -maxmatch -l {params.k} {input[0]} {input[1]} >> {output}"

rule benchmark_compute_mum_mummer:
    input: query, target
    output: "data/anchors/mummer_mum/read{r}.anchor"
    params: k = k_length
    shell: "./mummer/mummer -mum -l {params.k} {input[0]} {input[1]} >> {output}"


rule generate_config_file_for_bdbwt_extended_minimizers:
    input: target, query
    output: "bdbwt-mem/configs/ext_min/config{r}"
    params: k = k_length
    run:
        get_congif_bdbwt_mem(output[0], input[0], input[1], params.k, 1)

rule benchmark_compute_extended_minimizers_bdbwt:
    input: "bdbwt-mem/configs/ext_min/config{r}"
    output: "data/anchors/bdbwt-ext-mini/read{r}.anchor"
    shell: "./bdbwt-mem/main {input} >> {output}"

rule generate_config_file_for_bdbwt_MEM:
    input: target, query
    output: "bdbwt-mem/configs/mem/config{r}"
    params: k = k_length
    run:
        get_congif_bdbwt_mem(output[0], input[0], input[1], params.k, 0)


rule benchmark_compute_MEM_bdbwt:
    input: "bdbwt-mem/configs/mem/config{r}"
    output: "data/anchors/bdbwt-mem/read{r}.anchor"
    shell: "./bdbwt-mem/main {input} >> {output}"

rule benchmark_compute_minimap2_minimizers:
    input: query, target
    output: "data/anchors/minimap/read{r}.anchor"
    params: k = minimap_vairable_k_length
    shell: "./minimap2/minimap2 -k {params.k} {input[0]} {input[1]} --print-seeds &>> {output}"

#format the anchors
rule process_mem_mummer_anchors:
    input: "data/anchors/mummer_mem/read{r}.anchor"
    output: "data/anchors/mummer_mem-tidy/read{r}.anchor"
    run:
        f2 = open(output[0], 'w')
        with open(input[0]) as f:
            write2 = False
            for line in f:
                if line[0] != ">":
                    parts = line.split()
                    f2.write(f"{int(parts[1])-1},{int(parts[0])-1},{int(parts[2])}\n")
        f2.close()

use rule process_mem_mummer_anchors as process_mem_mummer_anchors_constant_k with:
    input: "data/anchors/constant_k/mummer_mem/read{r}.anchor"
    output: "data/anchors/constant_k/mummer_mem-tidy/read{r}.anchor"

rule process_mum_mummer_anchors:
    input: "data/anchors/mummer_mum/read{r}.anchor"
    output: "data/anchors/mummer_mum-tidy/read{r}.anchor"
    run:
        f2 = open(output[0], 'w')
        with open(input[0]) as f:
            write2 = False
            for line in f:
                if line[0] != ">":
                    parts = line.split()
                    f2.write(f"{int(parts[1])-1},{int(parts[0])-1},{int(parts[2])}\n")
        f2.close()

use rule process_mum_mummer_anchors as process_mum_mummer_anchors_constant_k with:
    input: "data/anchors/constant_k/mummer_mum/read{r}.anchor"
    output: "data/anchors/constant_k/mummer_mum-tidy/read{r}.anchor"

rule process_bdbwt_extended_min_anchors:
    input: "data/anchors/bdbwt-ext-mini/read{r}.anchor"
    output: "data/anchors/bdbwt-ext-mini-tidy/read{r}.anchor"
    run:
        f2 = open(output[0], 'w')
        with open(input[0]) as f:
            write2 = False
            for line in f:
                if write2:
                    parts = line.split(',')
                    f2.write(f"{int(parts[1])},{int(parts[0])},{int(parts[2])}\n")
                if line.strip() == "MEMs:":
                    write2 = True
        f2.close()

use rule process_bdbwt_extended_min_anchors as process_bdbwt_extended_min_anchors_constant_k with:
    input: "data/anchors/constant_k/bdbwt-ext-mini/read{r}.anchor"
    output: "data/anchors/constant_k/bdbwt-ext-mini-tidy/read{r}.anchor"

rule process_bdbwt_mem_anchors:
    input: "data/anchors/bdbwt-mem/read{r}.anchor"
    output: "data/anchors/bdbwt-mem-tidy/read{r}.anchor"
    run:
        f2 = open(output[0], 'w')
        with open(input[0]) as f:
            write2 = False
            for line in f:
                if write2:
                    parts = line.split(',')
                    f2.write(f"{int(parts[1])},{int(parts[0])},{int(parts[2])}\n")
                if line.strip() == "MEMs:":
                    write2 = True
        f2.close()

rule process_minimap_anchors:
    input: "data/anchors/minimap/read{r}.anchor"
    output: "data/anchors/minimap-tidy/read{r}.anchor"
    run:
        f2 = open(output[0], 'w')
        with open(input[0]) as f:
            for j,line in enumerate(f):
                if line[0:2] == 'SD':
                    parts = line.split()
                    x = int(parts[2])
                    y = int(parts[4])
                    k = int(parts[5])

                    f2.write(f"{y-k+1},{x-k+1},{k}\n")
        f2.close()

#run the chainX with the anchors
rule run_chainX_with_mem_mummer:
    input: target, query, "data/anchors/mummer_mem-tidy/read{r}.anchor"
    output: "data/chains/mummer-mem/chain{r}.txt"
    shell: "./ChainX/chainX -m sg -q {input[1]} -t {input[0]} --anchors {input[2]} >> {output}"

rule benchmark_run_chainX_with_mum_mummer:
    input: target, query, "data/anchors/mummer_mum-tidy/read{r}.anchor"
    output: "data/chains/mummer-mum/chain{r}.txt"
    shell: "./ChainX/chainX -m sg -q {input[1]} -t {input[0]} --anchors {input[2]} >> {output}"

rule benchmark_run_chainX_with_extended_minim:
    input: target, query, "data/anchors/bdbwt-ext-mini-tidy/read{r}.anchor"
    output: "data/chains/bdbwt-ext-mini/chain{r}.txt"
    shell: "./ChainX/chainX -m sg -q {input[1]} -t {input[0]} --anchors {input[2]} >> {output}"

rule benchmark_run_chainX_with_bdbwt_mem:
    input: target, query, "data/anchors/bdbwt-mem-tidy/read{r}.anchor"
    output: "data/chains/bdbwt-mem/chain{r}.txt"
    shell: "./ChainX/chainX -m sg -q {input[1]} -t {input[0]} --anchors {input[2]} >> {output}"

rule benchmark_run_chainX_with_minimizers:
    input: target, query, "data/anchors/minimap-tidy/read{r}.anchor"
    benchmark: "benchmarks/chaining/minimap_mm{r}.tsv"
    output: "data/chains/minimap/chain{r}.txt"
    shell: "./ChainX/chainX -m sg -q {input[1]} -t {input[0]} --anchors {input[2]} >> {output}"

#anchor stats
rule anchor_stats_mummer_mem:
    input: expand("data/anchors/mummer_mem-tidy/read{r}.anchor", r = [_ for _ in range(0,n)])
    output: "results/anchors/mummer-mem.csv"
    params: anchor_type = "mummer-mem"
    run:
        anchor_stats(input, output[0], params.anchor_type)

rule anchor_stats_mummer_mum:
    input: expand("data/anchors/mummer_mum-tidy/read{r}.anchor", r = [_ for _ in range(0,n)])
    output: "results/anchors/mummer-mum.csv"
    params: anchor_type = "mummer-mum"
    run:
        anchor_stats(input, output[0], params.anchor_type)

rule anchor_stats_extended_minimizers:
    input: expand("data/anchors/bdbwt-ext-mini-tidy/read{r}.anchor", r = [_ for _ in range(0,n)])
    output: "results/anchors/bdbwt-ext-mini.csv"
    params: anchor_type = "bdbwt-ext-mini"
    run:
        anchor_stats(input, output[0], params.anchor_type)

#use rule anchor_stats_extended_minimizers as anchor_stats_extended_minimizers_constant_k with:
#    input: expand("data/anchors/contant_k/bdbwt-ext-mini-tidy/read{r}.anchor", r = [_ for _ in range(0,n)])
#    output: "results/anchors/contant_k/bdbwt-ext-mini.csv"

rule anchor_stats_bdbwt_mems:
    input: expand("data/anchors/bdbwt-mem-tidy/read{r}.anchor", r = [_ for _ in range(0,n)])
    output: "results/anchors/bdbwt-mem.csv"
    params: anchor_type = "bdbwt-mem"
    run:
        anchor_stats(input, output[0], params.anchor_type)

rule anchor_stats_minimap_mm:
    input: expand("data/anchors/minimap-tidy/read{r}.anchor", r = [_ for _ in range(0,n)])
    output: "results/anchors/minimap.csv"
    params: anchor_type = "minimap"
    run:
        anchor_stats(input, output[0], params.anchor_type)

#draw some figures, get the results
rule results:
    input: expand("data/chains/{anchor_type}/chain{r}.txt", anchor_type=anchor_types, r = [_ for _ in range(0,n)])
    output: "results/{anchor_type}-summary.txt"
    params: anchor_type = "{anchor_type}"
    run:
        summary(output[0], params.anchor_type)

def summary(summary_output_path, anchor_type):
        input_folder = f'data/chains/{anchor_type}/'
        jaccards = []
        number_of_chains = 0
        total_chain_length = 0
        read_length_sum = 0
        for root, dirs, files in os.walk(input_folder):
            for fi in files:
                read_number = re.findall(r"\d+", fi)[0]
                with open(f'{input_folder}{fi}') as f:
                    chain = []
                    for line in f:
                        try:
                            parts = line[1:len(line)-2].split(',')
                            x = int(parts[0])
                            y = int(parts[1])
                            length = int(parts[2])
                            chain.append((x,y,length))
                        except:
                            continue
                
                chain_length = chain_properties(chain)
                if chain_length > 0:
                    number_of_chains += 1
                    total_chain_length += chain_length
                else:
                    continue
                
                path_to_read_file = f"data/test/reads/read{read_number}.fasta"
                
                with open(path_to_read_file) as f:
                    read_properties = f.readline().split(';')
                    read_start_position = int(re.findall("\d+",read_properties[2])[0])
                    read_length = int(re.findall("\d+", read_properties[1])[0])
                read_length_sum += read_length
                similarity = jaccard_index(read_start_position, read_length, chain)
                jaccards.append(similarity)
        try:
            avg_chain_length = total_chain_length/number_of_chains
        except:
            avg_chain_length = 0
        try:
            avg_accuracy = sum(jaccards)/len(jaccards)
        except:
            avg_accuracy = 0
        try:
            avg_read_length = read_length_sum/number_of_chains
        except:
            avg_read_length = 0
        #the correct alignment \cap the chain alignment/thecorrectalignment + the chain alignment
        f2 = open(f"results/{anchor_type}-summary.txt", 'w')
        f2.write(f'avg chain length: {avg_chain_length}\n')
        f2.write(f'number of chans: {number_of_chains}\n')
        f2.write(f'total chain lenght: {total_chain_length}\n')
        f2.write(f'accuracy: {avg_accuracy}\n')
        f2.write(f'avg read length: {avg_read_length}\n')
        f2.close()

def jaccard_index(read_start_position, read_length, chain):
    read_end_position = read_start_position+read_length
    chain_coverage_of_read = 0
    total_coverage = 0
    total_coverage += read_length
    for (x,y,l) in chain:
        if x >= read_start_position and x+l <= read_end_position: #the anchor is within the orginal read
            chain_coverage_of_read += l
        else:
            total_coverage += l
    return chain_coverage_of_read/total_coverage

def chain_properties(chain):
    chain_length = 0 
    for x,y,l in chain:
        chain_length += l
    return chain_length

def get_congif_bdbwt_mem(output_file_path, target_file_path, query_file_path, k, mode, ):
    f2 = open(output_file_path, 'w')
    f2.write(f"Verbosity > 4\n")
    f2.write(f"Text1	> {query_file_path}\n")
    f2.write(f"Index1	> 1\n")
    f2.write(f"Text2	> {target_file_path}\n")
    f2.write(f"Index2	> 1\n")
    f2.write(f"Mode 	> {mode}\n")
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


def get_read_properties(read_file_path):
    with open(read_file_path) as f:
        read_properties = f.readline().split(';')
        print(read_properties)
        read_length = int(re.findall("\d+", read_properties[1])[0])
        read_start_position = int(re.findall("\d+",read_properties[2])[0])
        read_chromosome = read_properties[3]
        number_of_errors = int(re.findall("\d+",read_properties[4])[0])
        total_error_probability = float(re.findall("\d+.\d+", read_properties[5])[0])
        return read_length, read_start_position, read_chromosome, number_of_errors, total_error_probability

def anchor_stats(input, output, anchor_type):
    sum_of_bases = 0
    number_of_anchors = 0
    
    for file in input:
        with open(file) as f:
            for line in f:
                try:
                    x,y,l = [int(x) for x in line.split(',')]
                    sum_of_bases += l
                    number_of_anchors += 1
                except:
                    continue

    with open(output, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(['type', 'sum of bases', 'number of anchors','avg sum of bases', 'avg number of anchors'])
        csvwriter.writerow([anchor_type, sum_of_bases, number_of_anchors, sum_of_bases/n, number_of_anchors/n])
