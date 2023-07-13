import re
import math
from collections import deque
import os

target = "data/test/ecoli.fasta"
query = "data/test/reads/read{r}.fasta"
anchor_types=['mummer-mum', 'mummer-mem', "bdbwt-ext-mini", "bdbwt-mem", 'minimap']
k=22 #minumum mem length, k-mer length, hox:minimap2 only allows max 28 as a k-mer size
read_path = 'data/test/reads/read{r}.fasta'
r_number = glob_wildcards(read_path).r
minimap_k = 28
target_length = 0
with open(target) as f:
    for j,line in enumerate(f):
        if j != 0:
            target_length += len(line.strip())

rule all:
    input:
        expand("results/{anchor_type}-summary.txt", anchor_type = anchor_types),
        expand("results/anchors/{anchor_type}.txt", anchor_type = anchor_types)

#generate reads with simlord
rule generate_reads:
    input: target
    output: "data/test/ecoli_reads.fastq"
    shell: "simlord --read-reference {input} -c 1.0 --no-sam data/test/reads"

rule parse_reads:
    input: "data/test/reads.fastq"
    output: "data/test/reads/read{r}.fasta"
    run:
        j=0
        with open('data/test/reads.fastq') as f:
            for line in f:
                if line[0] == "@":
                    f2 = open(f"data/test/reads/read{j}.fasta", 'w')
                    f2.write(f">{line[1:]}")
                    f2.write(next(f))
                    j += 1
                    f2.close()

def k_size(wildcards):
    if config['constant_k']:
        return k
    read_properties = get_read_properties(f'data/test/reads/read{wildcards.r}.fasta')
    total_error_probability = read_properties[4]
    alpha = -math.log(1-total_error_probability)
    C = (2+total_error_probability)/(1-2*alpha)
    return int(C*math.log(target_length,4))

def minimap_k_size(wildcards):
    if config['constant_k']:
        if k > 28:
            return 28
        return k
    read_properties = get_read_properties(f'data/test/reads/read{wildcards.r}.fasta')
    total_error_probability = read_properties[4]
    alpha = -math.log(1-total_error_probability)
    C = (2+total_error_probability)/(1-2*alpha)
    if int(C*math.log(target_length,4)) < 28:
        return int(C*math.log(target_length,4))
    else:
        return 28

#anchor computing
rule benchmark_compute_mem_mummer:
    input: target, query
    output: "data/anchors/mummer-mem/read{r}.anchor"
    params: k = k_size
    shell: "./mummer/mummer -maxmatch -l {params.k} {input[1]} {input[0]} >> {output}"

rule benchmark_compute_mum_mummer:
    input: target, query
    output: "data/anchors/mummer-mum/read{r}.anchor"
    params: k = k_size
    shell: "./mummer/mummer -mum -l {params.k} {input[1]} {input[0]} >> {output}"


rule generate_config_file_for_bdbwt_extended_minimizers:
    input: target, query
    output: "bdbwt-mem/configs/ext_min/config{r}"
    params: k = k_size
    run:
        get_congif_bdbwt_mem(output[0], input[0], input[1], params.k, 1)

rule benchmark_compute_extended_minimizers_bdbwt:
    input: "bdbwt-mem/configs/ext_min/config{r}"
    output: "data/anchors/bdbwt-ext-mini/read{r}.anchor"
    shell: "./bdbwt-mem/main {input} >> {output}"

rule generate_config_file_for_bdbwt_MEM:
    input: target, query
    output: "bdbwt-mem/configs/mem/config{r}"
    params: k = k_size
    run:
        get_congif_bdbwt_mem(output[0], input[0], input[1], params.k, 0)

rule benchmark_compute_MEM_bdbwt:
    input: "bdbwt-mem/configs/mem/config{r}"
    output: "data/anchors/bdbwt-mem/read{r}.anchor"
    shell: "./bdbwt-mem/main {input} >> {output}"

rule benchmark_compute_minimap2_minimizers:
    input: query, target
    output: "data/anchors/minimap/read{r}.anchor"
    params: k = minimap_k_size
    shell: "./minimap2/minimap2 -k {params.k} {input[0]} {input[1]} --print-seeds &>> {output}"

#format the anchors
rule process_mem_mummer_anchors:
    input: "data/anchors/mummer-mem/read{r}.anchor"
    output: "data/anchors/mummer-mem-tidy/read{r}.anchor"
    run:
        f2 = open(output[0], 'w')
        with open(input[0]) as f:
            write2 = False
            for line in f:
                if line[0] != ">":
                    parts = line.split()
                    f2.write(f"{int(parts[1])-1},{int(parts[0])-1},{int(parts[2])}\n")
        f2.close()

use rule process_mem_mummer_anchors as process_mum_mummer_anchors with:
    input: "data/anchors/mummer-mum/read{r}.anchor"
    output: "data/anchors/mummer-mum-tidy/read{r}.anchor"

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
rule benchmark_run_chainX_with_mem_mummer:
    input: target, query, "data/anchors/mummer-mem-tidy/read{r}.anchor"
    benchmark: "benchmarks/chaining/mem_mummer{r}.tsv"
    output: "data/chains/mummer-mem/chain{r}.txt"
    shell: "./ChainX/chainX -m sg -q {input[1]} -t {input[0]} --anchors {input[2]} >> {output}"

use rule benchmark_run_chainX_with_mem_mummer as benchmark_run_chainX_with_mum_mummer with:
    input: target, query, "data/anchors/mummer-mum-tidy/read{r}.anchor"
    benchmark: "benchmarks/chaining/mum_mummer{r}.tsv"
    output: "data/chains/mummer-mum/chain{r}.txt"

use rule benchmark_run_chainX_with_mem_mummer as benchmark_run_chainX_with_extended_minim with:
    input: target, query, "data/anchors/bdbwt-ext-mini-tidy/read{r}.anchor"
    benchmark: "benchmarks/chaining/bdbwt_ext_mini{r}.tsv"
    output: "data/chains/bdbwt-ext-mini/chain{r}.txt"

use rule benchmark_run_chainX_with_mem_mummer as benchmark_run_chainX_with_bdbwt_mem with:
    input: target, query, "data/anchors/bdbwt-mem-tidy/read{r}.anchor"
    benchmark: "benchmarks/chaining/bdbwt_mem{r}.tsv"
    output: "data/chains/bdbwt-mem/chain{r}.txt"

use rule benchmark_run_chainX_with_mem_mummer as benchmark_run_chainX_with__minimizers with:
    input: target, query, "data/anchors/minimap-tidy/read{r}.anchor"
    benchmark: "benchmarks/chaining/minimap_mm{r}.tsv"
    output: "data/chains/minimap/chain{r}.txt"

#anchor stats
rule anchor_stats_mummer_mem:
    input: expand("data/anchors/mummer-mem-tidy/read{r}.anchor", r=r_number)
    output: "results/anchors/mummer-mem.txt"
    params: anchor_type = "mummer-mem"
    run:
        anchor_stats(input, output[0], params.anchor_type)

rule anchor_stats_mummer_mum:
    input: expand("data/anchors/mummer-mum-tidy/read{r}.anchor", r=r_number)
    output: "results/anchors/mummer-mum.txt"
    params: anchor_type = "mummer-mum"
    run:
        anchor_stats(input, output[0], params.anchor_type)

rule anchor_stats_extended_minimizers:
    input: expand("data/anchors/bdbwt-ext-mini-tidy/read{r}.anchor", r=r_number)
    output: "results/anchors/bdbwt-ext-mini.txt"
    params: anchor_type = "bdbwt-ext-mini"
    run:
        anchor_stats(input, output[0], params.anchor_type)

rule anchor_stats_bdbwt_mems:
    input: expand("data/anchors/bdbwt-mem-tidy/read{r}.anchor", r=r_number)
    output: "results/anchors/bdbwt-mem.txt"
    params: anchor_type = "bdbwt-mem"
    run:
        anchor_stats(input, output[0], params.anchor_type)

rule anchor_stats_minimap_mm:
    input: expand("data/anchors/minimap-tidy/read{r}.anchor", r=r_number)
    output: "results/anchors/minimap.txt"
    params: anchor_type = "minimap"
    run:
        anchor_stats(input, output[0], params.anchor_type)

rule results:
    input: expand("data/chains/{anchor_type}/chain{r}.txt", anchor_type=anchor_types, r=r_number)
    output: "results/{anchor_type}-summary.txt"
    params: anchor_type = "{anchor_type}"
    run:
        summary(output[0], params.anchor_type)

def summary(summary_output_path, anchor_type):
        input_folder = f'data/chains/{anchor_type}/'
        chains = get_tuple_list_from_file(input_folder)
        chains_with_empty = get_tuple_list_from_file(input_folder, True)

        reads =  get_reads("data/test/reads/")
        avg_read_length_total = average_length(reads)
        avg_read_length_tidy = average_length({i:reads[i] for i in chains.keys()})

        avg_number_of_anchor_per_chain_total = average_number(chains_with_empty)
        avg_number_of_anchor_per_chain_tidy = average_number(chains)

        avg_number_of_chain_bases_total = average_length(chains_with_empty)
        avg_number_of_chain_bases_tidy = average_length(chains)

        avg_chain_coverage_of_read_total = avg_coverage(chains_with_empty, reads)
        avg_chain_coverage_of_read_tidy = avg_coverage(chains, reads)

        
        avg_jaccard_index_total = average_jaccard_index(chains_with_empty, reads)
        avg_jaccard_index_tidy = average_jaccard_index(chains, reads)

        f2 = open(f"{summary_output_path}", 'w')
        f2.write('total\n')
        f2.write(f'average read length: {avg_read_length_total}\n')
        f2.write(f'average number of anchors per chain: {avg_number_of_anchor_per_chain_total}\n')
        f2.write(f'average number of chain bases: {avg_number_of_chain_bases_total}\n')
        f2.write(f'average chain coverage of read: {avg_chain_coverage_of_read_total}\n')
        f2.write(f'average jaccard index: {avg_jaccard_index_total}\n')
        f2.write('tidy\n')
        f2.write(f'average read length: {avg_read_length_tidy}\n')
        f2.write(f'average number of anchors per chain: {avg_number_of_anchor_per_chain_tidy}\n')
        f2.write(f'average number of chain bases: {avg_number_of_chain_bases_tidy}\n')
        f2.write(f'average chain coverage of read: {avg_chain_coverage_of_read_tidy}\n')
        f2.write(f'average jaccard index: {avg_jaccard_index_tidy}\n')
        f2.close()

def avg_coverage(chains, reads):
    coverage_sum = 0
    for i,chain in chains.items():
        coverage_sum += coverage(chain,1)/reads[i][0][2]
    try:
        return coverage_sum/len(chains)
    except:
        return 0

def average_jaccard_index(chains, reads):
    jaccard_sum=0
    for i, chain in chains.items():
        jaccard_sum += jaccard_index(reads[i][0][2], chain)
    try:
        return jaccard_sum/len(chains)
    except:
        return 0

def jaccard_index(read_length, chain):
    chain_coverage_of_read = coverage(chain, 1)
    total_coverage = coverage(chain, 0) - chain_coverage_of_read + read_length
    return chain_coverage_of_read/total_coverage

def chain_properties(chain):
    chain_length = 0 
    for x,y,l in chain:
        chain_length += l
    return chain_length

def get_read_properties(read_file_path):
    with open(read_file_path) as f:
        read_properties = f.readline().split(';')
        read_length = int(re.findall("\d+", read_properties[1])[0])
        read_start_position = int(re.findall("\d+",read_properties[2])[0])
        read_chromosome = read_properties[3]
        number_of_errors = int(re.findall("\d+",read_properties[4])[0])
        total_error_probability = float(re.findall("\d+.\d+", read_properties[5])[0])
        return read_length, read_start_position, read_chromosome, number_of_errors, total_error_probability

def anchor_stats(input, output, anchor_type):
    anchors_total = get_tuple_list_from_file(f'data/anchors/{anchor_type}-tidy/', True)
    number_of_anchors = sum([len(x) for _,x in anchors_total.items()])
    average_number_of_anchors_per_read = average_number(anchors_total)
    average_length_of_anchors = average_length(anchors_total, True)
    number_of_no_anchor_reads = len([x for _,x in anchors_total.items() if len(x)==0])
    
    anchors_tidy = get_tuple_list_from_file(f'data/anchors/{anchor_type}-tidy/')
    number_of_anchors_tidy = sum([len(x) for _,x in anchors_tidy.items()])
    average_number_of_anchors_per_read_tidy = average_number(anchors_tidy)
    average_length_of_anchors_tidy = average_length(anchors_tidy, True)

    f2 = open(output, 'w')
    f2.write('total:\n')
    f2.write(f'total number of anchors: {number_of_anchors}\n')
    f2.write(f'average number of anchor per read: {average_number_of_anchors_per_read}\n')
    f2.write(f'average base length of anchors: {average_length_of_anchors}\n')
    f2.write(f'number of no anchor reads: {number_of_no_anchor_reads}\n')
    f2.write('tidy:\n')
    f2.write(f'total number of anchors: {number_of_anchors_tidy}\n')
    f2.write(f'average number of anchor per read: {average_number_of_anchors_per_read_tidy}\n')
    f2.write(f'average base length of anchors: {average_length_of_anchors_tidy}\n')
        
    
def get_reads(input_folder):
    dictionary_of_tuple_lists = {}
    for root, dirs, files in os.walk(input_folder):
        for fi in files:
            read_number = re.findall(r"\d+", fi)[0]
            read_properties = get_read_properties(f'{input_folder}{fi}') 
            dictionary_of_tuple_lists[int(read_number)] = [(read_properties[1], 0, read_properties[0])]
    return dictionary_of_tuple_lists

def get_tuple_list_from_file(input_folder, include_empty = False):
    dictionary_of_tuple_lists = {}
    for root, dirs, files in os.walk(input_folder):
        for fi in files:
            read_number = re.findall(r"\d+", fi)[0]
            with open(f'{input_folder}{fi}') as f:
                tuple_list = []
                for line in f:
                    try:
                        parts = re.findall(r"\d+", line)
                        x = int(parts[0])
                        y = int(parts[1])
                        length = int(parts[2])
                        tuple_list.append((x,y,length))
                    except:
                        continue
                if  len(tuple_list)>0 and not include_empty:
                    dictionary_of_tuple_lists[int(read_number)] = tuple_list
                elif include_empty:
                    dictionary_of_tuple_lists[int(read_number)] = tuple_list
    return dictionary_of_tuple_lists

#return average length of tuple list (a,b,l)
#average length of reads, average length of anchors
def average_length(input_dic, anchors=False):
    sum_of_chain_lengths = 0
    for _,chain in input_dic.items():
        sum_of_chain_lengths += sum([l for a,b,l in chain])
    try:
        if anchors:
            return sum_of_chain_lengths/sum([len(x) for _,x in input_dic.items()])
        else:  
            return sum_of_chain_lengths/len(input_dic)
    except:
        return 0

def average_number(input_dic):
    try:
        return sum([len(x) for _,x in input_dic.items()]) / len(input_dic)
    except:
        return 0

# a position in target, b position in query, get coverage of target sequence, assume that index is sorted by the start position
# if last end position < cur start position, pop the anchor
def coverage(input_list, index=0):
    if len(input_list) == 0:
        return 0
    if index==0:
        input_list.sort(key=lambda x: x[0])
    else:
        input_list.sort(key=lambda x: x[1])
    previous_positions = deque([input_list[0]])

    coverage = input_list[0][2]
    for a,b,l in input_list[1:]:
        
        if index == 0:
            cur_pos = a
        else:
            cur_pos = b

        pre_a, pre_b, pre_l = previous_positions[-1]
        if index == 0:
            pre_pos = pre_a
        else:
            pre_pos = pre_b

        if pre_pos+pre_l > cur_pos+l:
            continue
        if pre_pos+pre_l-1 > cur_pos:
            overlap_size = pre_pos + pre_l - cur_pos
            coverage += l - overlap_size
            previous_positions.append((a,b,l))
        else:
            coverage += l
        #remove all that are not in range
        while len(previous_positions)>0:
            a2, b2, l2 = previous_positions[0]
            if index == 0:
                pos2 = a2
            else:
                pos2 = b2
            if pos2 + l2-1 <= cur_pos:
                previous_positions.popleft()
            else:
                break
        if len(previous_positions) == 0:
            previous_positions.append((a,b,l))
        

    return coverage


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
