import os
import re
from constant import *
from collections import deque

def summary(type_of):
        input_folder = DIRS['chain'][type_of]
        chains = get_tuple_list_from_file(input_folder)
        chains_with_empty = get_tuple_list_from_file(input_folder, True)
        reads =  get_reads(READS_DIR)
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

        f2 = open(SUMMARY_PATH.format(type_of), 'a')
        f2.write(f'{type_of}\n')
        f2.write('chain summary \n')
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
        runtime_summary(DIRS['benchmarks-chains'][type_of], type_of)

def runtime_summary(input_folder, type_of):
    times = []
    for _, _, files in os.walk(input_folder):
        for fi in files:
            try:
                read_number = re.findall(r"\d+", fi)[0]
            except:
                continue
            with open(f'{input_folder}read{read_number}.txt') as f:
                time = re.findall(r"\d+.\d+", f.readline())[0]
            times.append(float(time))
    avg_time = sum(times)/len(times)
    writer =  open(SUMMARY_PATH.format(type_of), 'a')
    writer.write(f'average runtime {avg_time} seconds\n')
    writer.close()

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

def anchor_stats(anchor_type):
    anchors_total = get_tuple_list_from_file(DIRS['anchor-tidy'][anchor_type], True)

    number_of_anchors = sum([len(x) for _,x in anchors_total.items()])
    average_number_of_anchors_per_read = average_number(anchors_total)
    average_length_of_anchors = average_length(anchors_total, True)
    anchors_tidy = get_tuple_list_from_file(DIRS['anchor-tidy'][anchor_type])
    number_of_anchors_tidy = sum([len(x) for _,x in anchors_tidy.items()])
    average_number_of_anchors_per_read_tidy = average_number(anchors_tidy)
    average_length_of_anchors_tidy = average_length(anchors_tidy, True)
    f22 = open(SUMMARY_PATH.format(anchor_type), 'a')
    f22.write('anchor summary\n')
    f22.write('total:\n')
    f22.write(f'total number of anchors: {number_of_anchors}\n')
    f22.write(f'average number of anchor per read: {average_number_of_anchors_per_read}\n')
    f22.write(f'average base length of anchors: {average_length_of_anchors}\n')
    f22.write(f'number of no anchor reads: {len(anchors_total)-len(anchors_tidy)}\n')
    f22.write('tidy:\n')
    f22.write(f'total number of anchors: {number_of_anchors_tidy}\n')
    f22.write(f'average number of anchor per read: {average_number_of_anchors_per_read_tidy}\n')
    f22.write(f'average base length of anchors: {average_length_of_anchors_tidy}\n')
    f22.close()
    runtime_summary(DIRS['benchmarks-anchors'][anchor_type], anchor_type)

def get_reads(input_folder):
    dictionary_of_tuple_lists = {}
    for _, _, files in os.walk(input_folder):
        for fi in files:
            try:
                read_number = re.findall(r"\d+", fi)[0]
            except:
                continue
            read_properties = get_read_properties(f'{input_folder}{fi}') 
            dictionary_of_tuple_lists[int(read_number)] = [(read_properties[1], 0, read_properties[0])]
    return dictionary_of_tuple_lists

def get_tuple_list_from_file(input_folder, include_empty = False):
    dictionary_of_tuple_lists = {}
    for _, _, files in os.walk(input_folder):
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
    sum_of_lengths = 0
    for _,chain in input_dic.items():
        sum_of_lengths += sum([l for a,b,l in chain])
    try:
        if anchors:
            return sum_of_lengths/sum([len(x) for _,x in input_dic.items()])
        else:  
            return sum_of_lengths/len(input_dic)
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

def main():
    for type_of in ANCHOR_ALGOS.keys():
        summary(type_of)
        anchor_stats(type_of)

if __name__ == '__main__':

    main()