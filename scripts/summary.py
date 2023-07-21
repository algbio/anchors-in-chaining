from constant import *
from collections import deque
import csv
import shutil


def runtime_summary(input_folder):
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
    return avg_time


def avg_coverage(chains, reads):
    coverage_sum = 0
    for i, chain in chains.items():
        coverage_sum += coverage(chain, 1)/reads[i][0][2]
    try:
        return coverage_sum/len(chains)
    except:
        return 0


def average_jaccard_index(chains, reads):
    jaccard_sum = 0
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
    for x, y, l in chain:
        chain_length += l
    return chain_length


def anchor_stats(anchor_type, include_empty):
    anchors = get_tuple_list_from_file(
        DIRS['anchor-tidy'][anchor_type], include_empty)

    number_of_anchors = sum([len(x) for _, x in anchors.items()])
    average_number_of_anchors_per_read = average_number(anchors)
    average_length_of_anchors = average_length(anchors, True)
    run_time = runtime_summary(DIRS['benchmarks-anchors'][anchor_type])
    return [anchor_type, 'total' if include_empty else 'tidy', number_of_anchors, average_number_of_anchors_per_read, average_length_of_anchors, run_time]


# return average length of tuple list (a,b,l)
# average length of reads, average length of anchors
def average_length(input_dic, anchors=False):
    sum_of_lengths = 0
    for _, chain in input_dic.items():
        sum_of_lengths += sum([l for a, b, l in chain])
    try:
        if anchors:
            return sum_of_lengths/sum([len(x) for _, x in input_dic.items()])
        else:
            return sum_of_lengths/len(input_dic)
    except:
        return 0


def average_number(input_dic):
    try:
        return sum([len(x) for _, x in input_dic.items()]) / len(input_dic)
    except:
        return 0

# a position in target, b position in query, get coverage of target sequence, assume that index is sorted by the start position
# if last end position < cur start position, pop the anchor


def coverage(input_list, index=0):
    if len(input_list) == 0:
        return 0
    if index == 0:
        input_list.sort(key=lambda x: x[0])
    else:
        input_list.sort(key=lambda x: x[1])
    previous_positions = deque([input_list[0]])

    coverage = input_list[0][2]
    for a, b, l in input_list[1:]:

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
            previous_positions.append((a, b, l))
        else:
            coverage += l
        # remove all that are not in range
        while len(previous_positions) > 0:
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
            previous_positions.append((a, b, l))

    return coverage


def results(type_of, include_empty_anchors, reads):
    input_folder = DIRS['chain'][type_of]
    chains = get_tuple_list_from_file(input_folder, include_empty_anchors)
    avg_read_length = average_length(
        {i: reads[i] for i in chains.keys()})
    avg_number_of_anchors_per_chain = average_number(chains)
    avg_number_of_chain_bases = average_length(chains)
    avg_chain_coverage_of_read = avg_coverage(chains, reads)
    avg_jaccard_index = average_jaccard_index(chains, reads)
    avg_runtime = runtime_summary(DIRS['benchmarks-chains'][type_of])
    return [type_of, 'total' if include_empty_anchors else 'tidy', avg_read_length, avg_number_of_anchors_per_chain, avg_number_of_chain_bases, avg_chain_coverage_of_read, avg_jaccard_index, avg_runtime]


def main():
    if os.path.isfile(f'{RESULT_FOLDER}info.txt'):
        with open(f'{RESULT_FOLDER}info.txt') as f:
            for line in f:
                if line[0] == 'k':
                    k_int = int(re.findall(r"-*\d+", line)[0])
                    k = 'var' if k_int < 0 else k_int
                if line[:len('genome')] == 'genome':
                    genome = line.split(':')[1].split('/')[-1]
        os.remove(f'{RESULT_FOLDER}info.txt')
    chain_summary_values = [['type', 'mode', 'avg_read_length', 'avg_number_of_anchors_per_chain',
                             'avg_number_of_chain_bases', 'avg_chain_coverage_of_read', 'avg_jaccard_index', 'avg_runtime']]
    reads = get_reads(READS_DIR)
    anchor_summary_values = [['type', 'mode', 'number_of_anchors',
                              'average_number_of_anchors_per_read', 'average_length_of_anchors', 'run_time']]
    for type_of in ANCHOR_ALGOS.keys():
        chain_summary_values.append(results(type_of, True, reads))
        chain_summary_values.append(results(type_of, False, reads))
        anchor_summary_values.append(anchor_stats(type_of, True))
        anchor_summary_values.append(anchor_stats(type_of, False))
    with open(CHAIN_SUMMARY_PATH.format(k=k, genome=genome), 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerows(chain_summary_values)
    with open(ANCHOR_SUMMARY_PATH.format(k=k, genome=genome), 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerows(anchor_summary_values)


if __name__ == '__main__':

    main()
