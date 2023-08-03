#!/usr/bin/python3
from constant import *
from collections import deque
import csv
import pandas as pd


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
    try:
        avg_time = sum(times)/len(times)
    except:
        avg_time = 0
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
    # print(chains)
    jaccard_sum = 0
    for i, chain in chains.items():
        value = jaccard_index(reads[i][0][0], reads[i][0][2], chain)
        jaccard_sum += value
        print(value)

    print(jaccard_sum)
    try:
        return jaccard_sum/len(chains)
    except:
        return 0


def jaccard_index(read_start_pos, read_length, chain):
    # print(read_start_pos)
    # print(read_length)
    # print(chain)
    intersection_val = intersection(chain, read_start_pos, read_length)

    #print(intersection_val)
    chain_length = coverage(chain, 1)
    # print(chain_length)
    # print(read_length+1)
    return intersection_val/(read_length+1 + chain_length - intersection_val)


def chain_properties(chain):
    chain_length = 0
    for x, y, l in chain:
        chain_length += l
    return chain_length


def anchor_stats(anchor_type, include_empty):
    anchors = get_tuple_list_from_file(
        DIRS['anchor-tidy'][anchor_type], include_empty)

    number_of_anchors = sum([len(x) for _, x in anchors.items()])
    average_length_of_anchors = average_length(anchors, True)
    run_time = runtime_summary(DIRS['benchmarks-anchors'][anchor_type])
    if include_empty:
        df = pd.read_csv(ANCHOR_STATS_PATH)
    else:
        df = pd.read_csv(ANCHOR_STATS_PATH)
        df = df[df['number_of_anchors'] > 0]
    averages = df.groupby('mode').mean()

    return [anchor_type,  # type
            'total' if include_empty else 'tidy',  # mode
            number_of_anchors,  # number_of_anchors
            # average_number_of_anchors_per_read
            averages['number_of_anchors'][anchor_type],
            average_length_of_anchors,  # average_length_of_anchors
            averages['precision'][anchor_type],  # average_precision
            df[(df['mode'] == anchor_type) & (df['number_of_anchors'] == 0)
               ].shape[0],  # 'number_of_reads_without_anchors',
            run_time,  # run_time
            averages['k_value'][anchor_type]  # 'k_value'
            ]


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


def intersection(chain, read_start_pos, length) -> int:
    # compute an intersection of chain and string
    # a position in target, b position in query,
    read_end_pos = read_start_pos + length
    if len(chain) == 0:
        return 0
    chain.sort(key=lambda x: x[0])
    intersection_size = 0
    previous_positions = deque()

    for a, b, l in chain:
        if a + l-1 < read_start_pos:
            continue
        if a > read_end_pos:
            break
        offset = 0
        if a < read_start_pos:
            offset += abs(a-read_start_pos)
        if a+l-1 > read_end_pos:
            offset += abs(a+l-1 - read_end_pos)

        if len(previous_positions) == 0:
            previous_positions.append((a, b, l))
            intersection_size += l - offset
            continue

        cur_pos = a
        pre_a, pre_b, pre_l = previous_positions[-1]
        pre_pos = pre_a

        if pre_pos+pre_l > cur_pos+l:
            continue
        if pre_pos+pre_l-1 > cur_pos:
            # overlap
            overlap_size = pre_pos + pre_l - cur_pos
            intersection_size += l - overlap_size - offset
            previous_positions.append((a, b, l))
        else:
            # no overlap
            intersection_size += l - offset
        # remove all that are not in range
        while len(previous_positions) > 0:
            a2, b2, l2 = previous_positions[0]
            pos2 = a2
            if pos2 + l2-1 <= cur_pos:
                previous_positions.popleft()
            else:
                break
        if len(previous_positions) == 0:
            previous_positions.append((a, b, l))
    print(intersection_size)
    return intersection_size


def coverage(input_list, index=0) -> int:
    # get coverage (intersection) of target sequence,
    # assume that index is sorted by the start position
    # if last end position < cur start position, pop the anchor
    if len(input_list) == 0:
        return 0
    if index == 0:
        input_list.sort(key=lambda x: x[0])
    else:
        input_list.sort(key=lambda x: x[1])
    previous_positions = deque([input_list[0]])

    coverage = input_list[0][2]
    for a, b, l in input_list[1:]:
        # print((a,b,l))
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
    else:
        genome = 'test'
        k = 'test'
    chain_summary_values = [['type',
                             'mode',
                             'avg_read_length',
                             'avg_number_of_anchors_per_chain',
                             'avg_number_of_chain_bases',
                             'avg_chain_coverage_of_read',
                             'avg_jaccard_index',
                             'avg_runtime']]
    reads = get_reads(READS_DIR)
    anchor_summary_values = [['type',
                              'mode',
                              'number_of_anchors',
                              'average_number_of_anchors_per_read',
                              'average_length_of_anchors',
                              'average_precision',
                              'number_of_reads_without_anchors',
                              'run_time',
                              'k_value']]
    for type_of in ANCHOR_TYPES:
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
