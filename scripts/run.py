import re
import math
import os
import argparse
import subprocess
from constant import *
import time


RUN_PATH = ''


def generate_congif_bdbwt_mem(output_file_path, target_file_path, query_file_path, k, mode):
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
    f2.write("BWTThrd > 0\n")
    f2.write("VerbCA  > 0\n")
    f2.write("VerbED > 0\n")
    f2.write("RawChain > 0\n")
    f2.write("strChain > 0\n")
    f2.write("recombAbs > 0\n")
    f2.write("linearRMQ > 0")
    f2.close()


def init_variable_k(target, read_properties):
    print('generating variable k')
    target_length = len(target)
    total_error_probability = read_properties[4]
    alpha = -math.log(1-total_error_probability)
    C = (2+total_error_probability)/(1-2*alpha)
    k = int(C*math.log(target_length, 4))
    return k, k if k < 28 else 28


def anchor_stats(k, mode, target, id, anchors):
    precision_val = precision(anchors, get_read(id))
    recall_val = recall(anchors, get_read(id), target)
    return


def recall(anchors, read, target, target_path, mode, k_value, read_id):
    # run the anchoring algorithm for the read and the substring of read in the query
    read_start_pos = read[0]
    read_end_pos = read_start_pos + read[2]
    temp = target[read_start_pos:read_end_pos]
    writer = open(ANCHOR_STATS_PATH[mode].format(read_id), 'w')
    writer.write(">testing \n")
    writer.write(temp)
    writer.close()
    if mode == BDBWT_EXT_MINI or mode == BDBWT_MEM:
        run_anchoring_algo(mode, read_id, k_value, target_path, ANCHOR_STATS_PATH[mode].format(
            read_id), ANCHOR_STATS_PATH[mode].format(read_id), False, config_path=ANCHOR_STATS_PATH[mode].format(f'config{read_id}'))
    else:
        run_anchoring_algo(mode, read_id, k_value, target_path, ANCHOR_STATS_PATH[mode].format(
            read_id), ANCHOR_STATS_PATH[mode].format(read_id))


# list of tuples
def precision(anchors, read):
    return len(true_positive_anchors(anchors, read))/len(anchors)


def true_positive_anchors(anchors, read):
    read_start_pos = read[0]
    read_length = read[2]
    true_positives = []
    for (t, q, l) in anchors:
        if read_start_pos <= t and t+l <= read_start_pos+read_length:
            true_positives += (t, q, l)
    return true_positives


def get_read_properties(id):
    with open(READ_PATH.format(id)) as f:
        for line in f:
            if line[0] == ">":
                read_properties = line.split(';')
                read_length = int(re.findall("\d+", read_properties[1])[0])
                read_start_position = int(
                    re.findall("\d+", read_properties[2])[0])
                read_chromosome = read_properties[3]
                number_of_errors = int(re.findall(
                    "\d+", read_properties[4])[0])
                total_error_probability = float(
                    re.findall("\d+.\d+", read_properties[5])[0])
                read_properties = read_length, read_start_position, read_chromosome, number_of_errors, total_error_probability
                return read_properties


def save_time(start_time, end_time, path):
    writer = open(path, 'w')
    writer.write("--- %s seconds ---" % (end_time - start_time))
    writer.close()


def run_anchoring_algo(mode, read_id, k_value, target_path, query_path, output_path, measure_time=False, config_path=''):
    if os.path.isfile(output_path):
        return
    if config_path != '' and os.path.isfile(config_path):
        return
    print(f'running {mode} for read {read_id}')

    if mode == BDBWT_EXT_MINI or mode == BDBWT_MEM:
        generate_congif_bdbwt_mem(
            config_path, target_path, query_path, k_value, 0 if mode == BDBWT_MEM else 1)
        start_time = time.time()
        subprocess.run(ANCHOR_ALGOS[mode].format(
            config_path, output_path, RUN_PATH), shell=True)
        end_time = time.time()
    else:
        start_time = time.time()
        subprocess.run(ANCHOR_ALGOS[mode].format(
            k_value, query_path, target_path, output_path, RUN_PATH), shell=True)
        end_time = time.time()
    if measure_time:
        save_time(start_time, end_time,
                  BENCHMARK_ANCHOR_PATH[mode].format(read_id))


def run_chainx(read_id, target_path):
    print('running chainX')
    for type_of, path in CHAIN_PATH.items():
        i = read_id
        if not os.path.isfile(path.format(i)):
            start_time = time.time()
            subprocess.run(
                f"./{RUN_PATH}ChainX/chainX -m sg -q {READ_PATH.format(i)} -t {target_path} --anchors {TIDY_ANCHOR_PATH[type_of].format(i)} >> {path.format(i)}", shell=True)
            end_time = time.time()
            save_time(start_time, end_time,
                      BENCHMARK_CHAIN_PATH[type_of].format(read_id))


def parse_anchors(read_id):
    print('parsing anchors')
    i = read_id
    for type_of, path in ANCHOR_PATH.items():
        writer = open(TIDY_ANCHOR_PATH[type_of].format(i), 'w')
        with open(path.format(i)) as f:
            write = False
            for line in f:
                if type_of == BDBWT_EXT_MINI or type_of == BDBWT_MEM:
                    if write:
                        parts = line.split(',')
                        writer.write(
                            f"{int(parts[1])},{int(parts[0])},{int(parts[2])}\n")
                    if line.strip() == "MEMs:":
                        write = True
                if type_of == MUMMER_MEM or type_of == MUMMER_MUM:
                    if line[0] != ">":
                        parts = line.split()
                        writer.write(
                            f"{int(parts[1])-1},{int(parts[0])-1},{int(parts[2])}\n")
                if type_of == MINIMAP:
                    if line[0:2] == 'SD':
                        parts = line.split()
                        x = int(parts[2])
                        y = int(parts[4])
                        k = int(parts[5])
                        writer.write(f"{y-k+1},{x-k+1},{k}\n")
            writer.close()


def main(target, k, read_id):
    target_str = get_target(target)
    read_properties = get_read_properties(read_id)
    if k > 0:
        print('costant k')
        k_value = k
        if k > 28:
            minimap_k_value = 28
        else:
            minimap_k_value = k
    else:
        k_value, minimap_k_value = init_variable_k(target_str, read_properties)

    run_anchoring_algo(BDBWT_EXT_MINI, read_id, k_value, target, READ_PATH.format(
        read_id), ANCHOR_PATH[BDBWT_EXT_MINI].format(read_id), True, CONFIG_PATH[BDBWT_EXT_MINI].format(read_id))
    run_anchoring_algo(BDBWT_MEM, read_id, k_value, target, READ_PATH.format(
        read_id), ANCHOR_PATH[BDBWT_MEM].format(read_id), True, CONFIG_PATH[BDBWT_MEM].format(read_id))
    run_anchoring_algo(MUMMER_MEM, read_id, k_value, target,
                       READ_PATH.format(read_id), ANCHOR_PATH[MUMMER_MEM], True)
    run_anchoring_algo(MUMMER_MUM, read_id, k_value, target,
                       READ_PATH.format(read_id), ANCHOR_PATH[MUMMER_MUM], True)
    run_anchoring_algo(MINIMAP, read_id, k_value, target,
                       READ_PATH.format(read_id), ANCHOR_PATH[MINIMAP], True)

    parse_anchors(read_id)

    run_chainx(read_id, target)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--target", help="path to target file")
    parser.add_argument("-i", "--read_id", help="id of the input read")
    parser.add_argument("-k", "--k-size",
                        help="Size of constant k.", default=-1, type=int)
    parser.add_argument("-p", "--path",
                        help="path of run files", default='')

    args = parser.parse_args()
    target = args.target
    read_id = args.read_id
    k_size = args.k_size
    RUN_PATH = args.path

    main(target, k_size, read_id)
