import math
import argparse
import subprocess
from constant import *
import time
import csv


RUN_PATH = ''


def generate_congif_bdbwt_mem(output_file_path, target_file_path, query_file_path, k, mode) -> str:
    # generate a text file for bdbwt mem tool
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


def init_variable_k(target, read_properties) -> int:
    # generates the k variable according to Shawn, Yu
    print('generating variable k')
    target_length = len(target)
    total_error_probability = read_properties[4]
    alpha = -math.log(1-total_error_probability)
    C = (2+total_error_probability)/(1-2*alpha)
    k = int(C*math.log(target_length, 4))
    return k, k if k < 28 else 28


def anchor_stats(mode, id, anchors, k_value):
    precision_val = precision(anchors, get_read(id))
    return id, len(anchors), precision_val, mode, k_value


def precision(anchors, read):
    try:
        return true_positive_anchors(anchors, read)/len(anchors)
    except:
        return 0


def true_positive_anchors(anchors, read):
    read_start_pos = read[0]
    read_length = read[2]
    true_positives = 0
    for (t, _, l) in anchors:
        if read_start_pos <= t and t+l <= read_start_pos+read_length:
            true_positives += 1
    return true_positives


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
    elif mode == BR_INDEX_MEM:
        if not os.path.isfile(target_path + '.bri'):
            subprocess.run(R_INDEX_ALGO.format(
                target_path, RUN_PATH), shell=True)
        subprocess.run(R_INDEX_ALGO.format(query_path, RUN_PATH), shell=True)
        start_time = time.time()
        subprocess.run(ANCHOR_ALGOS[mode].format(
            k_value, output_path, target_path+'.bri', query_path+'.bri', RUN_PATH), shell=True)
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


def parse_anchors(mode, input_path, output_path='', read_str='', target_str=''):
    # returns a list of tuples if output path is not given
    print('parsing anchors')
    anchors = []
    with open(input_path) as f:
        write = False
        for line in f:
            if mode == BDBWT_EXT_MINI or mode == BDBWT_MEM:
                if write:
                    parts = line.split(',')
                    anchors.append(
                        ((int(parts[1]), int(parts[0]), int(parts[2]))))
                if line.strip() == "MEMs:":
                    write = True
            if mode == MUMMER_MEM or mode == MUMMER_MUM:
                if line[0] != ">":
                    parts = line.split()
                    anchors.append(
                        (int(parts[1])-1, int(parts[0])-1, int(parts[2])))
            if mode == MINIMAP:
                if line[0:2] == 'SD':
                    parts = line.split()
                    x = int(parts[2])
                    y = int(parts[4])
                    k = int(parts[5])
                    if read_str[x-k+1:x] == target_str[y-k+1:y]:
                        anchors.append((y-k+1, x-k+1, k))
            if mode == BDBWT_MUM or mode == EXTENDED_MINIMAP:
                parts = line.split(',')
                anchors.append((int(parts[0]), int(parts[1]), int(parts[2])))
        if output_path == '':
            return anchors
        write_anchors(anchors, output_path)


def write_anchors(anchors, output_path):
    writer = open(output_path, 'w')
    for x, y, l in anchors:
        writer.write(f'{x},{y},{l}\n')
    writer.close()


def generate_ext_minimizers(anchors, target, read, output_path):
    print('generate extended minimizers')
    new_anchors = []
    for x, y, length in anchors:
        new_x = x
        new_y = y
        # try extend left
        while True:
            if target[new_x-1:x+length] == read[new_y-1:y+length]:
                # match
                new_x -= 1
                new_y -= 1
            else:
                break
        # try extend right
        new_length = length
        while True:
            if target[new_x:x+new_length+1] == read[new_y:y+new_length+1]:
                new_length += 1
            else:
                break
        new_anchors.append((new_x, new_y, new_length))
    write_anchors(new_anchors, output_path)


def generate_MUMs(anchors, output_path, read):
    print('generating MUMs')
    anchor_counter = {}
    for x, y, l in anchors:
        if read[y:y+l] in anchor_counter.keys():
            anchor_counter[read[y:y+l]
                           ] = (anchor_counter[read[y:y+l]][0]+1, (x, y, l))
        else:
            anchor_counter[read[y:y+l]] = (1, (x, y, l))
    new_anchors = []
    for _, (n, position) in anchor_counter.items():
        if n == 1:
            new_anchors.append(position)
    if output_path == '':
        return new_anchors

    write_anchors(new_anchors, output_path)
    

def main(target_path, k, read_id, test_mode):
    target_str = get_target(target_path)
    read_str = get_read_str(READ_PATH.format(read_id))
    read_properties = get_read_properties(READ_PATH.format(read_id))

    # generate k i variable k
    if k > 0:
        print('costant k')
        k_value = k
        if k > 28:
            minimap_k_value = 28
        else:
            minimap_k_value = k
    else:
        k_value, minimap_k_value = init_variable_k(target_str, read_properties)

    # run anchoring algos
    run_anchoring_algo(BDBWT_EXT_MINI, read_id, k_value, target_path, READ_PATH.format(
        read_id), ANCHOR_PATH[BDBWT_EXT_MINI].format(read_id), True, CONFIG_PATH[BDBWT_EXT_MINI].format(read_id))
    run_anchoring_algo(BDBWT_MEM, read_id, k_value, target_path, READ_PATH.format(
        read_id), ANCHOR_PATH[BDBWT_MEM].format(read_id), True, CONFIG_PATH[BDBWT_MEM].format(read_id))
    run_anchoring_algo(MUMMER_MEM, read_id, k_value, target_path,
                       READ_PATH.format(read_id), ANCHOR_PATH[MUMMER_MEM].format(read_id), True)
    run_anchoring_algo(MUMMER_MUM, read_id, k_value, target_path,
                       READ_PATH.format(read_id), ANCHOR_PATH[MUMMER_MUM].format(read_id), True)
    run_anchoring_algo(MINIMAP, read_id, minimap_k_value, target_path,
                       READ_PATH.format(read_id), ANCHOR_PATH[MINIMAP].format(read_id), True)
    run_anchoring_algo(BR_INDEX_MEM, read_id, k_value, target_path, READ_PATH.format(
        read_id), TIDY_ANCHOR_PATH[BR_INDEX_MEM].format(read_id), True)
    # parse anchors
    for mode, _ in ANCHOR_ALGOS.items():
        if mode == MINIMAP:
            parse_anchors(mode, input_path=ANCHOR_PATH[mode].format(
                read_id), output_path=TIDY_ANCHOR_PATH[mode].format(read_id), read_str=read_str, target_str=target_str)
        elif mode != BR_INDEX_MEM:
            parse_anchors(mode, input_path=ANCHOR_PATH[mode].format(
                read_id), output_path=TIDY_ANCHOR_PATH[mode].format(read_id))

    # generate the custom anchors
    generate_MUMs(get_tuple_list(TIDY_ANCHOR_PATH[BDBWT_MEM].format(
        read_id)), TIDY_ANCHOR_PATH[BDBWT_MUM].format(read_id), read_str)
    generate_MUMs(get_tuple_list(TIDY_ANCHOR_PATH[BR_INDEX_MEM].format(
        read_id)), TIDY_ANCHOR_PATH[BR_INDEX_MUM].format(read_id), read_str)
    generate_ext_minimizers(get_tuple_list(TIDY_ANCHOR_PATH[MINIMAP].format(
        read_id)), target_str, read_str, TIDY_ANCHOR_PATH[EXTENDED_MINIMAP].format(read_id))

    # compute anchor stats
    if not os.path.isfile(ANCHOR_STATS_PATH):
        with open(ANCHOR_STATS_PATH, 'a', newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow(
                ['read_number', 'number_of_anchors', 'precision', 'mode', 'k_value'])

    for mode, path in TIDY_ANCHOR_PATH.items():
        stats = anchor_stats(mode, read_id,
                             get_tuple_list(path.format(read_id)), k_value)
        with open(ANCHOR_STATS_PATH, 'a', newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow(stats)

    # check the anchor correctness if test mode
    if test_mode:
        for mode in ANCHOR_TYPES:
            anchors = get_tuple_list(TIDY_ANCHOR_PATH[mode].format(read_id))
            for a, b, l in anchors:
                if target_str[a:a+l] != read_str[b:b+l]:
                    print(mode)
                    print(mode)
                    print('mismatch in anchors')
                    print(f'target:{target_str[a:a+l]}')
                    print(f'read: {read_str[b:b+l]}')

    run_chainx(read_id, target_path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--target", help="path to target file")
    parser.add_argument("-tm", "--test_mode",
                        help="use test mode", action='store_true')
    parser.add_argument("-i", "--read_id", help="id of the input read")
    parser.add_argument("-k", "--k-size",
                        help="Size of constant k.", default=-1, type=int)
    parser.add_argument("-p", "--path",
                        help="path of run files", default='')
    parser.add_argument(
        "--read_range", help='run on multiple reads', default=-1, type=int)
    
    args = parser.parse_args()
    target = args.target
    read_id = args.read_id
    k_size = args.k_size
    RUN_PATH = args.path
    read_range = args.read_range
    test_mode = args.test_mode


    
    if not os.path.isfile(f'{RESULT_FOLDER}info.txt'):
        f2 = open(f'{RESULT_FOLDER}info.txt', 'w')
        f2.write(f'k:{k_size}\n')
        f2.write(f'genome: {target.split(".")[0]}')
        f2.close()
        
    if read_range > 0:
        for i in range(read_range):
            main(target, k_size, f'{i}', test_mode)
    else:
        main(target, k_size, read_id, test_mode)
