#!/usr/bin/python3
import re
import math
from collections import deque
import os
import argparse
import subprocess
import shutil

BDBWT_MEM = 'bdbwt-mem'
BDBWT_EXT_MINI = 'bdbwt-ext-mini'
MUMMER_MUM = 'mummer-mum'
MUMMER_MEM = 'mummer_mem'
MINIMAP = 'minimap'

DATA_FOLDER = "data"
ANCHOR_DIR = f'{DATA_FOLDER}/anchors/'
TIDY_ANCHOR_DIR = f'{DATA_FOLDER}/anchors-tidy/'
CHAIN_DIR = f'{DATA_FOLDER}/chains/'
READS_DIR = f'{DATA_FOLDER}/reads/'
RESULT_FOLDER = 'results/'
ANCHOR_ALGOS = {BDBWT_MEM: "./bdbwt-mem/main {0} >> {1}",
                BDBWT_EXT_MINI: "./bdbwt-mem/main {0} >> {1}",
                MUMMER_MUM: "./mummer/mummer -mum -l {0} {1} {2} >> {3}",
                MUMMER_MEM: "./mummer/mummer -maxmatch -l {0} {1} {2} >> {3}",
                MINIMAP: "./minimap2/minimap2 -k {0} {1} {2} --print-seeds > {3} 2>&1"}
CONFIG_DIR_MEM = f'bdbwt-mem/configs/{BDBWT_MEM}/'
CONFIG_DIR_EXT_MINI = f'bdbwt-mem/configs/{BDBWT_EXT_MINI}/'
DIRS = {'read': READS_DIR,
        f'config {BDBWT_MEM}': CONFIG_DIR_MEM,
        f'config {BDBWT_EXT_MINI}': CONFIG_DIR_EXT_MINI,
        'chain': {x: f'{CHAIN_DIR}{x}/' for x in ANCHOR_ALGOS.keys()},
        'anchor': {x: f'{ANCHOR_DIR}{x}/' for x in ANCHOR_ALGOS.keys()},
        'anchor-tidy': {x: f'{ANCHOR_DIR}{x}/' for x in ANCHOR_ALGOS.keys()}}

READ_PATH = f'{READS_DIR}'+'read{0}.fasta'
CONFIG_PATH = {BDBWT_MEM: DIRS[f"config {BDBWT_MEM}"]+'config{0}',
               BDBWT_EXT_MINI: DIRS[f"config {BDBWT_EXT_MINI}"]+'config{0}'}
ANCHOR_PATH = {anchor_type: DIRS['anchor'][anchor_type] +
               'read{0}.txt' for anchor_type in ANCHOR_ALGOS.keys()}
TIDY_ANCHOR_PATH = {anchor_type: DIRS['anchor-tidy'][anchor_type] +
                    'read{0}.txt' for anchor_type in ANCHOR_ALGOS.keys()}
CHAIN_PATH = {anchor_type: DIRS['chain'][anchor_type] +
              'read{0}.txt' for anchor_type in ANCHOR_ALGOS.keys()}

def fresh_run():
    for type_of, dir in DIRS.items():
        # print(type_of)
        # print(dir)
        if type_of == 'chain' or type_of == 'anchor' or type_of == 'anchor-tidy':
            for _, sub_dir in dir.items():
                shutil.rmtree(sub_dir, ignore_errors=True)
        else:
            shutil.rmtree(dir, ignore_errors=True)

def init_directories():
    for type_of, dir in DIRS.items():
        # print(type_of)
        # print(dir)
        if type_of == 'chain' or type_of == 'anchor' or type_of == 'anchor-tidy':
            for _, sub_dir in dir.items():
                if not os.path.exists(sub_dir):
                    os.makedirs(sub_dir)
        else:
            if not os.path.exists(dir):
                os.makedirs(dir)


def generate_congif_bdbwt_mem(output_file_path, target_file_path, query_file_path, k, mode, ):
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


def init_variable_k(target, read_properties_dic):
    print('generating variable k')
    variable_k = []
    variable_k_minimap = []
    target_length = 0
    with open(target) as f:
        for j, line in enumerate(f):
            if j != 0:
                target_length += len(line.strip())
    for _,read_properties in read_properties_dic.items():
        total_error_probability = read_properties[4]
        alpha = -math.log(1-total_error_probability)
        C = (2+total_error_probability)/(1-2*alpha)
        k = int(C*math.log(target_length, 4))

        variable_k.append(k)

        if k < 28:
            variable_k_minimap.append(k)
        else:
            variable_k_minimap.append(28)

    return variable_k, variable_k_minimap


def get_read_properties(fastq_file_path):
    dic_read_properties = {}
    i=0
    with open(fastq_file_path) as f:
        for line in f:
            if line[0] == "@":
                read_properties = line.split(';')
                read_length = int(re.findall("\d+", read_properties[1])[0])
                read_start_position = int(re.findall("\d+", read_properties[2])[0])
                read_chromosome = read_properties[3]
                number_of_errors = int(re.findall("\d+", read_properties[4])[0])
                total_error_probability = float(re.findall("\d+.\d+", read_properties[5])[0])
                dic_read_properties[i] = read_length, read_start_position, read_chromosome, number_of_errors, total_error_probability
                i+=1
    return dic_read_properties


def parse_reads(fast=True):
    print('parsing reads')
    if fast:
        if os.path.isfile(READ_PATH.format('s')): return
        i = 0
        read=""
        with open(f'{DATA_FOLDER}/reads.fastq') as f:
            for line in f:
                if line[0] == "@":
                    read += f'{next(f).strip()}$'
            i += 1
        writer = open(READ_PATH.format('s'), 'w')
        writer.write('>all the reads\n')
        writer.write(read)
        return
    i = 0
    with open('data/reads.fastq') as f:
        for line in f:
            if line[0] == "@":
                f2 = open(READ_PATH.format(i), 'w')
                f2.write(f">{line[1:]}")
                f2.write(next(f))
                i += 1
                f2.close()
    return i

def run_bdbwt(mode, number_of_reads, k_values, target, fast=True):
    print(f'running bdbwt-mem for {mode}')
    if fast:
        if not os.path.isfile(CONFIG_PATH[mode].format('')):
            generate_congif_bdbwt_mem(
                CONFIG_PATH[mode].format(''), target, READ_PATH.format('s'), min(k_values), 0 if mode == BDBWT_MEM else 1)
        if not os.path.isfile(ANCHOR_PATH[mode].format('s')):
            subprocess.run(ANCHOR_ALGOS[mode].format(
                CONFIG_PATH[mode].format(''), ANCHOR_PATH[mode].format('s')), shell=True)
        return
    
    for i in range(number_of_reads):
        print(f'read {i}')
        if not os.path.isfile(CONFIG_PATH[mode].format(i)):
            generate_congif_bdbwt_mem(
                CONFIG_PATH[mode].format(i), target, READ_PATH.format(i), k_values[i], 0 if mode == BDBWT_MEM else 1)
        if not os.path.isfile(ANCHOR_PATH[mode].format(i)):
            subprocess.run(ANCHOR_ALGOS[mode].format(
                CONFIG_PATH[mode].format(i), ANCHOR_PATH[mode].format(i)), shell=True)


def run_mummer(mode, number_of_reads, k_values, target_path, fast=True):
    print(f'running mummer for {mode}')
    if fast:
        if not os.path.isfile(ANCHOR_PATH[mode].format('s')):
            subprocess.run(ANCHOR_ALGOS[mode].format(
                min(k_values), READ_PATH.format('s'), target_path, ANCHOR_PATH[mode].format('s')), shell=True)
        return
    
    for i in range(number_of_reads):
        print(f'read {i}')
        if not os.path.isfile(ANCHOR_PATH[mode].format(i)):
            subprocess.run(ANCHOR_ALGOS[mode].format(
                k_values[i], READ_PATH.format(i), target_path, ANCHOR_PATH[mode].format(i)), shell=True)


def run_minimap(number_of_reads, k_values, target_path, fast=True):
    print(f'running minimap ')
    if fast:
        if not os.path.isfile(ANCHOR_PATH[MINIMAP].format('s')):
            subprocess.run(ANCHOR_ALGOS[MINIMAP].format(
                min(k_values), READ_PATH.format('s'), target_path, ANCHOR_PATH[MINIMAP].format('s')), shell=True)
            
        return
    
    for i in range(number_of_reads):
        print(f'read {i}')
        if not os.path.isfile(ANCHOR_PATH[MINIMAP].format(i)):
            subprocess.run(ANCHOR_ALGOS[MINIMAP].format(
                k_values[i], READ_PATH.format(i), target_path, ANCHOR_PATH[MINIMAP].format(i)), shell=True)


def run_chainx(number_of_reads, target_path):
    print('running chainX')
    for type_of, path in CHAIN_PATH.items():
        for i in range(number_of_reads):
            if not os.path.isfile(path.format(i)):
                subprocess.run(
                    f"./ChainX/chainX -m sg -q {READ_PATH.format(i)} -t {target_path} --anchors {TIDY_ANCHOR_PATH[type_of].format(i)} &>> {path.format(i)}", shell=True)


def parse_anchors(number_of_reads):
    print('parsing anchors')
    for type_of, path in ANCHOR_PATH.items():
        for i in range(number_of_reads):
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


def main(target="data/ecoli.fasta", constant_k=False, k=0, number_of_reads=10, coverage_of_reads=1.0):
    fresh_run()
    # generate anchors
    init_directories()
    if number_of_reads > 0:
        if not os.path.isfile(f"{DATA_FOLDER}/reads.fastq"):
            subprocess.run(
                f"simlord --read-reference {target} -n {number_of_reads} --no-sam data/reads", shell=True)
    else:
        if not os.path.isfile(f"{DATA_FOLDER}/reads.fastq"):
            subprocess.run(
                f"simlord --read-reference {target} -c {coverage_of_reads} --no-sam data/reads", shell=True)
    parse_reads()
    read_properties = get_read_properties(f'{DATA_FOLDER}/reads.fastq')
    number_of_reads = len(read_properties)
    print(number_of_reads)
    if constant_k:
        print('costant k')
        k_values = [k for _ in range(number_of_reads)]
        if k > 28:
            minimap_k_values = [28 for _ in range(number_of_reads)]
        else:
            minimap_k_values = [k for _ in range(number_of_reads)]
    else:
        k_values, minimap_k_values = init_variable_k(target, read_properties)
    
    run_bdbwt(BDBWT_EXT_MINI, number_of_reads, k_values, target)
    run_bdbwt(BDBWT_MEM, number_of_reads, k_values, target)
    run_mummer(MUMMER_MEM, number_of_reads, k_values, target)
    run_mummer(MUMMER_MUM, number_of_reads, k_values, target)
    run_minimap(number_of_reads, minimap_k_values, target)
    print('äend')
    # parse_anchors(number_of_reads)

    #run_chainx(number_of_reads, target)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--target", help="path to target file")
    parser.add_argument("-m", "--mode",
                        help="Constant k mode, if True use constant k, if false use variable k. Default: True", default=True)
    parser.add_argument("-k", "--k-size",
                        help="Size of constant k. Default: 22", default=22)
    parser.add_argument("-rn", "--number_of_reads",
                        help="number of reads", default=None)
    parser.add_argument("-rc", "--read_coverage",
                        help="Default function, coverage of the reads", default=1.0)

    args = parser.parse_args()

    target = args.target
    constant_k = args.mode
    k = args.k_size
    number_of_reads = args.number_of_reads
    read_coverage = args.read_coverage

    main()
