import re
import os

BDBWT_MEM = 'bdbwt-mem'
BDBWT_EXT_MINI = 'bdbwt-ext-mini'
BDBWT_MUM = 'bdbwt-mum'
MUMMER_MUM = 'mummer-mum'
MUMMER_MEM = 'mummer-mem'
MINIMAP = 'minimap'
EXTENDED_MINIMAP = 'ext-minimap'
ANCHOR_TYPES = [BDBWT_MEM, BDBWT_EXT_MINI, BDBWT_MUM,
                MUMMER_MEM, MUMMER_MUM, MINIMAP, EXTENDED_MINIMAP]
TARGET_STR = ''

DATA_FOLDER = "data/"
ANCHOR_DIR = f'{DATA_FOLDER}anchors/'
ANCHOR_STATS_DIR = f'{DATA_FOLDER}anchors-stats/'
TIDY_ANCHOR_DIR = f'{DATA_FOLDER}anchors-tidy/'
CHAIN_DIR = f'{DATA_FOLDER}chains/'
READS_DIR = f'{DATA_FOLDER}reads/'
RESULT_FOLDER = 'results/'
BENCHMARK_FOLDER = 'benchmarks/'
BENCHMARK_ANCHORS_DIR = f'{BENCHMARK_FOLDER}anchors/'
BENCHMARK_CHAINS_DIR = f'{BENCHMARK_FOLDER}chains/'
ANCHOR_ALGOS = {BDBWT_MEM: "./{2}bdbwt-mem/main {0} > {1}",
                BDBWT_EXT_MINI: "./{2}bdbwt-mem/main {0} > {1}",
                MUMMER_MUM: "./{4}mummer/mummer -mum -l {0} {1} {2} > {3}",
                MUMMER_MEM: "./{4}mummer/mummer -maxmatch -l {0} {1} {2} > {3}",
                MINIMAP: "./{4}minimap2/minimap2 -k {0} {1} {2} --print-seeds > {3} 2>&1"}
CONFIG_DIR_MEM = f'bdbwt-mem/configs/{BDBWT_MEM}/'
CONFIG_DIR_EXT_MINI = f'bdbwt-mem/configs/{BDBWT_EXT_MINI}/'
DIRS = {'read': READS_DIR,
        f'config {BDBWT_MEM}': CONFIG_DIR_MEM,
        f'config {BDBWT_EXT_MINI}': CONFIG_DIR_EXT_MINI,
        'chain': {x: f'{CHAIN_DIR}{x}/' for x in ANCHOR_TYPES},
        'anchor': {x: f'{ANCHOR_DIR}{x}/' for x in ANCHOR_TYPES},
        'anchor-stats': ANCHOR_STATS_DIR,
        'anchor-tidy': {x: f'{TIDY_ANCHOR_DIR}{x}/' for x in ANCHOR_TYPES},
        'results': RESULT_FOLDER,
        'benchmarks-anchors': {x: f'{BENCHMARK_ANCHORS_DIR}{x}/' for x in ANCHOR_TYPES},
        'benchmarks-chains': {x: f'{BENCHMARK_CHAINS_DIR}{x}/' for x in ANCHOR_TYPES}}

READ_PATH = f'{READS_DIR}'+'read{0}.fasta'
CONFIG_PATH = {BDBWT_MEM: DIRS[f"config {BDBWT_MEM}"]+'config{0}',
               BDBWT_EXT_MINI: DIRS[f"config {BDBWT_EXT_MINI}"]+'config{0}'}
ANCHOR_PATH = {anchor_type: DIRS['anchor'][anchor_type] +
               'read{0}.txt' for anchor_type in ANCHOR_TYPES}
ANCHOR_STATS_PATH = f'{ANCHOR_STATS_DIR}stat.csv'
TIDY_ANCHOR_PATH = {anchor_type: DIRS['anchor-tidy'][anchor_type] +
                    'read{0}.txt' for anchor_type in ANCHOR_TYPES}
CHAIN_PATH = {anchor_type: DIRS['chain'][anchor_type] +
              'read{0}.txt' for anchor_type in ANCHOR_TYPES}
CHAIN_SUMMARY_PATH = f'{RESULT_FOLDER}chain-summary'+'-{k}-{genome}.csv'
ANCHOR_SUMMARY_PATH = f'{RESULT_FOLDER}anchor-summary'+'-{k}-{genome}.csv'
BENCHMARK_ANCHOR_PATH = {anchor_type: DIRS['benchmarks-anchors']
                         [anchor_type] + 'read{0}.txt' for anchor_type in ANCHOR_TYPES}
BENCHMARK_CHAIN_PATH = {anchor_type: DIRS['benchmarks-chains']
                        [anchor_type] + 'read{0}.txt' for anchor_type in ANCHOR_TYPES}


def get_target(target_path: str) -> str:
    target = ''
    with open(target_path) as f:
        for j, line in enumerate(f):
            if j != 0:
                target += line.strip()
    return target


def get_read(id):
    read_properties = get_read_properties(READ_PATH.format(id))
    return (read_properties[1], 0, read_properties[0])


def get_reads(input_folder):
    dictionary_of_tuple_lists = {}
    for _, _, files in os.walk(input_folder):
        for fi in files:
            try:
                read_number = re.findall(r"\d+", fi)[0]
            except:
                continue
            read_properties = get_read_properties(f'{input_folder}{fi}')
            dictionary_of_tuple_lists[int(read_number)] = [(
                read_properties[1], 0, read_properties[0])]
    return dictionary_of_tuple_lists


def get_tuple_list_from_file(input_folder, include_empty=False):
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
                        tuple_list.append((x, y, length))
                    except:
                        continue
            if len(tuple_list) > 0 and not include_empty:
                dictionary_of_tuple_lists[int(read_number)] = tuple_list
            elif include_empty:
                dictionary_of_tuple_lists[int(read_number)] = tuple_list
    return dictionary_of_tuple_lists


def get_read_properties(read_file_path):
    with open(read_file_path) as f:
        read_properties = f.readline().split(';')
        read_length = int(re.findall("\d+", read_properties[1])[0])
        read_start_position = int(re.findall("\d+", read_properties[2])[0])
        read_chromosome = read_properties[3]
        number_of_errors = int(re.findall("\d+", read_properties[4])[0])
        total_error_probability = float(
            re.findall("\d+.\d+", read_properties[5])[0])
        return read_length, read_start_position, read_chromosome, number_of_errors, total_error_probability
