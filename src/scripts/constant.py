import re
import os

BDBWT_MEM = 'bdbwt-mem'
BDBWT_EXT_MINI = 'bdbwt-ext-mini'
BDBWT_MUM = 'bdbwt-mum'
MUMMER_MUM = 'mummer-mum'
MUMMER_MEM = 'mummer-mem'
MINIMAP = 'minimap'
EXTENDED_MINIMAP = 'ext-minimap'
BR_INDEX_MEM = 'br-index-mem'
BR_INDEX_MUM = 'br-index-mum'
ANCHOR_TYPES = [BDBWT_MEM, BDBWT_EXT_MINI, BDBWT_MUM,
                MUMMER_MEM, MUMMER_MUM, MINIMAP, EXTENDED_MINIMAP,
                BR_INDEX_MEM, BR_INDEX_MUM]

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
R_INDEX_ALGO = './{1}br-index-mems/build/bri-build {0}'
ANCHOR_ALGOS = {BDBWT_MEM: "./{2}bdbwt-mem/main {0} > {1}",
                BDBWT_EXT_MINI: "./{2}bdbwt-mem/main {0} > {1}",
                MUMMER_MUM: "./{4}mummer/mummer -mum -l {0} {1} {2} > {3}",
                MUMMER_MEM: "./{4}mummer/mummer -maxmatch -l {0} {1} {2} > {3}",
                MINIMAP: "./{4}minimap2/minimap2 -k {0} {1} {2} --print-seeds > {3} 2>&1",
                BR_INDEX_MEM: "./{4}br-index-mems/build/bri-mem -k {0} -o {1} {2} {3}"}
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
    # return target string from target file
    target = ''
    with open(target_path) as f:
        for j, line in enumerate(f):
            if j != 0:
                target += line.strip()
    return target


def get_read(id) -> tuple:
    # get tuple (position at the target, position at the read, read length)
    read_properties = get_read_properties(READ_PATH.format(id))
    return (read_properties[1], 0, read_properties[0])


def get_read_str(read_file_path) -> str:
    # get read string
    with open(read_file_path) as f:
        for line in f:
            if line[0] != '>':
                return line.strip()


def get_reads(input_folder) -> dict:
    # get all the reads from the read folder
    # dictionary structure {read_id: [(position at the target, position at the read, read length)]..}
    dictionary_of_tuple_lists = {}
    for _, _, files in os.walk(input_folder):
        for fi in files:
            if 'bri' not in fi:
                try:
                    read_number = re.findall(r"\d+", fi)[0]
                except:
                    continue
                read_properties = get_read_properties(f'{input_folder}{fi}')
                dictionary_of_tuple_lists[int(read_number)] = [(
                    read_properties[1], 0, read_properties[0])]
    return dictionary_of_tuple_lists


def get_tuple_list_from_file(input_folder, include_empty=False) -> dict:
    # reads tidy anchor file or chain file
    # returns {read_id: [(position at the target, position at the read, anchor length)..]..}
    dictionary_of_tuple_lists = {}
    for _, _, files in os.walk(input_folder):
        for fi in files:
            read_number = re.findall(r"\d+", fi)[0]
            tuple_list = get_tuple_list(f'{input_folder}{fi}')
            if len(tuple_list) > 0 and not include_empty:
                dictionary_of_tuple_lists[int(read_number)] = tuple_list
            elif include_empty:
                dictionary_of_tuple_lists[int(read_number)] = tuple_list
    return dictionary_of_tuple_lists


def get_tuple_list(file_path) -> list:
    tuple_list = []
    with open(file_path) as f:
        for line in f:
            try:
                parts = re.findall(r"\d+", line)
                x = int(parts[0])
                y = int(parts[1])
                length = int(parts[2])
                tuple_list.append((x, y, length))
            except:
                continue
    return tuple_list


def get_read_properties(read_file_path) -> tuple:
    # returns read_length, read_start_position, read_chromosome, number_of_errors, total_error_probability
    with open(read_file_path) as f:
        read_properties = f.readline().split(';')
        read_length = int(re.findall("\d+", read_properties[1])[0])
        read_start_position = int(re.findall("\d+", read_properties[2])[0])
        read_chromosome = read_properties[3]
        number_of_errors = int(re.findall("\d+", read_properties[4])[0])
        total_error_probability = float(
            re.findall("\d+.\d+", read_properties[5])[0])
        return read_length, read_start_position, read_chromosome, number_of_errors, total_error_probability
