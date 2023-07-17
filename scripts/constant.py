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
ANCHOR_ALGOS = {BDBWT_MEM: "./{4}bdbwt-mem/main {0} >> {1}",
                BDBWT_EXT_MINI: "./{4}bdbwt-mem/main {0} >> {1}",
                MUMMER_MUM: "./{4}mummer/mummer -mum -l {0} {1} {2} >> {3}",
                MUMMER_MEM: "./{4}mummer/mummer -maxmatch -l {0} {1} {2} >> {3}",
                MINIMAP: "./{4}minimap2/minimap2 -k {0} {1} {2} --print-seeds > {3} 2>&1"}
CONFIG_DIR_MEM = f'bdbwt-mem/configs/{BDBWT_MEM}/'
CONFIG_DIR_EXT_MINI = f'bdbwt-mem/configs/{BDBWT_EXT_MINI}/'
DIRS = {'read': READS_DIR,
        f'config {BDBWT_MEM}': CONFIG_DIR_MEM,
        f'config {BDBWT_EXT_MINI}': CONFIG_DIR_EXT_MINI,
        'chain': {x: f'{CHAIN_DIR}{x}/' for x in ANCHOR_ALGOS.keys()},
        'anchor': {x: f'{ANCHOR_DIR}{x}/' for x in ANCHOR_ALGOS.keys()},
        'anchor-tidy': {x: f'{TIDY_ANCHOR_DIR}{x}/' for x in ANCHOR_ALGOS.keys()}}

READ_PATH = f'{READS_DIR}'+'read{0}.fasta'
CONFIG_PATH = {BDBWT_MEM: DIRS[f"config {BDBWT_MEM}"]+'config{0}',
               BDBWT_EXT_MINI: DIRS[f"config {BDBWT_EXT_MINI}"]+'config{0}'}
ANCHOR_PATH = {anchor_type: DIRS['anchor'][anchor_type] +
               'read{0}.txt' for anchor_type in ANCHOR_ALGOS.keys()}
TIDY_ANCHOR_PATH = {anchor_type: DIRS['anchor-tidy'][anchor_type] +
                    'read{0}.txt' for anchor_type in ANCHOR_ALGOS.keys()}
CHAIN_PATH = {anchor_type: DIRS['chain'][anchor_type] +
              'read{0}.txt' for anchor_type in ANCHOR_ALGOS.keys()}
