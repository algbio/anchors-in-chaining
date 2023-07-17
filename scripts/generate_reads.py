import os
import argparse
import shutil
import subprocess
from constant import *

def fresh_run():
    for type_of, dir in DIRS.items():
        if type_of == 'chain' or type_of == 'anchor' or type_of == 'anchor-tidy':
            for _, sub_dir in dir.items():
                shutil.rmtree(sub_dir, ignore_errors=True)
        else:
            shutil.rmtree(dir, ignore_errors=True)


def init_directories():
    for type_of, dir in DIRS.items():
        if type_of == 'chain' or type_of == 'anchor' or type_of == 'anchor-tidy':
            for _, sub_dir in dir.items():
                if not os.path.exists(sub_dir):
                    os.makedirs(sub_dir)
        else:
            if not os.path.exists(dir):
                os.makedirs(dir)

def parse_reads():
    if os.path.isfile(READ_PATH.format('s')):
        return
    i = 0
    read = ""
    with open(f'{DATA_FOLDER}/reads.fastq') as f:
        for line in f:
            if line[0] == "@":
                read += f'{next(f).strip()}$'
        i += 1
    writer = open(READ_PATH.format('s'), 'w')
    writer.write('>all the reads\n')
    writer.write(read)

    i = 0
    with open('data/reads.fastq') as f:
        for line in f:
            if line[0] == "@":
                f2 = open(READ_PATH.format(i), 'w')
                f2.write(f">{line[1:]}")
                f2.write(next(f))
                i += 1
                f2.close()
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fresh_run",
                        help="Boolean, Clear directories", action='store_true')
    parser.add_argument("-t", "--target", help="path to target file")
    parser.add_argument("-rc", "--coverage",
                        help="coverage of reads", default=-1, type=float)
    parser.add_argument("-rn", "--number", help="number of reads", default=-1, type = int)

    args = parser.parse_args()

    boolean_fresh_run = args.fresh_run
    target = args.target
    coverage_of_reads = args.coverage
    number_of_reads = args.number

    if boolean_fresh_run:
        fresh_run()
    init_directories()
    if number_of_reads > 0:
        subprocess.run(
            f"simlord --read-reference {target} -n {number_of_reads} --no-sam data/reads", shell=True)
    if coverage_of_reads > 0:
        subprocess.run(
            f"simlord --read-reference {target} -c {coverage_of_reads} --no-sam data/reads", shell=True)
    # parse reads
    parse_reads()