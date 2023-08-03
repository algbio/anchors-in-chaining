import numpy as np

from generate_reads import *
from run import *
from constant import *
#sort of sanity check
#generate reads without errors and see that results are correct
def main():
    target = get_target('data/ecoli.fasta')
    read_length = 8000
    rng = np.random.default_rng()
    start_pos = rng.integers(0,len(target)-read_length, 25)
    for i, pos in enumerate(start_pos):
        read = target[pos:pos+read_length]
        writer = open('data/reads.fastq', 'a')
        writer.write(f'@m0/13111/CCS Read={i};length={read_length}bp;startpos={pos};chromosome=U00096.3;numberOfErrors=0;totalErrorProb=0.0;passes=0;passesLeft=0;passesRight=0;cutPosition=0\n')
        writer.write(f'{read}\n')
        writer.write('x\n')
        writer.write('XXX\n')

if __name__ == '__main__':
    main()