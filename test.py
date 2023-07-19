import subprocess


#read0
TARGET_STR = get_target('data/ecoli.fasta')
start_pos = 350317
length = 13111 
sub = TARGET_STR[start_pos:start_pos+length]

