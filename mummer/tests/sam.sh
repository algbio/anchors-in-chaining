nucmer --maxmatch --sam-long /dev/stdout $D/seed_reads_1.fa $D/seed_reads_0.fa | check_cigar /dev/stdin $D/seed_reads_1.fa $D/seed_reads_0.fa
nucmer --maxmatch --sam-short /dev/stdout $D/seed_fq_genome.fa $D/seed_fq_reads_0.fq | tail -n +3 | sort | tee tee.sam | test_md5 4a3a273f3c543e2f9eaf2edf95a60821
nucmer --sam-long /dev/stdout -l 10 <(echo -e ">101\nggtttatgcgctgttatgtctatggacaaaaaggctacgagaaactgtagccccgttcgctcggacccgcgtcattcgtcggcccagctctacccg") <(echo -e ">21\nggtttatgcgctgttttgtctatggaaaaaaggctacgagaaactgtagccccgttcgctcggtacccgcgtcattcgtcggcccatctctacccg") | tail -n +3 | test_md5 f656b26b59de04e7c94c7c0c0f7e3a0c