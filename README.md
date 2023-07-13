Comparing Minimizers and MEMS as anchors in chaining

Chaining algorithm: 
* modified ChainX, [algbio-ChainX](https://github.com/nrizzo/algbio-ChainX)
* original ChainX, [ChainX](https://github.com/at-cg/ChainX)

Anchors:
* asymmetric MEM
* symmetric MEM, r-index, [br-index-mem](https://github.com/algbio/br-index-mems)
* MUM
* Minimizers

Existing implementation of chaining + different anchors [bdbwt-mem-chaining](https://github.com/algbio/bdbwt-mem-chaining)

Snakemake file is run with command 
```
snakemake --cores 1
```

TODO:
* ~~run on e coli sequence~~
* ~~generate reads~~
    > `simlord  --read-reference data/test/ecoli.fasta -n 10000 --no-sam  data/test/ecoli_reads`
    > https://academic.oup.com/bioinformatics/article/32/17/2704/2450740
* confirm that anchors are matches
    > mostly matches
    > **what is the fastest method to compute the anchors?**
        > compare: 
        1. MUMMER, mem and mum, 
        2. bdbwt, extended minimizer, MEMs 
        3. r-index, MEM, MUM?
        4. minimap2 minimizers, extended minimizers with minimap2 minimizers
    > bdbwt is horrible computing mems when compared to chainX:s mummer, with e-coli.
* compute MUMS
    > use the ChainX mummer implementation or implement ur own
* ~~k mer size dependent on the length of the sequence as in the paper. 22 for e coli 31 for huma~~
* ~~report number of anchors and the sum of anchor bases.~~
* ~~compute the Jaccard index and draw siome figures~~

Stats to report:
* average length of reads
* average number of anchor per chain
* average length of anchor bases
* average number of anchors per read
* anchor coverage of read
* chain coverage of read
* Jaccard index (the shared bases between the chain and the read)/union of the chain and the read