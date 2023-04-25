Comparing Minimizers and MEMS as anchors in chaining

Chaining algorithm: 
* modified ChainX, [algbio-ChainX](https://github.com/nrizzo/algbio-ChainX)
* original ChainX, [ChainX](https://github.com/at-cg/ChainX)
```
./chainX -m g -q data/time_global/mutated_80_perc.fasta -t data/time_global/Chromosome_2890043_3890042_0.fasta
```

Anchors:
* asymmetric MEM
* symmetric MEM, r-index, [br-index-mem](https://github.com/algbio/br-index-mems)
* MUM
* Minimizers

Existing implementation of chaining + different anchors [bdbwt-mem-chaining](https://github.com/algbio/bdbwt-mem-chaining)
