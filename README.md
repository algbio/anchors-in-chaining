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

Snaemake file is run with command 
```
snakemake --cores 1
```