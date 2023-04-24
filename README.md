Comparing Minimizers and MEMS as anchors in chaining

Chaining algorithm: 
* modified ChainX, [algbio-ChainX](https://github.com/nrizzo/algbio-ChainX)
* original ChainX, [ChainX](https://github.com/at-cg/ChainX)

Anchors:
* asymmetric MEM
* symmetric MEM, r-index, [br-index-mem](https://github.com/algbio/br-index-mems)
```    
mkdir inputs
mkdir outputs
cp ../patterns.txt inputs/
cp ../text.txt inputs/
./br-build inputs/patterns.txt
./br-build inputs/text.txt
./bri-mem -k 4 -o outputs/MEMs.txt inputs/text.txt.bri inputs/patterns.txt.bri
cat outputs/MEMs.txt
```
* MUM
* Minimizers

Existing implementation of chaining + different anchors [bdbwt-mem-chaining](https://github.com/algbio/bdbwt-mem-chaining)
