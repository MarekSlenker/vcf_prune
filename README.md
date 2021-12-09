# VCF_prune

This code was created in collaboration with [Jordan Koch](https://github.com/jmkoch) and [Patrick Monnahan](https://github.com/pmonnahan)
 
## Description
`VCF_prune.py` is a python2 script that takes **multi-sample (any ploidy or a mix of ploidies) vcf files as input**, and outputs files that are formatted for use in **STRUCTURE**.  
This script works by walking along the chromosome in each input vcf file. It then randomly selects a single SNP in each window (passed filtering; width of window and distances between windows are defined by the user in the `–w` and `-d`  arguments) and translates the nucleotide bases to integers. It aligns these integers in a matrix and labels and arranges them according to the necessary STRUCTURE input data format.


## Usage
**Dependencies**: numpy, scipy, and transposer  
To see the required arguments for this script, type `python2 VCF_prune.py` into the command line.  

```
VCF_prune.py [<args>]

The VCF_prune.py args are:
  -v        STR     path to vcfs (required)
  -w        INT     size of scaffold window (required)
  -d        INT     distance between any 2 windows (required)

  -m        FLOAT   amount of missing data to allow per site (between 0-1; required)
  -minf     FLOAT   minimum allele frequency, ALT or REF (between 0-1) [0]
  -maxf     FLOAT   maximum allele frequency, ALT or REF (between 0-1) [1]
  -minSnps  INT     minimal amount of SNPs in window. If less, window will be skipped. [1]
  
  -o        STR     output prefix (required)
  -r        INT     number of replicate data sets [1]
  
  -p        INT     length of population name (required)
  -n        INT     ploidy of output Structure file (required)

  -s                subsample polyploid data to create psuedo-diploid data
  -gz               use if vcfs are gzipped
  -vcf              if used, pruned VCF files will be printed
```
This script reads each SNP and evaluates if it passed the desired threshold (set by arguments `-m`, `-minf`, and `-maxf`). If it does, the SNP is recorded. It records each SNP passing threshold within a single window, whose width is set by `-w` argument. When passing the window border, the script randomly selects a single SNP for the STRUCTURE file. Searching for SNPs continues in the next window. Distance in between two windows is set by `-d` argument.  

POPDATA: `VCF_prune.py` writes POPDATA column to STRUCTURE file. Sample's membership in the population is inferred according to the sample names in the header line of vcf files. Let's consider the following header,  
`#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT  AL004_25  AL004_28  AL004_29  AL004_30  AL005_11  AL005_13`  
where underscore character delimits population name and individual name. In this case, the length of the population name is 5 (`-p 5`).

PLOIDY: `VCF_prune.py` works with mixed ploidy vcf files. Set highest ploidy by `-n` argument, and the absence of an allele (e.g.,  if ploidy is 4, diploids has 2 alleles absenting ) will be coded as -9 in the resulted STRUCTURE file. 

`-r` stands for the number of randomly replicated pruned datasets you want to create.  
`-s` creates pseudo-diploid (‘diploidized’) data. This is required for fastStructure, as fastStructure can only handle diploid data. So, if you are running tetraploid or mixed_ploidy datasets, you will need to specify that you want diploidized data. 


## Examples
To run `VCF_prune.py` with test_data, execute the following command:
```
python2 VCF_prune.py -v test_data -w 100 -d 0 -m 0.4 -minf 0.05 -maxf 0.95 -p 5 -n 4 -o test -s
```
results will appear in test_data directory. The resulted structure files are ready to run.  

&nbsp;  
### How to cite
**Use the following formula:** ... vcf files were pruned using the vcf_prune.py script available at https://github.com/MarekSlenker/vcf_prune.  

&nbsp;  
### Other questions not covered here and reporting problems
If you have a question or you encounter a problem, feel free to email me at marek.slenker@savba.sk, or use [issues](https://github.com/MarekSlenker/vcf_prune/issues). I will do my best to help you.
