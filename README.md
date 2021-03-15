# VCF_prune

This code was created in collaboration with [Jordan Koch](https://github.com/jmkoch)
 
## Description
`VCF_prune.py` is a python2 script that takes **multi-sample (any ploidy or mix of ploidies) vcf files as input**, and outputs files that are formatted for use in **STRUCTURE**.  
This script works by walking along the chromosome in each input vcf file. It then randomly selects a single SNP in each window (passed filtering; width of window and distances between windows are defined by the user in the `–w` and `-d`  arguments), and translates the nucleotide bases to integers. It aligns these integers in a matrix, and labels and arranges them according to the necessary STRUCTURE input data format.


## Usage
**Dependencies**: numpy, scipy, and transposer  
To see the required arguments for this script, type `python2.7 VCF_prune.py` into the command line.  

```
VCF_prune.py [<args>]

The VCF_prune.py args are:
  -v     STR     path to vcfs (required)
  -w     INT     size of scaffold window (required)
  -d     INT     distance between any 2 windows (required)

  -m     FLOAT   amount of missing data to allow per site (between 0-1; required)
  -minf  FLOAT   minimum allele frequency, ALT or REF (between 0-1) [0]
  -maxf  FLOAT   maximum allele frequency, ALT or REF (between 0-1) [1]

  -o     STR     output prefix (required)
  -r     INT     number of replicate data sets [1]
  
  -p     INT     length of population name (required)
  -n     INT     ploidy of output Structure file (required)

  -s             subsample polyploid data to create psuedo-diploid data
  -gz            use if vcfs are gzipped
```
This script reads each SNP and evaluates if it passed the desired threshold (seted by arguments `-m`, `-minf`, and `-maxf`). If it does, the SNP is recorded. It records each SNP passing threshold within a single window, whose width is set by `-w` argument. When passing the window border, the script randomly selects a single SNP for the STRUCTURE file. Searching for SNPs continues in the next window. Distance in between two windows is set by `-d` argument.


POPDATA

-n
-r
-s
-gz

