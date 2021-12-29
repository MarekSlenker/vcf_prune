# VCF_prune v 1.1.0

The initial version of this script was created by [Jordan Koch](https://github.com/jmkoch) and [Patrick Monnahan](https://github.com/pmonnahan). 
 
## Description
`VCF_prune.py` is a python2 script that takes **multi-sample (any ploidy or a mix of ploidies) vcf files as input**, and outputs files that are formatted for use in **STRUCTURE**.  

This script works in two modes: 
* by walking along the scafold of each input vcf file. A random SNP is selected from each window (passed filtering; width of window and distances between windows are defined by the user by the `–-winSize` and `--winDist`  arguments). Ideal for data retrieved from **whole-genome resequencing**,  where the SNPs are +- continuously spread over the scaffolds.
* select random SNP from each of provided loci (defined in '--regions' argument). This option is ideal for the SNPs produced by **RadSeq**, as those SNPs are not spread continuously, but grouped in RAD loci.  

SNPs are translated to integers and outputed in the STRUCTURE input data format.



## Usage
**Dependencies**: numpy, and transposer  
To see the required arguments for this script, type `python2 VCF_prune.py` into the command line.  

```
VCF_prune.py [<args>]

The VCF_prune.py args are:
  --inVcfs    STR     path to vcfs (required)
  --winSize   INT     size of scaffold window (required, if --regions is not specified)
  --winDist   INT     distance between any 2 windows (required, if --regions is not specified)
  --regions   STR     path to file defining regions from which to take SNPs (1 SNP per region) format: "CHR:FROM-TO" (required, if --winSize and --winDist is not specified)

  --missing   FLOAT   amount of missing data to allow per site (between 0-1) [1]
  --minf      FLOAT   minimum allele frequency, ALT or REF (between 0-1) [0]
  --maxf      FLOAT   maximum allele frequency, ALT or REF (between 0-1) [1]
  --minSnps   INT     Minimal amount of SNPs in window. If less, window will be skipped [1]

  --prefix    STR     output prefix
  --reps      INT     number of replicates to produce (data sets) [1]
  
  --popFlagLength  INT     length of population name
  --ploidy         INT     ploidy of output Structure file [2]

  --subsample        subsample polyploid data to create psuedo-diploid data
  --gz               use if vcfs are gzipped
  --vcf              use if you want to print also pruned VCF files 
```



**--popFlagLength**: If this argument is used, `VCF_prune.py` writes POPDATA column to STRUCTURE file. Sample's membership in the population is inferred according to the sample names in the header line of vcf files. Let's consider the following header,  
`#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT  AL004_25  AL004_28  AL004_29  AL004_30  AL005_11  AL005_13`  
where underscore character delimits population name and individual name. In this case, the length of the population name is 5 (`--popFlagLength 5`).

**--ploidy**: `VCF_prune.py` works with mixed ploidy vcf files. Set highest ploidy by `--ploidy` argument, and the absence of an allele (e.g.,  if ploidy is 4, diploids has 2 alleles absenting) will be coded as `-9` in the resulted STRUCTURE file.  

**--regions**: Regions are provided in text file, one region per row, format: "CHR:FROM-TO". To find region in RadSeq data, use `identifiRadLoci.workflow`. See comments inside! You will need bam files for this step. Plot the distances and set the `betweenLociDist` threshold (if the distance between two potential loci is lower (or equal) than betweenLociDist, loci are considered as a single locus. If the distance is greater than betweenLociDist, loci are considered to be two different loci). See Rplots.pdf for example of histogram. 


**--reps** stands for the number of randomly replicated pruned datasets you want to create.  

**--subsample** creates pseudo-diploid (‘diploidized’) data. This is required for fastStructure, as fastStructure can only handle diploid data. So, if you are running tetraploid or mixed_ploidy datasets, you will need to specify that you want diploidized data. 


## Examples
To run `VCF_prune.py` with test_data, execute the following command:
```
python2 VCF_prune.py --inVcfs test_data --winSize 10000 --winDist 1000 --missing 0.1 --minf 0.05 --maxf 0.95 --popFlagLength 5 --ploidy 4 --prefix test --gz

python2 VCF_prune.py --inVcfs test_data --regions test_data/regions1000apart.txt --ploidy 4 --minSnps 10 --reps 20 --gz

```
results will appear in test_data directory. The resulted structure files are ready to run.  




&nbsp;  
### How to cite
**Use the following formula:** ... vcf files were pruned using the vcf_prune.py script available at https://github.com/MarekSlenker/vcf_prune.  

&nbsp;  
### Other questions not covered here and reporting problems
If you have a question or you encounter a problem, feel free to email me at marek.slenker@savba.sk, or use [issues](https://github.com/MarekSlenker/vcf_prune/issues). I will do my best to help you.
