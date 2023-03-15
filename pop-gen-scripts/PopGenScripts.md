## Pop Gen Scripts Markdown
Katherine Carbeck


> ***Table of Contents***
> 1. [PCA](#1-pca-on-final-319-individual-vcf-file-no-ld-pruning-yet)
> 2. [LD pruning](#2-ld-pruning)



### 1. PCA on final 319 individual VCF file

On terminal:

```bash
cp *319Samples* /workdir/kcarbeck/pca
cd /workdir/kcarbeck/pca
module load R/4.2.1-r9
R
```

Enter R environment on terminal:
```R
setwd("/workdir/kcarbeck/pca")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install(c("SNPRelate"))
install.packages("scales")
library(SNPRelate)
library(gdsfmt)


###new PCA with vcf done with --very-sensitive-local aligned data.
#reformats the vcf file to gds file for use in futher analysis
snpgdsVCF2GDS(vcf.fn="/workdir/kcarbeck/pca/filtered_SOSP_319Samples_final.vcf.gz", out.fn="319Samples_out.gds",
              method = c("biallelic.only"),compress.annotation="ZIP.max", 
              snpfirstdim=FALSE, verbose=TRUE)
     # Start file conversion from VCF to SNP GDS ...
     # Method: extracting biallelic SNPs
     # Number of samples: 319
     # Parsing "/workdir/kcarbeck/pca/filtered_SOSP_319Samples_final.vcf.gz" ...
     #         import 12504847 variants.
     # + genotype   { Bit2 319x12504847, 951.1M } *
     # Optimize the access efficiency ...
     # Clean up the fragments of GDS file:
     #     open the file '319Samples_out.gds' (1.0G)
     #     # of fragments: 1970
     #     save to '319Samples_out.gds.tmp'
     #     rename '319Samples_out.gds.tmp' (1.0G, reduced: 22.9K)
     #     # of fragments: 20

```

Then copy files out of workdir to local computer to run below code in RStudio:

```bash
cp 319Samples_out.gds /lustre2/home/lc736_0001/song_sparrow/final_vcf/pca
```

[Plot PCA in R](01a.plotPCA.R)






### 2. LD pruning 
[LD Pruning Script outlined below](02.LDpruning.sh)

One of the major assumptions of PCA/ADMIXTURE is that the data we use is indpendent - i.e. there are no spurious correlations among the measured variables. This is obviously not the case for most genomic data as allele frequencies are correlated due to physical linkage and linkage disequilibrium. So as a first step, we need to prune our dataset of variants that are in linkage.


#### 1. generate the input file in plink format:

```bash
cd pca
VCF=/workdir/kcarbeck/pca/filtered_SOSP_319Samples_final.vcf.gz
```


```bash
/programs/plink2_linux_x86_64_20221024/plink2 --vcf $VCF --allow-extra-chr \
 --double-id --rm-dup --new-id-max-allele-len 51 missing \
 --set-missing-var-ids @:#\$r:\$a --indep-pairwise 50 10 0.1 \
 --out SOSP_319Samples_LDpruning
```
*commands:*
     **--vcf** - specified the location of our VCF file.
     **--double-id** - told plink to duplicate the id of our samples (this is because plink typically expects a family and individual id - i.e. for pedigree data - this is not necessary for us.
    **--allow-extra-chr** - allow additional chromosomes beyond the human chromosome set. This is necessary as otherwise plink expects chromosomes 1-22 and the human X chromosome.
     **--new-id-max-allele-len** - whole-genome sequencing frequently contain variants which have not been assigned standard IDs. If you don't want to throw out all of that data, you'll usually want to assign them chromosome-and-position-based IDs.
    **--indep-pairwise** - finally we are actually on the command that performs our linkage pruning! The first argument, 50 denotes we have set a window of 50 Kb. The second argument, 10 is our window step size - meaning we move 10 bp each time we calculate linkage. Finally, we set an r2 threshold - i.e. the threshold of linkage we are willing to tolerate. Here we prune any variables that show an r2 of greater than 0.1.
     **--out** Produce the prefix for the output data.


#### 2. prune and create pca
   
``` bash
/programs/plink2_linux_x86_64_20221024/plink2 --vcf $VCF --double-id --allow-extra-chr --extract SOSP_319Samples_LDpruning.prune.in --make-bed --pca --out ldPruned_out
```


commands:
     **--extract** - this just lets plink know we want to extract only these positions from our VCF - in other words, the analysis will only be conducted on these.
     **--make-bed** - this is necessary to write out some additional files for another type of population structure analysis - a model based approach with admixture.
     **--pca** - fairly self explanatory, this tells plink to calculate a principal components analysis.


