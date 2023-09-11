## Pop Gen Scripts Markdown


### 1. LD pruning 
One of the major assumptions of PCA/ADMIXTURE is that the data we use is indpendent - i.e. there are no spurious correlations among the measured variables. This is obviously not the case for most genomic data as allele frequencies are correlated due to physical linkage and linkage disequilibrium. So as a first step, we need to prune our dataset of variants that are in linkage.

Required scripts: 
* [01.admixture.sh](01.admixture.sh)
* [ldpruning.sh](ldpruning.sh)

```bash
# no missing data
bcftools view -e 'GT[*]="mis"' SOSP_316Samples_final.vcf.gz > SOSP_316Samples_final.NoMiss.vcf.gz
# check num of snnps
bcftools stats SOSP_316Samples_final.NoMiss.vcf.gz > snps.stats &
    # 36,445
#bgzip SOSP_316Samples_final.NoMiss.vcf.gz &
#tabix SOSP_316Samples_final.NoMiss.vcf.gz.gz 

# Run ldpruning script from Joana Meier
bash ldpruning.sh SOSP_316Samples_final.NoMiss.vcf.gz & 
  # After filtering, kept 316 out of 316 Individuals
  # Outputting VCF file...
  # After filtering, kept 7412265 out of a possible 13561802 Sites
  # Run Time = 4927.00 seconds
bcftools stats SOSP_316Samples_final.NoMiss.LDpruned.vcf.gz > snps.stats 
    #20,860

# Generate the input file in plink format
vcftools --gzvcf SOSP_316Samples_final.NoMiss.LDpruned.vcf.gz --plink --out SOSP_NoMiss_LDpruned_plink 
    # After filtering, kept 21844 out of a possible 21844 Sites
    # Run Time = 3.00 seconds

/programs/plink-1.9-x86_64-beta3.30/plink --file SOSP_NoMiss_LDpruned_plink --make-bed --allow-extra-chr --out SOSP_NoMiss_LDpruned_plink &
```




### 2. ADMIXRURE

Required scripts: 
* [01.admixture.sh](01.admixture.sh)
* [plotAdmixture.R](plotAdmixture.R)

Use output from LD pruning step: SOSP_NoMiss_LDpruned_plink files



### 3. PCA on final 316 individual VCF file

Required scripts:
* [02.plotPCA.R](02.plotPCA.R)

On terminal:

```bash
cp *316Samples* /workdir/kcarbeck
cd /workdir/kcarbeck
module load R/4.2.1-r9
R
```

Enter R environment on terminal:
```R
setwd("/workdir/kcarbeck")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install()
#BiocManager::install(c("SNPRelate"))
#install.packages("scales")
library(SNPRelate)
#library(gdsfmt)


###new PCA with vcf done with --very-sensitive-local aligned data.
#reformats the vcf file to gds file for use in futher analysis
snpgdsVCF2GDS(vcf.fn="/workdir/kcarbeckSOSP_316Samples_final.vcf.gz", out.fn="316Samples_032123_out.gds",
              method = c("biallelic.only"),compress.annotation="ZIP.max", 
              snpfirstdim=FALSE, verbose=TRUE)
  # Start file conversion from VCF to SNP GDS ...
  # Method: extracting biallelic SNPs
  # Number of samples: 316
  # Parsing "/workdir/kcarbeck/pca/SOSP_316Samples_final.vcf.gz" ...
  #         import 12483733 variants.
  # + genotype   { Bit2 316x12483733, 940.5M } *
  # Optimize the access efficiency ...
  # Clean up the fragments of GDS file:
  #     open the file '316Samples_032123_out.gds' (1.0G)
  #     # of fragments: 1963
  #     save to '316Samples_032123_out.gds.tmp'
  #     rename '316Samples_032123_out.gds.tmp' (1.0G, reduced: 22.8K)
  #    # of fragments: 20

```

Then copy files out of workdir to local computer to run below code in RStudio:

```bash
cp 316Samples_032123_out.gds /lustre2/home/lc736_0001/song_sparrow/final_vcf/pca
```

Now, bring to local computer and plot in R.



### 4. PCA on LD pruned VCF

Required scripts:
* [02.plotPCA.R](02.plotPCA.R)

On terminal:

```bash
module load R/4.2.1-r9
R
```

Enter R environment on terminal:
```R
setwd("/workdir/kcarbeck/pca")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install()
#BiocManager::install(c("SNPRelate"))
#install.packages("scales")
library(SNPRelate)
#library(gdsfmt)


###new PCA with vcf done with --very-sensitive-local aligned data.
#reformats the vcf file to gds file for use in futher analysis
snpgdsVCF2GDS(vcf.fn="/workdir/kcarbeck/pca/SOSP_316Samples_final.LDpruned.vcf.gz", out.fn="316Samples_LDpruned_out.gds",
              method = c("biallelic.only"),compress.annotation="ZIP.max", 
              snpfirstdim=FALSE, verbose=TRUE)
  # Start file conversion from VCF to SNP GDS ...
  # Method: extracting biallelic SNPs
  # Number of samples: 316
  # Parsing "/workdir/kcarbeck/pca/SOSP_316Samples_final.LDpruned.vcf.gz" ...
  #         import 6841816 variants.
  # + genotype   { Bit2 316x6841816, 515.5M } *
  # Optimize the access efficiency ...
  # Clean up the fragments of GDS file:
  #   open the file '316Samples_LDpruned_out.gds' (565.2M)
  #   # of fragments: 1080
  #   save to '316Samples_LDpruned_out.gds.tmp'
  #   rename '316Samples_LDpruned_out.gds.tmp' (565.2M, reduced: 12.4K)
  #   # of fragments: 20
```

Then copy files out of workdir to local computer to run below code in RStudio:

```bash
cp 316Samples_LDpruned_out.gds /lustre2/home/lc736_0001/song_sparrow/final_vcf/pca
```

Now, bring to local computer and plot in R.




