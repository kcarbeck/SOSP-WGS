# pop gen analyses

## 1. filter for LD 
Ran code contained in 01.pca.sh script using filtered VCF of all individuals. It will output two files needed to run PCA:
* .eigenval and .eigenvec 

And several files needed for downstream analyses:
* .bed = binary file necessary for admixture analysis. It is essentially the genotypes of the pruned dataset recoded as 1s and 0s.
* .bim = a map file (i.e. information file) of the variants contained in the bed file.
* .fam = a map file for the individuals contained in the bed file.

## 2. PCA
run 01.pca.sh then copy eigenval and eigenvec files into R to plot using plotPCA.R

## 3. ADMIXTURE
run 02.admixture.sh then copy .P and .Q files over to plot in R using plotAdmixture.R

## 4. LD Decay
run 03.ldDecay.sh using supporting script 03.ldDecayCalc.py 

