#!/bin/bash
#katherine carbeck
#28 feb 2023
#script to convert vcf to baypass input
# following https://github.com/Bio-protocol/Genome-Environment_Association_Analyses_with_BayPass


# baypass requires:
    # -gfile
    # -efile
## The genotyping data file is organized as a matrix with nsnp rows and 2 * npop columns. Each row corresponds to one marker and the number of columns is twice the number of populations because each pair of numbers corresponds to each allele (or read counts for PoolSeq experiment) counts in one population 

#what samples are in the VCF
bcftools query -l filtered_SOSP_352Samples_122322.LDpruned.vcf.gz
bcftools query -l filtered_SOSP_352Samples_122322.LDpruned.vcf.gz | wc -l
    # 352

#remove single ON sample 
vcftools --gzvcf filtered_SOSP_352Samples_122322.LDpruned.vcf.gz --remove-indv ON_melodia_S79 --recode --recode-INFO-all 
bgzip -c out.recode.vcf bgzip > filtered_SOSP_351Samples.LDpruned.vcf.gz

#now how many in vcf?
bcftools query -l filtered_SOSP_351Samples.LDpruned.vcf.gz | wc -l
    #351

########### Step 1: generate genotyping data file from the full VCF using vcf2baypass.pl ###########
zcat /workdir/kcarbeck/filtered_SOSP_351Samples.LDpruned.vcf.gz | perl vcf2baypass.pl /workdir/kcarbeck/popfile.tsv /workdir/kcarbeck/cache/SOSP.baypass 

    # This step will generate 3 files in cache/ with the SOSP.baypass prefix. The main file cache/SOSP.baypass.txt is the one required by BayPass. Two other files cache/SOSP.baypass.pop and cache/SOSP.baypass.snp record the order of populations and locations of SNPs, respectively, and will be used in the following steps. Check if the numbers of populations and SNPs are correct:
wc -l SOSP.baypass.pop # number of populations, which should be 27 in this case
wc -l SOSP.baypass.txt # number of SNPs, which should be 13,663,960
wc -l SOSP.baypass.snp # number of SNPs as well 13663960


############# Step 2: generate the covariate data file from environmental dataset using env2baypass.pl ###########
cat popenvPCA.tsv | perl env2baypass.pl SOSP.baypass.pop SOSP_env.baypass

    # The output covariate data file cache/SOSP_env.baypass.txt contains the values of the environmental variables for each population in the format required by BayPass. The associated file cache/SOSP_env.baypass.cov records the names of the covariates, which will be used for plotting.


########## Step 3: running BayPass ##########
#download and compile BayPass
cd program/
wget http://www1.montpellier.inra.fr/CBGP/software/baypass/files/baypass_2.3.tar.gz
tar -zxvf baypass_2.3.tar.gz
cd baypass_2.3/sources/
make clean all FC=gfortran
/workdir/kcarbeck/program/baypass_2.3/sources/g_baypass -help # check
cd ..

#run BayPass under the core model mode to generate covariance matrix
npop=$(wc -l /workdir/kcarbeck/cache/SOSP.baypass.pop | cut -d " " -f1)
# echo $npop #27

/workdir/kcarbeck/program/baypass_2.3/sources/g_baypass -npop $npop -gfile /workdir/kcarbeck/cache/SOSP.baypass.txt -outprefix /workdir/kcarbeck/cache/output/SOSP_core -nthreads 16
    # BayPass Version 2.3
    # Reading and checking input data
    # Reading and checking analysis parameters
    # Analysis started. It will consist of:
    #     i) 20 pilot runs of  500 iterations (to adjust proposal distributions)
    #    ii) a burn-in period of  5000 iterations
    #   iii) final MCMC sampling of 1000 parameter values sampled every  20 iterations (i.e.,   20000 iterations)
    #  Note: progress bars indicate the progression at each of the three steps
    #        while the given e.t.a. is the (estimated) remaining time until full completion of the analysis
    # estimated about 7-8 days to complete







