# Bioinformatics pipeline for song sparrow whole genome re-sequencing

## 1. FastQC & MultiQC

Run fastqc.sh script on terminal and view output using: 

    export LC_ALL=en_US.UTF-8
    export PATH=/programs/miniconda3/bin:$PATH
    source activate multiqc

run software using command: 

    multiqc .

after done running deactivate conda 
    
    conda deactivate
    
there will be a file called multiqc_report.html. Copy this to your computer and visualize in web browser.


## 2. Adapter removal

Resources: https://adapterremoval.readthedocs.io/en/latest/
Citation: Schubert, Lindgreen, and Orlando (2016). AdapterRemoval v2: rapid adapter trimming, identification, and read merging. BMC Research Notes, 12;9(1):88

## 3. Mapping

- I did this on medium 40 core machine at Cornell BioHPC. bowtie2 step completed about 19 individuals per day 
- samtools sort was the limiting factor – increasing RAM makes it faster
- there were issues with some files not finishing sorting -- would have to kill them and restart...still not sure why this happens. I think it's easier to run the code in 3 separate steps (first align, then index, then qualimap)

## 4. Prep for variant calling
 I did this on a medium 24 core machine at Cornell BioHPC.
 
Add or Replace Read Groups:
https://gatk.broadinstitute.org/hc/en-us/articles/360037226472-AddOrReplaceReadGroups-Picard-
- the awk script included was run on personal computer to assemble script for add or replace read groups and put in a text file to run in parallel on terminal

Mark Duplicates:
https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-
- uses a lot of memory (could probably run more than 8 at a time if I decreased -Xmx to ~10g?)

## 5. Variant discovery & filtering

split by chromosome/contig to run in parallel

mpileup: http://www.htslib.org/doc/bcftools.html#mpileup





