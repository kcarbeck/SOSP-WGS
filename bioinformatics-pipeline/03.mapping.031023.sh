#samtools align to ref genome new batch of samples
#10 march 2023
#authors:katherine carbeck & jen walsh

#load your sample names into an array.
#in this example I will get the sample names from the adapter removal files I generated after a previous step.

bowtie2-build -f /workdir/kcarbeck/SongSparrow_reference.fasta SOSPindex

## this was done on medium 24 core machine --
mkdir /workdir/kcarbeck/align

INDS=($(for i in /workdir/kcarbeck/fastqc/*.settings; do echo $(basename -s .settings $i); done))
#basename - will remove the directory path and returns the file name. -s tell it which suffix to remove from the end of the file name (in this case .settings)
#note: the variable INDS will now contain an array of the sample names extracted from the files names.

#If you want to see what is stored in this variable you can type:
#(@=all the elements in the array)
echo ${INDS[@]}
# AK_kenaiensis_S23 ON_melodia_81 ON_melodia_83 ON_melodia_88 ON_melodia_S77 OR_montana_S49




###
REFERENCE=SOSPindex

# Modified this loop to map files and pipe to a bam file
[ -f /workdir/kcarbeck/align/bowtie2Commands.txt ] && rm /workdir/kcarbeck/align/bowtie2Commands.txt

for SAMPLEID in ${INDS[@]}; do
  #declare variables. This makes it easier and neater to write your command line and you just have to change these for future projects.
  ONESEQ=/workdir/kcarbeck/fastqc/${SAMPLEID}.pair1.truncated
  TWOSEQ=/workdir/kcarbeck/fastqc/${SAMPLEID}.pair2.truncated
  USEQ=/workdir/kcarbeck/fastqc/${SAMPLEID}.collapsed,/workdir/kcarbeck/fastqc/${SAMPLEID}.collapsed.truncated,/workdir/kcarbeck/fastqc/${SAMPLEID}.singleton.truncated
  OUTPUT=/workdir/kcarbeck/align/${SAMPLEID}.bam

  # align with bowtie - the output is piped directly into samtools to avoid having the intermediate .sam file.
  echo "Aligning $SAMPLEID with bowtie"
  #this just writes a line telling you which sample is being worked on.
  echo "bowtie2 -p 6 --phred33 --very-sensitive-local -x $REFERENCE -I 149 -X 900 --rg-id $SAMPLEID --rg SM:$SAMPLEID -1 $ONESEQ -2 $TWOSEQ -U $USEQ|samtools view -bS > $OUTPUT" >> /workdir/kcarbeck/align/bowtie2Commands.txt

done

parallel -j 4 < /workdir/kcarbeck/align/bowtie2Commands.txt


### sambamba to sort
#https://www.basepairtech.com/blog/sorting-bam-files-samtools-vs-sambamba/
# loop for Sambamba sort
[ -f /workdir/kcarbeck/align/SambambaSortCommands.txt ] && rm /workdir/kcarbeck/align/SambambaSortCommands.txt

cd /workdir/kcarbeck/align/

for BAM in *.bam; do
  OUT=/workdir/kcarbeck/align/${BAM}_sorted.bam

  #create sort file
  echo "/programs/sambamba-0.7.1/sambamba sort -t 30 -m 45G -o $OUT $BAM" >> /workdir/kcarbeck/align/SambambaSortCommands.txt

done

parallel -j 4 < /workdir/kcarbeck/align/SambambaSortCommands.txt


# check if bam files are sorted 
# Will show SO:coordinate if sorted
samtools view -H *.bam | grep @HD


### index for loop
[ -f /workdir/kcarbeck/align/samtoolsIndexCommands.txt ] && rm /workdir/kcarbeck/align/samtoolsIndexCommands.txt

cd /workdir/kcarbeck/align/

for BAM in *.bam; do

  #create index file
  echo "samtools index $BAM" >> /workdir/kcarbeck/align/samtoolsIndexCommands.txt

done

# only 6 files in this batch
parallel -j 6 < /workdir/kcarbeck/align/samtoolsIndexCommands.txt

### qualimap for loop
# may have to increase mem limit:
#/programs/qualimap_v2.2.1/qualimap bamqc -bam SAMPLE.bam --java-mem-size=30G -outfile SAMPLE.sorted.pdf

[ -f /workdir/kcarbeck/align/qualimapCommands.txt ] && rm /workdir/kcarbeck/align/qualimapCommands.txt

cd /workdir/kcarbeck/align/

for BAM in *.bam; do

  #string replacement command
  PDF=${BAM/%.bam/.pdf}
  # qualimap
  echo "/programs/qualimap_v2.2.1/qualimap bamqc -bam $BAM -outfile $PDF --java-mem-size=128G" >> /workdir/kcarbeck/align/qualimapCommands.txt

done

parallel -j 4 < /workdir/kcarbeck/align/qualimapCommands.txt


# to run and get log output:
# bash qualimap.sh 2>&1 | tee qualimap_may1_$(date +%Y%m%d-%Hh%Mm%Ss).log

