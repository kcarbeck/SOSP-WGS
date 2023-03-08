# Script for prepping mapped reads for variant caller
# authors:katherine carbeck and jen walsh

# A new Bam file will be produced for each step. In the end, you only need to keep samplename_sortedRGmark.bam

# Reserved a medium machine (24 cores)
# Bam files from mapping script above are indexed and sorted already
# copy into workdir:
\ls *.bam | parallel -j 20 cp -a {} /workdir/jlw395 &



#### Prepare the genome: index it (fai and dict files)
# R:Refernece
# O:Output
java -Xmx48g -jar /programs/picard-tools-2.8.2/picard.jar CreateSequenceDictionary R=SongSparrow_reference.fasta O=SongSparrow_reference.dict
samtools faidx SongSparrow_reference.fasta

#Traditionally, I would add readgroup info, but I have SM info in the bam files, which is used for labeling individuals in VCF and am fine with that. I have no other use for RG info for this project

###### mark duplicates - identifies duplicate reads: https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-
# Metrics File: File name to write duplicate metrics to
# MAX_FILE_HANDLES_FOR_READ_ENDS_MAP: maximum number of file handles to keep open when spilling read ends to a desk. keep this set at 1000
# uses a lot of memory (could probably run more than 8 at a time if I decreased -Xmx to ~10g?)

for RGBAM in *_sorted.bam; do

  #string replacement command
  OUT=${RGBAM/%_sorted.bam/_sorted_mark.bam}
  METRICS=${RGBAM/%_sorted.bam/.metrics.txt}

  #create index file
  echo "java -Xmx15g -jar /programs/picard-tools-2.8.2/picard.jar MarkDuplicates INPUT=$RGBAM OUTPUT=$OUT METRICS_FILE=$METRICS MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000" >> /workdir/jlw395/markDuplicatesCommands.txt

done

parallel -j 8 < /workdir/jlw395/markDuplicatesCommands.txt



#### Validate files
# since running HaplotypeCaller, don't need to realign or fix mates unless there is an error in ValidateSamFile

[ -f /workdir/jlw395/validateCommands.txt ] && rm /workdir/jlw395/validateCommands.txt

cd /workdir/jlw395/

for MARKBAM in *mark.bam; do

  OUTFILE=${MARKBAM/%.bam/_validate}

  echo "java -Xmx48g -jar /programs/picard-tools-2.8.2/picard.jar ValidateSamFile I=$MARKBAM OUTPUT=$OUTFILE MODE=SUMMARY" >> /workdir/jlw395/ValidateCommands.txt

done

parallel -j 22 < /workdir/jlw395/ValidateCommands.txt

# concatenate all validate output into one text file to see if there are any errors
cat *validate > summary_validate.txt



#### index .bam files for HaplotypeCaller
#can run multiple at a time
[ -f /workdir/jlw395/IndexCommands.txt ] && rm /workdir/jlw395/IndexCommands.txt

cd /workdir/jlw395/

for MARKBAM in *mark.bam; do

  #create index file
  echo "java -Xmx10g -jar /programs/picard-tools-2.8.2/picard.jar BuildBamIndex I=$MARKBAM" >> /workdir/jlw395/IndexCommands.txt

done

parallel -j 22 < /workdir/jlw395/IndexCommands.txt
