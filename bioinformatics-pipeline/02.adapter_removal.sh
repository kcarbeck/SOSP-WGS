#Script for trimming and adpater removal
#author: katherine carbeck
#Modified September 2022

#Trimns: Trim consecutive Ns from the 5’ and 3’ termini. If quality trimming is also enabled (--trimqualities), then stretches of mixed low-quality bases and/or Ns are trimmed.
#Trimqualities: Trim consecutive stretches of low quality bases (threshold set by --minquality) from the 5’ and 3’ termini. If trimming of Ns is also enabled (--trimns), then stretches of mixed low-quality bases and Ns are trimmed.
#minquality: Set the threshold for trimming low quality bases using --trimqualities and --trimwindows. Default is 2
#collapse: In paired-end mode, merge overlapping mates into a single and recalculate the quality scores. In single-end mode, attempt to identify templates for which the entire sequence is available. In both cases, complete “collapsed” reads are written with a ‘M_’ name prefix, and “collapsed” reads which are trimmed due to quality settings are written with a ‘MT_’ name prefix.
#adapterlist: Read one or more adapter sequences from a table. The first two columns (separated by whitespace) of each line in the file are expected to correspond to values passed to –adapter1 and –adapter2.
#minlength: Reads shorter than this length are discarded following trimming. Defaults to 15Reads shorter than this length are discarded following trimming. Defaults to 15.

# Usage: /programs/adapterremoval_2.1.1/bin/AdapterRemoval --file1 Pathto/R1.fastq.gz --file2 Pathto/R2.fastq.gz --trimns --trimqualities --minquality 20 --minlength 25 --collapse --threads 8 --adapter-list for_adapter_removal_2.txt --basename LIBRARY_SAMPLEID


# create adapterRemovalCommands.txt (have to manually change "SAMPLENAME" to what you want to name each sample)
awk -v FS=" " '
{
    if (NR > 0){
        x=$0;
        getline y;
        print "/programs/adapterremoval_2.1.1/bin/AdapterRemoval --file1" " " x " " "--file2" " " y " " "--trimns --trimqualities --minquality 35 --minlength 25 --collapse --threads 8 --adapter-list for_adapter_removal.txt --basename" " " "SAMPLENAME"
    }
}
' listOfFiles.list > adapterRemovalCommands.txt




# output (adapterRemovalCommands.txt) looks like this:
/programs/adapterremoval_2.1.1/bin/AdapterRemoval --file1 S23_CKDL220033855-1A_HMLFCDSX5_L2_1.fq.gz --file2 S23_CKDL220033855-1A_HMLFCDSX5_L2_2.fq.gz --trimns --trimqualities --minquality 35 --minlength 25 --collapse --threads 8 --adapter-list for_adapter_removal.txt --basename SAMPLENAME
/programs/adapterremoval_2.1.1/bin/AdapterRemoval --file1 S49_CKDL230002535-1A_HMLFCDSX5_L2_1.fq.gz --file2 S49_CKDL230002535-1A_HMLFCDSX5_L2_2.fq.gz --trimns --trimqualities --minquality 35 --minlength 25 --collapse --threads 8 --adapter-list for_adapter_removal.txt --basename SAMPLENAME
/programs/adapterremoval_2.1.1/bin/AdapterRemoval --file1 S77_CKDL220033857-1A_HMLFCDSX5_L2_1.fq.gz --file2 S77_CKDL220033857-1A_HMLFCDSX5_L2_2.fq.gz --trimns --trimqualities --minquality 35 --minlength 25 --collapse --threads 8 --adapter-list for_adapter_removal.txt --basename SAMPLENAME
/programs/adapterremoval_2.1.1/bin/AdapterRemoval --file1 S81_CKDL230002536-1A_HMLFCDSX5_L2_1.fq.gz --file2 S81_CKDL230002536-1A_HMLFCDSX5_L2_2.fq.gz --trimns --trimqualities --minquality 35 --minlength 25 --collapse --threads 8 --adapter-list for_adapter_removal.txt --basename SAMPLENAME
/programs/adapterremoval_2.1.1/bin/AdapterRemoval --file1 S83_CKDL220033859-1A_HMLFCDSX5_L2_1.fq.gz --file2 S83_CKDL220033859-1A_HMLFCDSX5_L2_2.fq.gz --trimns --trimqualities --minquality 35 --minlength 25 --collapse --threads 8 --adapter-list for_adapter_removal.txt --basename SAMPLENAME
/programs/adapterremoval_2.1.1/bin/AdapterRemoval --file1 S88_CKDL220033860-1A_HMLFCDSX5_L2_1.fq.gz --file2 S88_CKDL220033860-1A_HMLFCDSX5_L2_2.fq.gz --trimns --trimqualities --minquality 35 --minlength 25 --collapse --threads 8 --adapter-list for_adapter_removal.txt --basename SAMPLENAME

# then edit SAMPLENAME field:
/programs/adapterremoval_2.1.1/bin/AdapterRemoval --file1 S23_CKDL220033855-1A_HMLFCDSX5_L2_1.fq.gz --file2 S23_CKDL220033855-1A_HMLFCDSX5_L2_2.fq.gz --trimns --trimqualities --minquality 35 --minlength 25 --collapse --threads 8 --adapter-list for_adapter_removal.txt --basename AK_kenaiensis_S23
/programs/adapterremoval_2.1.1/bin/AdapterRemoval --file1 S49_CKDL230002535-1A_HMLFCDSX5_L2_1.fq.gz --file2 S49_CKDL230002535-1A_HMLFCDSX5_L2_2.fq.gz --trimns --trimqualities --minquality 35 --minlength 25 --collapse --threads 8 --adapter-list for_adapter_removal.txt --basename OR_montana_S49
/programs/adapterremoval_2.1.1/bin/AdapterRemoval --file1 S77_CKDL220033857-1A_HMLFCDSX5_L2_1.fq.gz --file2 S77_CKDL220033857-1A_HMLFCDSX5_L2_2.fq.gz --trimns --trimqualities --minquality 35 --minlength 25 --collapse --threads 8 --adapter-list for_adapter_removal.txt --basename ON_melodia_S77
/programs/adapterremoval_2.1.1/bin/AdapterRemoval --file1 S81_CKDL230002536-1A_HMLFCDSX5_L2_1.fq.gz --file2 S81_CKDL230002536-1A_HMLFCDSX5_L2_2.fq.gz --trimns --trimqualities --minquality 35 --minlength 25 --collapse --threads 8 --adapter-list for_adapter_removal.txt --basename ON_melodia_81
/programs/adapterremoval_2.1.1/bin/AdapterRemoval --file1 S83_CKDL220033859-1A_HMLFCDSX5_L2_1.fq.gz --file2 S83_CKDL220033859-1A_HMLFCDSX5_L2_2.fq.gz --trimns --trimqualities --minquality 35 --minlength 25 --collapse --threads 8 --adapter-list for_adapter_removal.txt --basename ON_melodia_83
/programs/adapterremoval_2.1.1/bin/AdapterRemoval --file1 S88_CKDL220033860-1A_HMLFCDSX5_L2_1.fq.gz --file2 S88_CKDL220033860-1A_HMLFCDSX5_L2_2.fq.gz --trimns --trimqualities --minquality 35 --minlength 25 --collapse --threads 8 --adapter-list for_adapter_removal.txt --basename ON_melodia_88

# move updated txt file and "for_adapter_removal.txt" to server workdir and run in parallel 
parallel -j 22 < /workdir/kcarbeck/fastqc/adapterRemovalCommands.txt
    # Read 2 adapters / adapter pairs from 'for_adapter_removal.txt'...
    # Trimming paired end reads ...
    # Processed a total of 12,750,848 reads in 1:31.5s; 139,000 reads per second on average ...
    # Read 2 adapters / adapter pairs from 'for_adapter_removal.txt'...
    # Trimming paired end reads ...
    # Processed a total of 6,541,312 reads in 1:42.0s; 64,000 reads per second on average ...
    # Read 2 adapters / adapter pairs from 'for_adapter_removal.txt'...
    # Trimming paired end reads ...
    # Processed a total of 8,765,440 reads in 2:25.1s; 60,000 reads per second on average ...
    # Read 2 adapters / adapter pairs from 'for_adapter_removal.txt'...
    # Trimming paired end reads ...
    # Processed a total of 69,265,594 reads in 5:17.7s; 218,000 reads per second on average ...
    # Read 2 adapters / adapter pairs from 'for_adapter_removal.txt'...
    # Trimming paired end reads ...
    # Processed a total of 50,254,334 reads in 5:31.4s; 151,000 reads per second on average ...
    # Read 2 adapters / adapter pairs from 'for_adapter_removal.txt'...
    # Trimming paired end reads ...
    # Processed a total of 67,511,438 reads in 6:12.0s; 181,000 reads per second on average ...
