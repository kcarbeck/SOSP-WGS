# fastqc 
# author: katherine carbeck


# cp all files in directory and subdirs without dir structure using GNU coreutils:
find /lustre2/home/lc736_0001/song_sparrow/rawdata/ontario_feb23/usftp21.novogene.com/01.RawData -type f -exec cp -t /workdir/kcarbeck/fastqc {} +

# count files in current directory
ls -1 | wc -l

# find /workdir/kcarbeck/usftp21.novogene.com/01.RawData/*/*.gz -type f > listOfFiles.list

find *.gz -type f > listOfFiles.list

awk -v FS="\t" '
{
    if (NR > 0){
        samp=$1
        print "fastqc"  " " samp
    }
}
' listOfFiles.list > list.txt


parallel -j 22 < /workdir/kcarbeck/fastqc/list.txt


# Run fastqc.sh script on terminal and view output using: 
export LC_ALL=en_US.UTF-8
export PATH=/programs/miniconda3/bin:$PATH
source activate multiqc

# run software using command: 
multiqc .

# after done running deactivate conda 
conda deactivate
    
# there will be a file called multiqc_report.html. Copy this to your computer and visualize in web browser.
