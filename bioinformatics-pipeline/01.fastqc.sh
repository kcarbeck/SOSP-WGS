# fastqc 
# author: katherine carbeck

# find /workdir/kcarbeck/usftp21.novogene.com/01.RawData/*/*.gz -type f > listOfFiles.list

find *.gz -type f > listOfFiles.list

awk -v FS="\t" '
{
    if (NR > 1){
        samp=$1
        print "fastqc"  " " samp
    }
}
' listOfFiles.list > list.txt


parallel -j 22 < /workdir/kcarbeck/list.txt
