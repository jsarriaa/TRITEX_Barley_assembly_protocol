#!/bin/zsh

minimap='/opt/Bio/minimap2/2.17/bin/minimap2'
samtools='/opt/Bio/samtools/1.8/bin/samtools'
novosort="/opt/Bio/novocraft/V3.06.05/bin/novosort"
cutadapt='/opt/Bio/cutadapt/1.15/bin/cutadapt'

indexsize='50G'
batchmem='5G'

m=(-minlen_read 30)

zparseopts -D -K -- -minlen_read:=m -threads:=t -mem:=e -ref:=r -tmp:=p 

minlen=$m[2]
threads=$t[2]
mem=$e[2]
ref=$r[2]
tmp=$p[2]
i=$1

if [[ $threads -gt 16 ]]; then
 novosort_threads=16
else
 novosort_threads=$threads
fi

base=$i/$i:t
bam=${base}.bam
minimaperr=${base}_minimap.err
samtoolserr=${base}_samtools.err
sorterr1=${base}_novosort1.err
sorterr2=${base}_novosort2.err
cutadapterr=${base}_cutadapt.err

b=$base:t
rgentry="@RG\tID:$b\tPL:ILLUMINA\tPU:$b\tSM:$b"

#AGATCGGAAGAGC is the common the sequence of the Illumina read1 and read2 adapters. It need not be changed by the users.
#https://cutadapt.readthedocs.io/en/stable/guide.html#illumina-truseq
#The `find | egrep | grep` command need to be modified to match the format of the file name
$cutadapt -f fastq --interleaved -a AGATCGGAAGAGC \
  -O 1 -m $minlen 2> $cutadapterr \
  <(find $i | egrep '_R1' | grep 'f.*q.gz$' | sort | xargs zcat) \
  <(find $i | egrep '_R2' | grep 'f.*q.gz$' | sort | xargs zcat) \
 | $minimap -ax sr -R $rgentry -t $threads -2 -I $indexsize -K$batchmem $ref /dev/stdin 2> $minimaperr \
 | $samtools view -Su /dev/stdin 2> $samtoolserr \
 | $novosort -c 16 -t $tmp -m $mem --keepTags --md /dev/stdin 2> $sorterr1 \
 | $novosort -c 16 -t $tmp -m $mem -n -o $bam /dev/stdin 2> $sorterr2

echo $pipestatus > ${base}_pipestatus.txt

