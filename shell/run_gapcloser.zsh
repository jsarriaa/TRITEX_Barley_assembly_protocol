#!/bin/zsh

# Fill gaps with GapCloser.

g=(-gapcloser '/filer-dg/agruppen/seq_shared/mascher/source/GapCloser-bin-v1.12-r6/GapCloser')
s=(-samtools '/opt/Bio/samtools/1.7/bin/samtools')
e=(-seqtk '/opt/Bio/seqtk/1.0-r76/bin/seqtk')

t=(-threads 1)
k=(-kmer 31)
l=(-length 155)
n=(-name "")

zparseopts -D -E -K -- -kmer:=k -length:=l -config:=c -input:=i -threads:=t -outdir:=o -name:=n \
 -gapcloser:=g -samtools:=s -seqtk:=e

gap=$g[2]
samtools=$s[2]
seqtk=$e[2]
outdir=$o[2]
length=$l[2]
kmer=$k[2]
input=$i[2]
cfg=$c[2]
threads=$t[2]

if [[ -z $n[2] ]]; then
 dir=$outdir/${input:t:r}
else
 dir=$outdir/$n[2]
fi

fasta=$dir/${dir:t}.fasta

if ! mkdir $dir; then
 exit 1
fi

$gap -t $threads -l $length -p $kmer -a $input -b $cfg -o $fasta > $dir/log.out 2> $dir/log.err 

$samtools faidx $fasta 

# Calculate the N content

$seqtk cutN -n 1 -g $fasta | awk '{s+=$3-$2} END {print s}' | read a
awk '{s+=$2} END {print s"."}' $fasta.fai | read b
echo $a $b $[ $a / $b ] >> $dir/log.out
