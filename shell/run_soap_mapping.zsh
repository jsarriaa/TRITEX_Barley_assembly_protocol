#!/bin/zsh

# Prepare contigs assembled by Minia for use in SOAPDenovo and align paired-end and mate-pair reads to the assembly.

f=(-fusion '/filer-dg/agruppen/seq_shared/mascher/source/SOAPdenovo2-master/fusion/SOAPdenovo-fusion')
s=(-soap '/filer-dg/agruppen/seq_shared/mascher/source/SOAPdenovo2-bin-LINUX-generic-r240/SOAPdenovo-127mer')
t=(-threads 1)
u=(-suffix "")
k=(-kmer 127) # This is the maximum k-mer size in SOAPDenovo.

zparseopts -D -E -K -- -soap:=s -fusion:=f -kmer:=k -config:=c -input:=i -threads:=t -outdir:=o -suffix:=u

input=$i[2]
threads=$t[2]
out=$o[2]
suffix=$u[2]
cfg=$c[2]
kmer=$k[2]
soap=$s[2]
fusion=$f[2]

if [[ -n $suffix ]]; then
 suffix="_"$suffix
fi

if grep -qm1 'gz$' <<< $input; then
 dir=$out/${input:t:r:r}$suffix
else
 dir=$out/${input:t:r}$suffix
fi

prefix="$dir/$dir:t"

if ! mkdir $dir; then
 exit 1
fi

if grep -qm1 'gz$' <<< $input; then
 fasta="$dir/${input:t:r}"
 gzip -cd $input > $fasta
else
 fasta="$dir/${input:t}"
 ln -s $input $dir 
fi

$fusion -D -K $kmer -p $threads -c $fasta -g $prefix > ${prefix}_fusion.out 2> ${prefix}_fusion.err && \
$soap map -s $cfg -g $prefix -p $threads -k $kmer > ${prefix}_map.out 2> ${prefix}_map.err 

if grep -qm1 'gz$' <<< $input; then
 rm -f $fasta
fi
