#!/bin/zsh

# Run BFC on a set of overlapping reads. Despite the name, the scripts makes no assumption about the insert size, so it can be used for overlapping paired-end libraries of various insert sizes.

b=(-bfc '/filer-dg/agruppen/seq_shared/mascher/source/bfc-master/bfc')
k=(-ktrim 61)
m=(-bbmerge 1)
t=(-threads 1)
d=(-outdir '.')

zparseopts -D -E -K -- -bbmerge:=m -bfc:=b -ktrim:=k -input:=i -genomesize:=g -threads:=t -outdir:=d 

bfc=$b[2]
size=$g[2]
threads=$t[2]
prefix=$d[2]"/PE450_bfc"
folder=$i[2]
ktrim=$k[2]

if [[ $m[2] -eq 0 ]]; then
 regex='_pear.assembled.fastq$'
else
 regex='_bbmerge.fq.gz$'
fi

{ 
 # Create k-mer hash to be re-used for PE800 and MP correction
 find $folder -type f | sort | grep $regex | xargs cat \
  | $bfc -s $size -t $threads -Ed $prefix.hash /dev/stdin 
} 2> ${prefix}_hash.err && \
{ 
 # Correct PE450 reads
 find $folder -type f | sort | grep $regex | xargs cat \
  | $bfc -s $size -t $threads -r $prefix.hash /dev/stdin \
  | tee >(sed -n '2~4p' | awk '{print length}' | sort -n | uniq -c | awk '{print $2"\t"$1}' > ${prefix}_correct_length_dist.tsv) \
  | pigz -p $threads > ${prefix}_correct.fq.gz
} 2> ${prefix}_correct.err && \
{
 #  Trim corrected reads containing singleton k-mers
 $bfc -1 -s $size -k $ktrim -t $threads ${prefix}_correct.fq.gz \
  | tee >(sed -n '2~4p' | awk '{print length}' | sort -n | uniq -c | awk '{print $2"\t"$1}' > ${prefix}_correct_trim_length_dist.tsv) \
  | pigz -p $threads > ${prefix}_correct_trim.fq.gz 
} 2> ${prefix}_trim.err
