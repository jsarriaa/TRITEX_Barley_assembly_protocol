#!/bin/zsh

# Run the pipe `nxtrim | bfc`  on a set of mate-pair reads produced with the Nextera method.
 
n=(-nxtrim '/filer-dg/agruppen/seq_shared/mascher/source/NxTrim/nxtrim')
f=(-bfc '/filer-dg/agruppen/seq_shared/mascher/source/bfc-master/bfc')
m=(-minlen 100)

zparseopts -D -E -K -- -nxtrim:=n -bfc:=f -minlen:=m -base:=b -genomesize:=g -threads:=t -hash:=c -outdir:=d -lib:=l

minlen=$m[2]
size=$g[2]
threads=$t[2]
count=$c[2]
outdir=$d[2]
lib=$l[2]
base=$b[2]
bfc=$f[2]
nxtrim=$n[2]

{
 $nxtrim --justmp --stdout -l $minlen 2> $outdir/${lib}_nxtrim.err \
  -1 <(find $base/raw_reads/$lib | egrep '_R1|_1.fq.gz$' | grep 'f.*q.gz$' | sort | xargs cat) \
  -2 <(find $base/raw_reads/$lib | egrep '_R2|_2.fq.gz$' | grep 'f.*q.gz$' | sort | xargs cat) \
  | $bfc -s $size -t $threads -r $count /dev/stdin \
  | awk -v o="pigz -p $threads > $outdir/${lib}.mp_unknown_bfc_R1.fq.gz" \
	-v p="pigz -p $threads > $outdir/${lib}.mp_unknown_bfc_R2.fq.gz" \
	'((NR - 1) % 8) + 1 < 5 { print | o; next } {print | p}' 
} 2> $outdir/${lib}_bfc.err > $outdir/${lib}_bfc.out 
