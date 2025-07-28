#!/bin/zsh

# Run the pipe `cutadapt | bfc`  on a set of mate-pair reads produced with the Nextera method.

f=(-bfc '/filer-dg/agruppen/seq_shared/mascher/source/bfc-master/bfc')
u=(-cutadapt '/opt/Bio/cutadapt/1.15/bin/cutadapt')
a=(-adapter 'AGATCGGAAGAGC') 
# the common sequence of the Illumina 5' and 3' adapters is used,
# see https://cutadapt.readthedocs.io/en/stable/guide.html#illumina-truseq
m=(-minlen 100)
# discard reads shorter than 100 bp after adapter trimming
t=(-threads 1)

zparseopts -D -E -K -- -cutadapt:=u -bfc:=f -adapter:=a -minlen:=m -base:=b -genomesize:=g -threads:=t -hash:=c -outdir:=d -lib:=l

adapter=$a[2]
minlen=$m[2]
size=$g[2]
threads=$t[2]
count=$c[2]
outdir=$d[2]
base=$b[2]
bfc=$f[2]
cutadapt=$u[2]

lib='PE800'

{
 $cutadapt -f fastq -a $adapter -A $adapter --interleaved -m $minlen 2> $outdir/${lib}_cutadapt.err \
  <(find $base/raw_reads/$lib | egrep '_R1|_1.fq.gz$' | grep 'f.*q.gz$' | sort | xargs zcat) \
  <(find $base/raw_reads/$lib | egrep '_R2|_2.fq.gz$' | grep 'f.*q.gz$' | sort | xargs zcat) \
  | $bfc -s $size -t $threads -r $count /dev/stdin \
  | awk -v o="pigz -p $threads > $outdir/${lib}.bfc_R1.fq.gz" \
	-v p="pigz -p $threads > $outdir/${lib}.bfc_R2.fq.gz" \
	'((NR - 1) % 8) + 1 < 5 { print | o; next } {print | p}' 
} 2> $outdir/${lib}_bfc.err > $outdir/${lib}_bfc.out 
