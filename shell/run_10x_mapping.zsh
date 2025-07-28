#!/bin/zsh

# Align 10X reads to an assembly and assign reads to molecules

a=(-adapter 'AGATCGGAAGAGC')
m=(-minlen_read 30)
q=(-minq 20)
d=(-maxdist 500000)
n=(-minlen_mol 1000)
x=(-maxlen 800)
o=(-onlymol 0)
cut=(-cutadapt '/opt/Bio/cutadapt/1.15/bin/cutadapt')
mini=(-minimap '/opt/Bio/minimap2/2.14/bin/minimap2')
sam=(-samtools '/opt/Bio/samtools/1.9/bin/samtools')
stk=(-seqtk '/opt/Bio/seqtk/1.0-r76/bin/seqtk')
bed=(-bedtools '/opt/Bio/bedtools/2.26.0/bin/bedtools')
novo=(-novosort '/opt/Bio/novocraft/V3.06.05/bin/novosort')

zparseopts -D -K -- -adapter:=a -minlen_read:=m -threads:=t -mem:=e -ref:=r -tmp:=p \
 -minq:=q -maxdist:=d -minlen_mol:=n -maxlen:=x -onlymol:=o \
 -cutadapt:=cut -minimap:=mini -samtools:=sam -seqtk:=stk \
 -bedtools:=bed -novosort:=novo

novosort=$novo[2]
bedtools=$bed[2]
seqtk=$stk[2]
samtools=$sam[2]
minimap=$mini[2]
cutadapt=$cut[2]
adapter=$a[2]
minlen=$m[2]
minimap_threads=$t[2]
mem=$e[2]
ref=$r[2]
tmp=$p[2]
minq=$q[2]
maxdist=$d[2]
minlenm=$n[2]
maxlen=$x[2]
onlymol=$o[2]
i=$1

if [[ $minimap_threads -gt 16 ]]; then
 novosort_threads=16
else
 novosort_threads=$minimap_threads
fi

awkcmd='
NR % 8 == 1 {
 split($1, h, ":")
 next
}

NR % 8 == 2 {
 head="@"substr($1, 1, 16)":"h[2]":"h[3]":"h[4]":"h[5]":"h[6]":"h[7]
 print head"\n"substr($1, 24)
 next
}

NR % 8 == 4 {
 print substr($1, 24)
 next
}

NR % 8 == 5 {
 next
}

NR % 8 == 6 {
 print head
}

1
'

base=$i/$i:t
bam=${base}.bam
stats=${base}_mapping_stats.tsv

minimaperr=${base}_minimap.err
samtoolserr=${base}_samtools.err
sorterr=${base}_novosort.err
sorterr2=${base}_novosort2.err
cutadapterr=${base}_cutadapt.err
statserr=${base}_mapping_stats.err

b=$base:t
rgentry="@RG\tID:$b\tPL:ILLUMINA\tPU:$b\tSM:$b"

# Note the regular expression for find R1 and R2. It will not work if, for instance, sample names contain R1 and R2.
if [[ $onlymol -eq 0 ]]; then
 $seqtk mergepe <(cat ${i}/*_R1*f*q.gz) <(cat ${i}/*_R2*f*q.gz) \
  | awk "$awkcmd" /dev/stdin \
  | $cutadapt -f fastq --interleaved -a $adapter -A $adapter -m $minlen /dev/stdin 2> $cutadapterr \
  | $minimap -ax sr -R $rgentry -I '50G' -t $minimap_threads -2 -K '5G' $ref /dev/stdin 2> $minimaperr \
  | $samtools view -Su /dev/stdin 2> $samtoolserr \
  | $novosort -c $novosort_threads -t $tmp -m $mem --md --keepTags /dev/stdin 2> $sorterr \
  | $novosort -c $novosort_threads -t $tmp -m $mem -n -o $bam /dev/stdin 2> $sorterr2 

 echo $pipestatus > ${base}_pipestatus.txt
fi

out=${bam:r}_molecules.tsv.gz
err=${bam:r}_molecules.err
ncores=$minimap_threads

{
 $samtools view -f 2 -uF 3332 $bam \
  | $bedtools pairtobed -f 1 -type both -abam - -bedpe -b =(awk '{print $1"\t0\t"$2"\t"$1}' $ref.fai) \
  | awk 'old != $7 {old=$7; printf "\n"$0; next} {printf "\t"$0}' \
  | awk -v minq=$minq 'NF && $18 == $25 && $8 >= minq && $22 >= minq' \
  | cut -f 2,3,5,6,7,12,14 \
  | sed 's/\(\t[^\t:]\+\):[^\t]\+/\1/' \
  | awk '{ printf $7"\t"
	   if($1 <= $3){
	    printf $1 - $6"\t"
	   } else {
	    printf $3 - $6"\t"
	   }
	   if($2 >= $4){
	    printf $2 - $6"\t"
	   } else {
	    printf $4 - $6"\t"
	   }
	   print $5}' \
  | awk -v maxlen=$maxlen 'maxlen >= $3 - $2' \
  | sort --parallel=$ncores --buffer-size=$mem -k 4,4 -k 1,1 -k 2,2n -k 3,3n \
  | awk -v maxdist=$maxdist -v minlen=$minlenm 'NR == 1 {s=$1; b=$4; e=$3; start=-1; end=0}
	 s != $1 || b != $4 || $2 - e > maxdist {
	  if(n > 1 && end - start >= minlen){
	   print s"\t"start"\t"end"\t"b"\t"n
	  }
	  s=$1; b=$4; e=$3; n=0; start=-1; end=0
	 }
	 {
	  n++
	  if(start < 0 || $2 < start){
	   start = $2
	  }
	  if($3 > end){
	   end = $3
	  }
	 }
	 END {
	  print s"\t"start"\t"end"\t"b"\t"n
  	 }' \
  | pigz -p $ncores > $out
} 2> $err
