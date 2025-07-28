#!/bin/zsh

# Call Minia3 iteratively on the error-corrected PE450 reads and the assembly from the previous step. 
# The output are the unitig without any heuristic resultion of bulges.

k=(-kmers 100,200,300,350,400,450,500)
m=(-mem 50000)
s=(-restart 0)
e=(-minia '/filer-dg/agruppen/seq_shared/mascher/source/minia/build500/bin/minia')
i=(-scripts 0)
sam=(-samtools '/opt/Bio/samtools/1.9/bin/samtools')
n=(-n50 '/filer-dg/agruppen/seq_shared/mascher/code_repositories/triticeae.bitbucket.io/shell/n50')

zparseopts -D -K -- -minia:=e -kmers:=k \
 -reads:=r -mem:=m -threads:=t -outdir:=d -restart:=s \
 -samtools:=sam -n50:=n

samtools=$sam[2]
n50=$n[2]
kmers=$k[2]
input=$r[2]
mem=$m[2]
threads=$t[2]
dir=$d[2]
restart=$s[2]
minia="${e[2]} -no-bulge-removal -no-tip-removal -no-ec-removal -out-compress 9 -debloom original"

tr , '\n' <<< $kmers | read k 
prefix=$dir/minia_k${k}

if [[ $restart -eq 0 ]]; then

 if [[ -e $prefix ]]; then
  exit 1
 fi

 rm -rf $prefix
 mkdir -p $prefix

 wd=$PWD
 cd $prefix
 eval $minia -in $input -kmer-size $k -max-memory $mem \
  -out-dir $prefix -out $prefix:t -nb-cores $threads \
  > $prefix/log.out 2> $prefix/log.err 

 rm -f *glue* dummy* *h5 *contigs.fa
 fa=$prefix/$prefix:t.unitigs.fa
 $samtools faidx $fa && $n50 $fa.fai > $fa.n50
 pigz -p $threads $fa

 cd $wd
else
 fa=$prefix/${prefix:t}.unitigs.fa
fi

kk=$k

tr , '\n' <<< $kmers | tail -n +2 | while read k; do
 prefix=$dir/minia_k${k}

 if [[ -e $prefix ]]; then
  exit 1
 fi
 
 rm -rf $prefix
 mkdir -p $prefix

 list=$prefix/input_files.txt
 print -l $fa.gz $fa.gz $fa.gz $input > $list

 wd=$PWD
 cd $prefix
 eval $minia -in $list -kmer-size $k -max-memory $mem \
   -out-dir $prefix -out $prefix:t -nb-cores $threads \
   > $prefix/log.out 2> $prefix/log.err 

 rm -f *glue* dummy* *h5 *contigs.fa
 fa=$prefix/$prefix:t.unitigs.fa
 $samtools faidx $fa && $n50 $fa.fai > $fa.n50
 pigz -p $threads $fa

 cd $wd
 
 kk=$k
done
