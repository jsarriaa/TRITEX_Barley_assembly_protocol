#!/bin/zsh

# Run the SOAP scaffolding for different link count thresholds for a single library and collate summary statistics.
# GNU Parallel is used for parallelization.

t=(-threads 1)
v=(-interval 3:25)
p=(-outdir '.')
g=(-pegrads "")
s=(-soap '/filer-dg/agruppen/seq_shared/mascher/source/SOAPdenovo2-bin-LINUX-generic-r240/SOAPdenovo-127mer')
nn=(-n50 '/filer-dg/agruppen/seq_shared/mascher/code_repositories/triticeae.bitbucket.io/shell/n50')
pp=(-parallel '/opt/Bio/parallel/20150222/bin/parallel')
sam=(-samtools '/opt/Bio/samtools/1.7/bin/samtools')

zparseopts -D -E -K -- -outdir:=p -input:=i -name:=n -size:=l -interval:=v -pegrads:=g -threads:=t \
 -soap:=s -n50:=nn -parallel:=pp -samtools:=sam

samtools=$sam[2]
parallel=$pp[2]
n50=$nn[2]
soap=$s[2]
outdir=$p[2]
input=$i[2]
lib=$n[2]
size=$l[2]
interval=$v[2]
grads=$g[2]
threads=$t[2]

tr ':' ' ' <<< $interval | xargs seq | while read j; do
 dir=$outdir/$input:t'_'$lib'_l'$j
 prefix=$dir/$input:t

 if ! mkdir $dir; then
  echo "Directory $dir cannot be created." 1>&2
  exit 1
 fi

 ln -s $input/*preGraphBasic $input/*Arc $input/*readOnContig.gz $input/*readInGap.gz $input/*updated.edge \
  $input/*ContigIndex $input/*contig $dir
 if [[ -z $grads ]]; then
  cp $input/*peGrads $dir
 else
  cp $grads $dir
 fi
 sed -i '/^'$size'/s/[0-9]\+$/'$j'/' $dir/*peGrads
 echo $prefix
done | $parallel --will-cite -j $threads \
 $soap scaff -g '{}' \> '{}'_scaf.out 2\> '{}'_scaf.err 

base=$input:t'_'$lib'_l'
find $outdir | grep 'scafSeq$' | grep $base | $parallel -j $threads --will-cite $samtools faidx '{}' && \
find $outdir | grep 'scafSeq$' | grep $base | $parallel -j $threads --will-cite $n50 '{}'.fai \> '{}'.n50 && \
find $outdir | grep 'scafSeq$' | grep $base | xargs rm -f && \
find $outdir -maxdepth 1 -type d | grep $base | while read i; do
 rev <<< $i | cut -d _ -f 1 | rev | cut -d l -f 2 | read len
 find $i | grep 'n50$' | read n
 find $i | grep 'scafStatistics$' | read a
 grep '^N50' $n | cut -f 2 | read m
 grep '^N90' $n | cut -f 2 | read m9
 grep '^len1M' $n | cut -f 2 | read s
 grep '^csize' $n | cut -f 2 | read t
 grep -m 1 '^Size_withoutN' $a | cut -f 2 | read w
 echo $len $t $s $m $m9 $w | tr ' ' '\t'
done | sort -nk 1
