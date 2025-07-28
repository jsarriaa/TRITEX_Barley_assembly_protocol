#!/bin/zsh

# Run the final SOAP scaffolding (including gap-filling) after optimal threshold have been determined.

t=(-threads 1)
g=(-gapdiff 50) #Default for SOAPDenovo2 scaff -F
p=(-pegrads "")
d=(-outdir '.')
sp=(-soap '/filer-dg/agruppen/seq_shared/mascher/source/SOAPdenovo2-bin-LINUX-generic-r240/SOAPdenovo-127mer')
s=(-suffix 'SCAFF')

zparseopts -D -E -K -- -input:=i -threads:=t -outdir:=d -suffix:=s -gapdiff:=g -pegrads:=p -soap:=sp

soap=$sp[2]
outdir=$d[2]
threads=$t[2]
input=$i[2]
grads=$p[2]
suffix=$s[2]
gap=$g[2]

if [[ -z $suffix ]]; then
 echo "Please specify suffix." 1>&2
 exit 1
fi

dir=$outdir/$input:t'_'$suffix
prefix=$dir/$input:t

if ! mkdir $dir; then
 echo "Directory $dir cannot be created." 1>&2
 exit 1
fi

{
 ln -s $input/*preGraphBasic $input/*Arc $input/*readOnContig.gz $input/*readInGap.gz $input/*updated.edge \
  $input/*ContigIndex $input/*contig $dir

 if [[ -z $grads ]]; then
  cp $input/*peGrads $dir
 else
  cp $grads $dir
 fi

 $soap scaff -F -G $gap -g $prefix 
} > ${prefix}_scaf.out 2> ${prefix}_scaf.err 
