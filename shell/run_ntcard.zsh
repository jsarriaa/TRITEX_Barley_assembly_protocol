#!/bin/zsh

# Obtain assembly size estimates for different k-mer using ntCard and KmerStream

m=(-maxcov 100)
t=(-threads 1)
a=(-kmers '100:10:400')
n=(-ntcard '/filer-dg/agruppen/seq_shared/mascher/source/ntCard/ntcard')
k=(-estimate 0)

zparseopts -D -E -K -- -ntcard:=n -estimate:=k \
 -kmers:=a -reads:=r  -outdir:=d -maxcov:=m -threads:=t

range=$a[2]
input=$r[2]
dir=$d[2]
maxcov=$m[2]
threads=$t[2]
ntcard=$n[2]

# Note the Python executable in $PATH will be used.
if [[ $k[2] == 0 ]]; then
 export PYTHONPATH="/filer-dg/agruppen/seq_shared/mascher/source/python2_packages"
 estimate='python /filer-dg/agruppen/seq_shared/mascher/source/KmerStream/KmerStreamEstimate.py'
else
 estimate="python ${k[2]}"
fi

export PATH="${ntcard:h}:$PATH"

prefix=$dir/${input:t:r:r}_ntcard

{
 tr ':' '\t' <<< $range | xargs seq | xargs | tr ' ' , | read kmers

 $ntcard -k $kmers -t $threads -c $maxcov -p $prefix $input

 prefix=$prefix'_k'

 {
  print "k\tF0\tf1\tF1\tsize\terror\tcov"
  {
   echo ""
   tr ':' '\t' <<< $range | xargs seq | while read k; do
    grep F0 ${prefix}$k.hist | cut -f 2 | read F0
    grep F1 ${prefix}$k.hist | cut -f 2 | read F1
    grep f1 ${prefix}$k.hist | cut -f 2 | read f1
    echo 1 $k $F0 $f1 $F1
   done 
  } | eval $estimate /dev/stdin | tr ' ' '\t' | cut -f 2,3,4,5,7- | tail -n +2 \
    | awk '{print $1"\t"$2/1e9"\t"$3/1e9"\t"$4/1e9"\t"$5/1e9"\t"$6"\t"$7}'
 } > $dir/ntcard_summary.tsv
} > $prefix.out 2> $prefix.err 
