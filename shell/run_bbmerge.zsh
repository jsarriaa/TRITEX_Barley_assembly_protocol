#!/bin/zsh

# Call BBMerge to merge overlapping paired-end hist and create a fragment size histogram (*_ihist.txt)

b=(-bbmerge '/filer-dg/agruppen/seq_shared/mascher/source/BBMap_37.93/bbmerge.sh')
s=(-strictness 'verystrict')
t=(-threads 1)

zparseopts -K -D -- -threads:=t -strictness:=s -mem:=e -bbmerge:=b

bbmerge=$b[2]
threads=$t[2]
strictness=$s[2]

tr ':' '\t' <<< $1 | read r1 r2 out

$bbmerge in=$r1 in2=$r2 \
 out=${out}_bbmerge.fq.gz \
 ihist=${out}_ihist.txt \
 adapter=default \
 ${strictness}=t \
 t=$threads \
 > ${out}_bbmerge.out 2> ${out}_bbmerge.err 
