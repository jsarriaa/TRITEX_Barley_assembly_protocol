#!/bin/zsh

# Write super-scaffolds or pseudomolecules from an assembly FASTA file and an AGP BED file

t=(-tabix '/opt/Bio/bcftools/1.9/bin/tabix')
b=(-bgzip '/opt/Bio/bcftools/1.9/bin/bgzip')
e=(-pathbedtools '/opt/Bio/bedtools/2.26.0/bin/bedtools')
p=(-parallel '/opt/Bio/parallel/20150222/bin/parallel')
o=(-cores 1)

zparseopts -D -K -- -tabix:=t -bgzip:=b -pathbedtools:=e -parallel:=p \
 -bed:=d -cores:=o -assembly:=i

bed=$d[2]
psmolinput=$i[2]
ncores=$o[2]
parallel=$p[2]
bedtools=$e[2]
bgzip=$b[2]
tabix=$t[2]

mktemp | read tmp
bedgz=$tmp.gz

$bgzip -c $bed > $bedgz
$tabix -C -s 4 -b 5 $bedgz

gzip -cd $bedgz | cut -f 4 | uniq \
 | $parallel -k --will-cite -j $ncores \
    $tabix $bedgz '{}' \
    \| $bedtools getfasta -s -fi $psmolinput -bed /dev/stdin -fo /dev/stdout \
    \| "sed '/^>/d'" \| tr -d "'\n'" \
    \| "awk 'BEGIN {print \">{}\"} NF'" \
 | fold -w 60   

rm -f $tmp $bedgz
