#!/bin/zsh

# Check the concordance between the input assembly and the super-scaffolds/pseudomolecules
# The `diff` output should be empty.

t=(-tabix '/opt/Bio/bcftools/1.9/bin/tabix')
b=(-bgzip '/opt/Bio/bcftools/1.9/bin/bgzip')
e=(-bedtools '/opt/Bio/bedtools/2.26.0/bin/bedtools')
s=(-samtools '/opt/Bio/samtools/1.9/bin/samtools')

zparseopts -D -K -- -tabix:=t -bgzip:=b -bedtools:=e -samtools:=s \
 -assembly:=a -agp:=g -psmol:=p

samtools=$s[2]
bedtools=$e[2]
bgzip=$b[2]
tabix=$t[2]
agpbed=$g[2]
assembly=$a[2]
psmol=$p[2]

mktemp | read tmp
sort -k 4,4 -k 5,5n $agpbed | $bgzip > $tmp.gz
$tabix -s 4 -b 5 $tmp.gz

cut -f 4 $agpbed | grep -v gap | grep -Fwf - \
 <(cut -f 1 $assembly.fai) | xargs $tabix $tmp.gz \
 | $bedtools getfasta -fi $psmol -name -s -bed /dev/stdin -fo /dev/stdout \
 | sed -e 's/(-)$//' -e 's/(+)$//' \
 | sed 's/::.*$//' | diff /dev/stdin \
    <(cut -f 4 $agpbed | grep -v gap | grep -Fwf - <(cut -f 1 $assembly.fai) | xargs $samtools faidx $assembly \
       | awk '/^>/ {print "\n"$1; next} {printf $0}' | awk NF) 

rm -f $tmp $tmp.gz
