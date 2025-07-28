#!/bin/zsh

# Perform an in silico digest of an assembly FASTA


#######################
### OLD ORIGINAL PATHS

# e=(-restrict '/opt/Bio/EMBOSS/6.6.0/bin/restrict')
# t=(-bedtools '/opt/Bio/bedtools/2.26.0/bin/bedtools')
# b=(-rebase '/filer-dg/agruppen/DG/mascher/rebase_181023/link_emboss_e.txt')
#######################

### JOAN:
# To download the rebase file:
# wget https://rebase.neb.com/rebase/link_emboss_e
# mv link_emboss_e link_emboss_e.txt

e=(-restrict '/scratch/software-phgv2/miniconda3/envs/IPK-tritex/bin/restrict')
t=(-bedtools '/scratch/software-phgv2/miniconda3/envs/IPK-tritex/bin/bedtools')
b=(-rebase '/scratch/GDB136/IPK/link_emboss_e.txt')

zparseopts -D -K -- -bedtools:=t -restrict:=e -rebase:=b -name:=n -sitelen:=l -enzyme:=s -ref:=r -minlen:=m 

bedtools=$t[2]
ref=$r[2]
name=$n[2]
emboss=$e[2]
enzyme=$s[2]
minlen=$m[2]
len=$l[2]
rebase=$b[2]

while [[ $# -gt 0 ]]; do
  case $1 in
    -ref)
      ref="$2"
      shift 2
      ;;
    -enzyme)
      enzyme="$2"
      shift 2
      ;;
    -sitelen)
      len="$2"
      shift 2
      ;;
    -minlen)
      minlen="$2"
      shift 2
      ;;
    *)
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
done

echo "Parsed arguments:"
echo "bedtools: $bedtools"
echo "ref: $ref"
echo "name: $name"
echo "emboss: $emboss"
echo "enzyme: $enzyme"
echo "minlen: $minlen"
echo "len: $len"
echo "rebase: $rebase"

bed=${ref:r}_${enzyme}_fragments_${minlen}bp.bed

echo "output BED file: $bed"

{
 $emboss -enzymes $enzyme -sitelen $len -sequence $ref -datafile $rebase -outfile /dev/stdout | awk NF \
  | awk '/^# Sequence:/ {printf "\n"$3} /^ *[0-9]/ {printf "\t"0+$1}' | awk NF \
  | awk -v len=$len '{for(i = 2; i < NF; i++) print $1"\t"$i+len-1"\t"$(i+1)-1}' \
  | sort -S10G -k 1,1 -k 2,2n \
  | awk '$3 - $2 >= '$minlen > $bed && \
 awk 'BEGIN{OFS="\t"}
      $3 - $2 <= 400 {print $0, $1":"$2"-"$3; next} 
      {print $1, $2, $2+200, $1":"$2"-"$3;
       print $1, $3-200, $3, $1":"$2"-"$3}' $bed \
  | $bedtools nuc -bed - -fi $ref > ${bed:r}_split.nuc.txt 
} 2> ${ref:r}_${enzyme}_fragments_${minlen}bp_digest.err  
