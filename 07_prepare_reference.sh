# Obtain from a reference genome (in this case Morex), a single copy 100bp regions

#!/bin/bash

debug_verbose=true


projectdir='/scratch/GDB136/IPK'
projectname='GDB136'
bitbucket="$projectdir/tritexassembly.bitbucket.io"
zsh=$(which zsh)

mask="$bitbucket/miscellaneous/mask_assembly.zsh" 
fa="$projectdir/MorexV3.fa"  # Reference genome FASTA file
fa_base=$(basename "$fa" .fa)
out="${fa_base}_singlecopy_100bp"
bed="$out/${fa_base}_masked_noGaps.bed"

if [ "$debug_verbose" = true ]; then
  echo "Debug mode is ON"
  echo "Project directory: $projectdir"
  echo "Project name: $projectname"
  echo "Bitbucket path: $bitbucket"
  echo "Reference FASTA file path: $fa"
  echo "BED file path: $bed"
  echo "Output file prefix: $out"
else
  echo "Debug mode is OFF"
fi

if [ "$debug_verbose" = true ]; then
  echo "Running the command: zsh $mask --mem 300G --fasta $fa --mincount 2 --out ."
fi

# Remember to edit at the mask script the path to:
# BBDuk
# fatotwobit
# twobitinfo
# kmercounterexact (from BBmap)
# samtools
# bedtools

echo "###################"
"$zsh" "$mask" --mem 300G --fasta "$fa" --mincount 2 --out "$out"
echo "###################"

# Check if the out was properly created
if [ ! -f "${out}.bed" ]; then
  echo "ERROR: BED file was not created: ${out}.bed"
  exit 1
fi

if [ "$debug_verbose" = true ]; then
  echo "Running the command: awk '\$3 - \$2 >= 100' $bed | grep -v chrUn | awk '{print \$0\"\\tseq_\"NR}' | tee $out.bed | bedtools getfasta -fi $fa -bed /dev/stdin -name -fo $out.fasta"
fi

echo "###################"
awk '$3 - $2 >= 100' "$bed" | grep -v chrUn | awk '{print $0"\tseq_"NR}' \
 | tee "$out.bed" | bedtools getfasta -fi "$fa" -bed /dev/stdin -name -fo "$out.fasta"
echo "###################"

# Check if the FASTA file was created successfully
if [ ! -s "$out.fasta" ]; then
  echo "ERROR: Failed to create FASTA file from BED"
  exit 1
else
  echo "FASTA file created successfully: $out.fasta"
fi  

