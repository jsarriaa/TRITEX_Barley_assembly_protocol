# Hi-C mapping against the in silico digested assembly (from Hi-Fi reads)
#!/bin/bash


debug_verbose=true
projectdir='/scratch/GDB136/IPK'
projectname='GDB136'
bitbucket="$projectdir/tritexassembly.bitbucket.io"
ref="$projectdir/assembly/$projectname/$projectname.p_ctg.fa"
hic_directory="$projectdir/HiC"
minlen=30  # Minimum read length after adapter trimming
maxlen=500  # Maximum read length to consider for analysis

if [ ${debug_verbose} = true ]; then
  echo "Note that you are using a version of TRITEX without private software, so you are avoiding novosort and using an open-source alternative: samtools."
fi

# create the mapping directory if it does not exist
if [ ! -d "$projectdir/mapping" ]; then
  echo "Creating mapping directory: $projectdir/mapping"
  mkdir -p "$projectdir/mapping"
fi
cd $projectdir/mapping/

enzyme='DpnII'     # Adjust to your own enzyme
linker='GATCGATC'  # This is the linker for DpnII enzyme, adjust to yours if necessary

map="$projectdir/run_hic_mapping_samtools_joan.sh"
bed="$projectdir/assembly/$projectname/${projectname}.p_ctg_${enzyme}_fragments_30bp.bed"

# Generate the temp dir if it does not exist
if [ ! -d '/scratch/GDB136/tmp' ]; then
  echo "Creating temporary directory: /scratch/GDB136/tmp"
  mkdir -p '/scratch/GDB136/tmp'
fi
TMPDIR='/scratch/GDB136/tmp'  # Temporary directory for intermediate files

if [ "$debug_verbose" = true ]; then  echo "Debug mode is ON"
  echo "Project directory: $projectdir"
  echo "Project name: $projectname"
  echo "Reference FASTA file path: $ref"
  echo "Bitbucket path: $bitbucket"
  echo "Mapping script path: $map"
  echo "Linker sequence: $linker"
  echo "Restriction enzyme: $enzyme"
  echo "BED file path: $bed"
  echo "Temporary directory: $TMPDIR"
  echo "Hi-C directory: $hic_directory"
  echo "minlen: $minlen"
  echo "maxlen: $maxlen"

echo "Running the command: bash $map $hic_directory -threads 64 -mem '200G' -linker $linker -ref $ref -bed $bed -tmp $TMPDIR -minlen $minlen -maxlen $maxlen hic"    
echo "remember to adjust at the $map script the path to: cutadpt, bgzip, minimap2, novosort, samtools and bedtools"

else
  echo "Debug mode is OFF"
fi

echo "###################"

bash "$map" $hic_directory -threads 64 -mem '4G' -linker "$linker" -ref "$ref" -bed "$bed" -tmp "$TMPDIR" -minlen "$minlen" -maxlen "$maxlen" hic
# JOAN NOTE: 
# mem 200G is exploting, too much memory. Tryinng with 4G

cd $projectdir

echo "###################"


echo "Hi-C mapping completed. Check the mapping directory: $projectdir/mapping/"