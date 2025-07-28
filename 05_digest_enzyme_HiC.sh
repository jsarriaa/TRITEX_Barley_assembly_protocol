# digest in silico the hifiasm assembly with the restriction enzyme used for the preparation of the Hi-C libraries

#!/bin/bash

debug_verbose=true

projectdir='/scratch/GDB136/IPK'
projectname='GDB136'
bitbucket="$projectdir/tritexassembly.bitbucket.io"
ref="$projectdir/assembly/$projectname/$projectname.p_ctg.fa"

enzyme='DpnII'

echo "Digesting in silico the assembly with the restriction enzyme: $enzyme"
if [ ! -f "$ref" ]; then
  echo "ERROR: Reference FASTA file not found: $ref"
  return 1
fi

if [ "$debug_verbose" = true ]; then
  echo "Debug mode is ON"
  echo "Project directory: $projectdir"
  echo "Project name: $projectname"
  echo "Reference FASTA file path: $ref"

echo "In the script "digest_emboss.zsh", change the paths to the executables of restrict, bedtools and rebase."

else
  echo "Debug mode is OFF"
fi

echo "running the script: zsh $bitbucket/shell/digest_emboss.zsh -ref $ref -enzyme \"$enzyme\" -sitelen 4 -minlen 30"
if [ ! -x "$bitbucket/shell/digest_emboss.zsh" ]; then
  echo "ERROR: Script not found or not executable: $bitbucket/shell/digest_emboss.zsh"
  return 1
fi

zsh $bitbucket/shell/digest_emboss.zsh -ref $ref -enzyme "$enzyme" -sitelen 4 -minlen 30
