# Convert contig GFA files to FASTA format
# Change contig names

debug_verbose=true

projectdir='/scratch/GDB136/IPK'
projectname='GDB136'

gfa="$projectdir/assembly/$projectname/$projectname.bp.p_ctg.gfa"

if [ "$debug_verbose" = true ]; then
  echo "Debug mode is ON"
  echo "Project directory: $projectdir"
  echo "Project name: $projectname"
  echo "GFA file path: $gfa"
else
  echo "Debug mode is OFF"
fi

echo "Converting GFA file to FASTA format: $gfa"
if [ ! -f "$gfa" ]; then
  echo "ERROR: GFA file not found: $gfa"
  return 1
fi

gfatools gfa2fa "$gfa" | seqtk rename - contig_ > "$projectdir/assembly/$projectname/$projectname.p_ctg.fa"
if [ "$debug_verbose" = true ]; then
  echo "Command executed: gfatools gfa2fa $gfa | seqtk rename - contig_ > $projectdir/assembly/$projectname/$projectname.p_ctg.fa"
fi

# Check if the FASTA file was created successfully
if [ ! -s "$projectdir/assembly/$projectname/$projectname.p_ctg.fa" ]; then
  echo "ERROR: Failed to create FASTA file from GFA"
  return 1
fi

# Check if contig names were changed successfully
if ! grep -q '^>' "$projectdir/assembly/$projectname/$projectname.p_ctg.fa"; then
  echo "ERROR: Contig names were not changed successfully in the FASTA file."
  return 1
fi

echo "runned succesfully the command: gfatools gfa2fa $gfa | seqtk rename - contig_ > $projectdir/assembly/$projectname/$projectname.p_ctg.fa"