# Check the stats of contigs and untigs of hifiasm assembly

projectdir='/scratch/GDB136/IPK'
projectname='GDB136'
bitbucket="$projectdir/tritexassembly.bitbucket.io"

debug_verbose=True

cd "$projectdir/assembly/$projectname" || { echo "ERROR: no existe $projectdir/assembly/$projectname"; return 1; }

echo "Checking contig and untig stats for $projectname"

echo "starting with contigs running n50 script from bitbucket"

grep '^S' $projectname.bp.p_ctg.noseq.gfa | cut -f 4 | cut -d : -f 1,3 \
 | tr : '\t' | $bitbucket/shell/n50 > $projectname.p_ctg.noseq.gfa.n50

# Check if the file is properly created and has lines
if [ ! -s $projectname.p_ctg.noseq.gfa.n50 ]; then
    echo "ERROR: Failed to create n50 file for contigs"
    echo failed running the command: grep '^S' $projectname.bp.p_ctg.noseq.gfa | cut -f 4 | cut -d : -f 1,3 | tr : '\t' | $bitbucket/shell/n50 > $projectname.p_ctg.noseq.gfa.n50
    if [ "$debug_verbose" = True ]; then
        echo "Debug: Check if the input GFA files exist and are not empty."
        ls -l $projectname.bp.p_ctg.noseq.gfa $projectname.bp.p_utg.noseq.gfa
    fi
    cd $projectdir || { echo "ERROR: no existe $projectdir"; return 1; }
    return 1
fi


echo "starting with untigs running n50 script from bitbucket"

if ! grep '^S' $projectname.bp.p_utg.noseq.gfa | cut -f 4 | cut -d : -f 1,3 | tr : '\t' | grep -q .; then
    echo "ERROR: No valid sequence lines found in $projectname.bp.p_utg.noseq.gfa"
    return 1
fi

grep '^S' $projectname.bp.p_utg.noseq.gfa | cut -f 4 | cut -d : -f 1,3 \
 | tr : '\t' | $bitbucket/shell/n50 > $projectname.p_utg.noseq.gfa.n50

# Check if the file is properly created and has lines
if [ ! -s $projectname.p_utg.noseq.gfa.n50 ]; then
    echo "ERROR: Failed to create n50 file for untigs"
    echo failed running the command: grep '^S' $projectname.bp.p_utg.noseq.gfa | cut -f 4 | cut -d : -f 1,3 | tr : '\t' | $bitbucket/shell/n50 > $projectname.p_utg.noseq.gfa.n50
    if [ "$debug_verbose" = True ]; then
        echo "Debug: Check if the input GFA files exist and are not empty."
        ls -l $projectname.bp.p_ctg.noseq.gfa $projectname.bp.p_utg.noseq.gfa
    fi
    cd $projectdir || { echo "ERROR: no existe $projectdir"; return 1; }
    return 1
fi

echo "outputs are $projectdir/assembly/$projectname/${projectname}.p_ctg.noseq.gfa.n50 and ${projectdir}/assembly/$projectname/${projectname}.p_utg.noseq.gfa.n50"
cd $projectdir || { echo "ERROR: no existe $projectdir"; return 1; }