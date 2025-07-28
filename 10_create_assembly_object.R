rscript_path <- "/scratch/software-phgv2/miniconda3/envs/IPK-tritex/bin/R"

# Load the contig assembly, guide map and Hi-C alignment records into R
# to create the assembly object

# All subsequent steps should be carried out in the working directory $projectdir/pseudomolecules

# set a debug_verbose variable
debug_verbose <- TRUE


projectname <- "GDB136"
maindir <- "/scratch/GDB136/IPK/"
pseudodir <- file.path(maindir, "pseudomolecules")
bitbucketdir <- file.path(maindir, "tritexassembly.bitbucket.io")
species <- "barley"
cores <- 8 #cores to use for parallel processing

if (debug_verbose) {
  print(paste("Project name:", projectname))
  print(paste("Main directory:", maindir))
  print(paste("Pseudomolecules directory:", pseudodir))
  print(paste("Bitbucket directory:", bitbucketdir))
}

# move to the pseudomolecules directory
setwd(pseudodir)

if (debug_verbose) {
  print(paste("Current working directory:", getwd()))
}

# Load the R functions for pseudomolecule construction
source(file.path(bitbucketdir, "R", "pseudomolecule_construction.R"))

# Check if the R functions were loaded successfully
if (!exists("read_morexaln_minimap")) {
    stop("Error: The function 'read_morexaln_minimap' is not defined. Please check the R script path or the function definition.")
}

if (debug_verbose) {
  print(paste("R functions for pseudomolecule construction loaded from", file.path(bitbucketdir, "R", "pseudomolecule_construction.R")))
  print("\n")
}

# import the guide map
readRDS(file.path(maindir, paste0(projectname,"_pseudopopseq.Rds"))) -> popseq


if (debug_verbose) {
      print(names(popseq))
      print(head(popseq))
}

# Read contig lengths
f <- file.path(maindir, "assembly", projectname, paste0(projectname, ".p_ctg.fa.fai"))
fread(f, head=F, select=1:2, col.names=c("scaffold", "length")) -> fai

if (debug_verbose) {
  print(paste("Contig lengths read from", f))
  print(head(fai))
}

# Read guide map alignment.
f <- file.path(maindir, paste0(projectname, ".p_ctg_", projectname, ".paf.gz"))
read_morexaln_minimap(paf=f, popseq=popseq, minqual=30, minlen=500, prefix=FALSE) -> morexaln

# Ensure chromosome names match for downstream joins
morexaln[, popseq_alphachr := sub("^chr", "", popseq_alphachr)]
morexaln[, sorted_alphachr := sub("^chr", "", sorted_alphachr)]

if (debug_verbose) {
  print(paste("Guide map alignment read from", f))
  print(head(morexaln))
}

# Read HiC links
hicdir <- file.path(maindir, "HiC")
fread(paste('find', hicdir, '| grep "_pairs.tsv.gz$" | xargs zcat'),
       header=F, col.names=c("scaffold1", "pos1", "scaffold2", "pos2")) -> fpairs

if (debug_verbose) {
  print(paste("Hi-C links read from", hicdir))
  print(head(fpairs))
}

# Init and save assembly, replace species with the right one. Here we are going to use "barley"

print(paste("Initializing assembly for species:", species))

init_assembly(fai=fai, cssaln=morexaln, fpairs=fpairs) -> assembly
anchor_scaffolds(assembly=assembly, popseq=popseq, species=species) -> assembly 
add_hic_cov(assembly, binsize=1e4, binsize2=1e6, minNbin=50, innerDist=3e5, cores=cores) -> assembly   
saveRDS(assembly, file=paste0("assembly_", projectname, ".Rds")) 

if (debug_verbose) {
  print(paste("Assembly object created and saved as assembly_", projectname, ".Rds", sep=""))
  print(paste("Remember: you're working at the directory:", getwd()))
}