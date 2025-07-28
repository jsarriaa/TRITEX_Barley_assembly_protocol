# compile pseudomolecules
# fasta in output folder


rscript_path <- "/scratch/software-phgv2/miniconda3/envs/IPK-tritex/bin/R"

debug_verbose <- TRUE

projectname <- "GDB136"
maindir <- "/scratch/GDB136/IPK/"
outdir <- file.path(maindir, "output")
pseudodir <- file.path(maindir, "pseudomolecules")
species <- "barley"
cores <- 8
enzyme <- "DpnII"
fasta <- file.path(maindir, "assembly", projectname, paste0(projectname, ".p_ctg.fa"))
assembly_v2_file <- file.path(pseudodir, paste0("assembly_v2_", projectname, ".Rds"))
hic_map_v2_file <- file.path(pseudodir, paste0("hic_map_v2_500kb.Rds"))

assembly_v2 <- readRDS(assembly_v2_file)
hic_map_v3 <- readRDS(hic_map_v2_file)
# Doing it with the hic_map_v2 object, which contains the scaffolds >= 500 kb
# Instead of 3 as the protocol, due i have not check scaffolds and made the v3

# Load the R functions for pseudomolecule construction
bitbucketdir <- file.path(maindir, "tritexassembly.bitbucket.io")
source(file.path(bitbucketdir, "R", "pseudomolecule_construction.R"))
library(openxlsx)


# work in outdir
setwd(outdir)


print("Remember edit in pseudomolecule_construction.R the path to bedtools, samtools, tabix and bgzip if you have not done it yet.  ")
print("TAKE A LOT OF CARE ABOUT THESSAMTOOLS PATHS, THEY ARE NOT THE SAME AS IN THE PREVIOUS STEPS")

sink("pseudomolecules_v1.log") 
compile_psmol(fasta=fasta, output="pseudomolecules_v1",
	hic_map=hic_map_v3, assembly=assembly_v2, cores=cores)
sink() 