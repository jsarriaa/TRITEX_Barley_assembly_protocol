# Second round of scaffolding for Tritex assembly, checking for chimeras
# Broken scaffolds in the last step are renamed in the plots. 
# It might not be necessary to break any more scaffolds; in this case, 
# just check the plots and if everything is alright, proceed to decrease the contig size to 500 kb.

rscript_path <- "/scratch/software-phgv2/miniconda3/envs/IPK-tritex/bin/R"

debug_verbose <- TRUE

projectname <- "GDB136"
maindir <- "/scratch/GDB136/IPK/"
pseudodir <- file.path(maindir, "pseudomolecules")
species <- "barley"
cores <- 8
enzyme <- "DpnII"
assembly_v2_file <- file.path(pseudodir, paste0("assembly_v2_", projectname, ".Rds"))

assembly_v2 <- readRDS(assembly_v2_file)

setwd(pseudodir)

# Load the R functions for pseudomolecule construction
bitbucketdir <- file.path(maindir, "tritexassembly.bitbucket.io")
source(file.path(bitbucketdir, "R", "pseudomolecule_construction.R"))


### USE THIS CODE IF:
    # If you'd like to  plot specific contigs:
#data.table(scaffold=c('contig_corrected_v1_8')) -> ss

    # If you'd like to check them all again:
assembly_v2$info[length >= 1e6, .(scaffold, length)][order(-length)] -> ss
####

plot_chimeras(assembly=assembly_v2, scaffolds=ss,
        species=species, refname=projectname, autobreaks=F,
        mbscale=1, file="assembly_v2_chimeras.pdf", cores=cores)

print("Chimera plots saved successfully.")



# In case there are more chimeras, proceed like the previous step. Assign `i` with the PDF page and don't forget to change the bin. 
# Otherwise, plot the Hi-C map with contigs >= 500 kb.

# Don’t forget to change the output files and object names (assembly_vX, hic_map_vX).

# if no more chimeras are found, proceed to decreasing the size to 500 kb. (In next script)