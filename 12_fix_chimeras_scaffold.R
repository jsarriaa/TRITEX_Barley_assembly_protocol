# Check the plots and search for drops in Hi-C coverage, misalignments and chimeras. 
#For each case, assign the PDF page number from "assembly_1Mb.pdf" and the positions where it needs breaking.

rscript_path <- "/scratch/software-phgv2/miniconda3/envs/IPK-tritex/bin/R"

debug_verbose <- TRUE

projectname <- "GDB136"
maindir <- "/scratch/GDB136/IPK/"
pseudodir <- file.path(maindir, "pseudomolecules")
species <- "barley"
cores <- 8

setwd(pseudodir)

assembly_file <- file.path(pseudodir, paste0("assembly_", projectname, ".Rds"))
assembly <- readRDS(assembly_file)

# Load the R functions for pseudomolecule construction
bitbucketdir <- file.path(maindir, "tritexassembly.bitbucket.io")
source(file.path(bitbucketdir, "R", "pseudomolecule_construction.R"))

## Here you are suposed to check the diagnostic plots for chimeric scaffolds.

# ChecK:
# 8, 39, 45, 192 .... # not checking anything else

assembly$info[length >= 1e6, .(scaffold, length)][order(-length)] -> ss

i=8 # page number in PDF
ss[i]$scaffold -> s
# change the bin position observed in the plots
assembly$cov[s, on='scaffold'][bin >= 20e6 & bin <= 25e6][order(r)][1, .(scaffold, bin)] -> b
# First time b has to be assigned, so it is not a rbind

i=39 # page number in PDF
ss[i]$scaffold -> s
# change the bin position observed in the plots
rbind(b, assembly$cov[s, on='scaffold'][bin >= 10e6 & bin <= 15e6][order(r)][1, .(scaffold, bin)])->b

i=45 # page number in PDF
ss[i]$scaffold -> s
# change the bin position observed in the plots
rbind(b, assembly$cov[s, on='scaffold'][bin >= 7e6 & bin <= 12e6][order(r)][1, .(scaffold, bin)])->b

i=158 # page number in PDF
ss[i]$scaffold -> s
# change the bin position observed in the plots
rbind(b, assembly$cov[s, on='scaffold'][bin >= 3e6 & bin <= 5e6][order(r)][1, .(scaffold, bin)])->b

i=192 # page number in PDF
ss[i]$scaffold -> s
# change the bin position observed in the plots
rbind(b, assembly$cov[s, on='scaffold'][bin >= 1e6 & bin <= 4e6][order(r)][1, .(scaffold, bin)])->b

i=221 # page number in PDF
ss[i]$scaffold -> s
# change the bin position observed in the plots
rbind(b, assembly$cov[s, on='scaffold'][bin >= 2e6 & bin <= 3e6][order(r)][1, .(scaffold, bin)])->b 

i=409 # page number in PDF
ss[i]$scaffold -> s
# change the bin position observed in the plots
rbind(b, assembly$cov[s, on='scaffold'][bin >= 1e6 & bin <= 2e6][order(r)][1, .(scaffold, bin)])->b


# Rename column and plot again to double-check the breakpoints are correct
setnames(b, "bin", "br")
plot_chimeras(assembly=assembly, scaffolds=b, br=b, species=species,
	refname=projectname,  mbscale=1, file="assembly_chimeras_final.pdf", cores=cores)


if (debug_verbose) {
    message("Chimeras have been fixed. Please check the plots in assembly_chimeras_final.pdf")
}

# Check everything is ok, and proceed to the next step

# If the breakpoints are correct, you can proceed to the next step.

break_scaffolds(b, assembly, prefix="contig_corrected_v1_", slop=1e4, cores=cores, species=species) -> assembly_v2
saveRDS(assembly_v2, file=paste0("assembly_v2_", projectname, ".Rds"))
