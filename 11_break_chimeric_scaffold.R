# Generate diagnostic plots to check for chimeras in the input assembly. If sequences from different chromosomes are joined, the chimeric scaffolds will have guide map markers from more than one chromosome.
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

if (debug_verbose) {
  print(paste("Assembly object loaded from", assembly_file))
#  print(str(assembly))
}

# Load the R functions for pseudomolecule construction
bitbucketdir <- file.path(maindir, "tritexassembly.bitbucket.io")
source(file.path(bitbucketdir, "R", "pseudomolecule_construction.R"))

# Create diagnostic plots for contigs longer than 1 Mb.

if (debug_verbose) {
print("Creating diagnostic plots for chimeric scaffolds...")
print("Command: plot_chimeras(assembly=assembly, scaffolds=s, species=species, refname=projectname, autobreaks=F, mbscale=1, file=\"assembly_1Mb.pdf\", cores=1)")
}

assembly$info[length >= 1e6, .(scaffold, length)][order(-length)] -> s
plot_chimeras(assembly=assembly, scaffolds=s, species=species,
    refname=projectname, autobreaks=F, mbscale=1, file="assembly_1Mb.pdf", cores=1) 

if (debug_verbose) {
  print("Diagnostic plots for chimeric scaffolds created.")
  print("Check the file 'assembly_1Mb.pdf' for details.")
}

print("Checking chromosome assignments in assembly$info:")
print(head(assembly$info[, .(scaffold, popseq_chr, sorted_chr, hic_chr)]))

print("Chromosome mapping table:")
chrNames(species=species)->wheatchr
print(wheatchr)
