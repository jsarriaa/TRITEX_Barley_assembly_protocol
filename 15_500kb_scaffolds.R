# Decrease size to 500 kb # and plot Hi-C map

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
fbed <- file.path(maindir, "assembly", projectname, paste0(projectname, ".p_ctg_", enzyme, "_fragments_30bp.bed"))
snuc <- file.path(maindir, "assembly", projectname, paste0(projectname, ".p_ctg_", enzyme, "_fragments_30bp_split.nuc.txt"))

setwd(pseudodir)

# Load the R functions for pseudomolecule construction
bitbucketdir <- file.path(maindir, "tritexassembly.bitbucket.io")
source(file.path(bitbucketdir, "R", "pseudomolecule_construction.R"))

read_fragdata(info=assembly_v2$info, file=fbed) -> frag_data
frag_data$info[!is.na(hic_chr) & length >= 5e5, .(scaffold, nfrag, length, chr=hic_chr, cM=popseq_cM)] -> hic_info
# Note that we are now filtering for contigs >= 500 kb

hic_map(info=hic_info, assembly=assembly_v2, frags=frag_data$bed,
	species=species, ncores=cores, min_nfrag_scaffold=30,
	max_cM_dist = 1000, binsize=2e5, min_nfrag_bin=10, gap_size=100) -> hic_map_v2

saveRDS(hic_map_v2, file="hic_map_v2_500kb.Rds")

# Make the Hi-C plots, files will be placed in your working directory.
hic_plots(rds="hic_map_v2_500kb.Rds", assembly=assembly_v2,
	cores=cores, species=species, nuc=snuc) -> hic_map_v2

# Exporting to Excel the table that will be used for manual curation.
write_hic_map(rds="hic_map_v2_500kb.Rds", file="hic_map_v2_500kb.xlsx", species=species)

print("Hi-C map v2 with contigs >= 500 kb saved successfully.")
print("excel file: hic_map_v2_500kb.xlsx")
