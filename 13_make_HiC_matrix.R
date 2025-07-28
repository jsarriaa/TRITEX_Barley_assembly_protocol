# Now we make the Hi-C contact matrix, or Hi-C map. Three files will be generated: 
# a PDF with the intrachromosomal plots, a PDF with the interchromosomal plots, and 
# a file to be loaded into R Shiny Map inspector (*_export.Rds).

rscript_path <- "/scratch/software-phgv2/miniconda3/envs/IPK-tritex/bin/R"

debug_verbose <- TRUE

projectname <- "GDB136"
maindir <- "/scratch/GDB136/IPK/"
pseudodir <- file.path(maindir, "pseudomolecules")
species <- "barley"
cores <- 8
enzyme <- "DpnII"
assembly_v2_file <- file.path(pseudodir, paste0("assembly_v2_", projectname, ".Rds"))
snuc <- file.path(maindir, "assembly", projectname, paste0(projectname, ".p_ctg_", enzyme, "_fragments_30bp_split.nuc.txt"))

if (!file.exists(assembly_v2_file)) {
    stop("assembly_v2 file does not exist. Please run 12_fix_chimeras_scaffold.R first.")
}

assembly_v2 <- readRDS(assembly_v2_file)

setwd(pseudodir)

# Load the R functions for pseudomolecule construction
bitbucketdir <- file.path(maindir, "tritexassembly.bitbucket.io")
source(file.path(bitbucketdir, "R", "pseudomolecule_construction.R"))

fbed <- file.path(maindir, "assembly", projectname, paste0(projectname, ".p_ctg_", enzyme, "_fragments_30bp.bed"))
read_fragdata(info=assembly_v2$info, file=fbed) -> frag_data

# Debug chr data to avoid errors
if (debug_verbose) {
	print("Checking chromosome assignments in assembly_v2$info:")
	print(table(is.na(frag_data$info$hic_chr)))
    print(head(frag_data$info))
}

if (all(is.na(frag_data$info$hic_chr))) {
    stop("All hic_chr values are NA. No contigs were assigned to chromosomes!")
}

# Consider only contigs >= 1 Mb first
frag_data$info[!is.na(hic_chr) & length >= 1e6, 
    .(scaffold, nfrag, length, chr=paste0("chr", hic_chr), cM=popseq_cM)] -> hic_info
	
hic_map(info=hic_info, assembly=assembly_v2, frags=frag_data$bed,
	species=species, ncores=cores, min_nfrag_scaffold=30,
	max_cM_dist = 1000, binsize=2e5, min_nfrag_bin=10, gap_size=100) -> hic_map_v1
saveRDS(hic_map_v1, file="hic_map_v1.Rds")

if (debug_verbose) {
	print("Hi-C map v1 saved successfully.")
}

# Make the Hi-C plots, files will be placed in your working directory.
projectdir <- maindir
snuc <- file.path(projectdir, "assembly", projectname, paste0(projectname, ".p_ctg_", enzyme, "_fragments_30bp_split.nuc.txt"))

hic_plots(rds="hic_map_v1.Rds", assembly=assembly_v2,
	cores=8, species=species, nuc=snuc) -> hic_map_v2

saveRDS(hic_map_v2, file="hic_map_v2.Rds") 
print("Hi-C map v2 saved successfully.")

hic_plots("hic_map_v2.Rds", assembly=assembly_v2, cores=cores,
	species=species, nuc=snuc) -> hic_map_v2
print("Hi-C plots v2 generated successfully.")

# Show the map stats
if (debug_verbose) {
	print("Hi-C map v2 generated successfully.")
	print(hic_map_v2$stats)
}

# Exporting to Excel the table that will be used for manual curation.
write_hic_map(rds="hic_map_v2.Rds", file="hic_map_v2.xlsx", species=species)
print("Hi-C map v2 exported to Excel successfully.")

