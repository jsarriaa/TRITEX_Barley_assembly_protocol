# This is to map again after doing changes after Map Inspector Shiny
# Run 17_MapInspectorShiny.R 


rscript_path <- "/scratch/software-phgv2/miniconda3/envs/IPK-tritex/bin/R"

debug_verbose <- TRUE

projectname <- "GDB136"
maindir <- "/scratch/GDB136/IPK/"
pseudodir <- file.path(maindir, "pseudomolecules")
species <- "barley"
cores <- 8
enzyme <- "DpnII"
assembly_v2_file <- file.path(pseudodir, paste0("assembly_v2_", projectname, ".Rds"))
popseq_file <- file.path(maindir, paste0(projectname, "_pseudopopseq.Rds"))
hic_map_v2_file <- file.path(pseudodir, paste0("hic_map_v2_500kb.Rds"))
assembly_v2 <- readRDS(assembly_v2_file)
snuc <- file.path(maindir, "assembly", projectname, paste0(projectname, ".p_ctg_", enzyme, "_fragments_30bp_split.nuc.txt"))

setwd(pseudodir)

# Load the R functions for pseudomolecule construction
bitbucketdir <- file.path(maindir, "tritexassembly.bitbucket.io")
source(file.path(bitbucketdir, "R", "pseudomolecule_construction.R"))

read_hic_map(rds="hic_map_v2_500kb.Rds", file="hic_map_v2_500kb.xlsx") -> nmap
diff_hic_map(rds="hic_map_v2_500kb.Rds", nmap, species=species)
hic_map(species=species, agp_only=T, map=nmap) -> hic_map_v3

saveRDS(hic_map_v3, file="hic_map_v3.Rds")

# Final Hi-C plots. You can check if the inverted scaffolds are correct.
snuc <- '$projectdir/assembly/project_name.p_ctg_MboI_fragments_30bp_split.nuc.txt'
hic_plots(rds="hic_map_v3.Rds", assembly=assembly_v2,
	cores=1, species=species, nuc=snuc) -> hic_map_v3

