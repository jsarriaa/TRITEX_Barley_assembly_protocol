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

packages <- c("shiny", "data.table", "DT")
installed <- packages %in% rownames(installed.packages())
if (any(!installed)) {
	options(repos = c(CRAN = "https://cloud.r-project.org"))  # Set CRAN mirror	
	#install.packages("shiny", dependencies = TRUE, repos = "https://cloud.r-project.org")
	#install.packages("DT", dependencies = TRUE, repos = "https://cloud.r-project.org")
	#install.packages("httpuv", dependencies = TRUE, repos = "https://cloud.r-project.org")}
}
print("Running Map Inspector Shiny...")
library(rmarkdown)
library(shiny)
run("/scratch/GDB136/IPK/tritexassembly.bitbucket.io/R/map_inspector.Rmd")
print("Map Inspector Shiny completed.")



