# Get some stats about the pseudomolecules

rscript_path <- "/scratch/software-phgv2/miniconda3/envs/IPK-tritex/bin/R"

debug_verbose <- TRUE

projectname <- "GDB136"
maindir <- "/scratch/GDB136/IPK/"
pseudodir <- file.path(maindir, "pseudomolecules")
species <- "barley"
cores <- 8
enzyme <- "DpnII"
hic_map_v2_file <- file.path(pseudodir, paste0("hic_map_v2_500kb.Rds"))

hic_map_v3 <- readRDS(hic_map_v2_file)
# Doing it with the hic_map_v2 object, which contains the scaffolds >= 500 kb
# Instead of 3 as the protocol, due i have not check scaffolds and made the v3

setwd(pseudodir)

# Load the R functions for pseudomolecule construction
bitbucketdir <- file.path(maindir, "tritexassembly.bitbucket.io")
source(file.path(bitbucketdir, "R", "pseudomolecule_construction.R"))
library(openxlsx)


hic_map_v3$agp[agp_chr != "chrUn" & gap == F][, .("Ncontig" = .N,
	"N50 (Mb)"=round(n50(scaffold_length/1e6), 1),
	"max_length (Mb)"=round(max(scaffold_length/1e6),1),
	"min_length (kb)"=round(min(scaffold_length/1e3),1)), key=agp_chr] -> res
hic_map_v3$chrlen[, .(agp_chr, "length (Mb)"=round(length/1e6, 1))][res, on="agp_chr"] -> res
setnames(res, "agp_chr", "chr")

hic_map_v3$agp[gap == F & agp_chr != "chrUn"][, .("Ncontig" = .N,
	"N50 (Mb)"=round(n50(scaffold_length/1e6), 1),
	"max_length (Mb)"=round(max(scaffold_length/1e6),1),
	"min_length (kb)"=round(min(scaffold_length/1e3),1))][, agp_chr := "1-7"] -> res2  # 7, barley chromosomes
hic_map_v3$chrlen[, .(agp_chr="1-7", "length (Mb)"=round(sum(length)/1e6, 1))][res2, on="agp_chr"] -> res2
setnames(res2, "agp_chr", "chr")

hic_map_v3$agp[gap == F & agp_chr == "chrUn"][, .(chr="un", "Ncontig" = .N,
	"N50 (Mb)"=round(n50(scaffold_length/1e6), 1),
	"max_length (Mb)"=round(max(scaffold_length/1e6),1),
	"min_length (kb)"=round(min(scaffold_length/1e3),1),
	"length (Mb)"=sum(scaffold_length/1e6))] -> res3

rbind(res, res2, res3) -> res

write.xlsx(res, file="hic_map_v3_pseudomolecule_stats.xlsx")

print("Pseudomolecule stats saved as hic_map_v3_pseudomolecule_stats.xlsx")

