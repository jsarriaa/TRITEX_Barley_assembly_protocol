# Optional step to plot collinearities between the scaffolds, helping to spot inverted o missplaced stuff

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
popseq <- readRDS(popseq_file)
hic_map_v2 <- readRDS(hic_map_v2_file)

setwd(pseudodir)

# Load the R functions for pseudomolecule construction
bitbucketdir <- file.path(maindir, "tritexassembly.bitbucket.io")
source(file.path(bitbucketdir, "R", "pseudomolecule_construction.R"))

assembly_v2$cssaln[, .(scaffold, pos, guide_chr=paste0("chr", "", popseq_alphachr), guide_pos=popseq_cM)] -> assembly_v2
hic_map_v2$agp[, .(scaffold, agp_start, agp_end, orientation, agp_chr)][assembly_v2, on="scaffold"]-> assembly_v2

assembly_v2[orientation == 1, agp_pos :=  agp_start + pos]
assembly_v2[orientation == -1, agp_pos :=  agp_end - pos]

# Remove projectname_ prefix from agp_chr
assembly_v2[, agp_chr := sub(paste0("^", projectname, "_"), "", agp_chr)]

assembly_v2[agp_chr == guide_chr] -> assembly_v2



print("starting collinear plot...")
pdf("hic_map_v2_vs_assembly_v2.pdf", height=9, width=9)
par(mar=c(5,5,3,3))
invisible(lapply(chrNames(species=species, agp=TRUE)$agp_chr[1:7], function(i) {   # 7 chr in barley
 print(paste("Debugging chromosome:", i))
 assembly_v2[agp_chr == i & !is.na(agp_pos) & !is.na(guide_pos), 
            plot(pch=".", main=i, agp_pos/1e6, guide_pos/1e6, xlab="Hi-C AGP", ylab="guide", las=1, bty='l')] 
            print(paste("Number of points for chromosome", i, ":", nrow(assembly_v2[agp_chr == i & !is.na(agp_pos) & !is.na(guide_pos)])))
 hic_map_v2$agp[agp_chr == i & gap == T][, abline(v=agp_start/1e6, col="gray")] 
}))
dev.off()

print("Collinear plot saved as hic_map_v2_vs_assembly_v2.pdf")