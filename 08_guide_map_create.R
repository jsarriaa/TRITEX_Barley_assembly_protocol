# From the Morex reference create a guide map
rscript_path <- "/scratch/software-phgv2/miniconda3/envs/IPK-tritex/bin/R"

library(data.table)

# Set up variables
f<-'/scratch/GDB136/IPK/MorexV3_singlecopy_100bp.bed'
fread(f) -> p
setnames(p, c("agp_chr", "agp_start", "agp_end", "seq"))

# Change chromosome names to their respective chr numbers

# Create integer chromosome column (do not overwrite agp_chr)
p[, agp_chr_num := as.integer(gsub("[^0-9]", "", agp_chr))]

# This is for unplaced scaffolds, assigning them as "NA"
p[grepl("NW", agp_chr), agp_chr:=as.integer(NA)]

# Build the guide map
p[, .(css_contig = seq,
      popseq_cM = (agp_start + agp_end) / 2 / 1e6,
      sorted_genome = "H",
      css_contig_length = agp_end - agp_start,
      popseq_alphachr = sub("^chr", "", agp_chr),      # REMOVE "chr" prefix
      sorted_arm = NA,
      popseq_chr = agp_chr_num,       # integer version
      sorted_alphachr = sub("^chr", "", agp_chr),      # REMOVE "chr" prefix
      sorted_chr = agp_chr_num,     
      sorted_lib = sub("^chr", "", agp_chr))] -> pp    # REMOVE "chr" prefix

# Save as RDS file for use in TRITEX pipeline
saveRDS(pp, "GDB136_pseudopopseq.Rds")

print("Guide map columns relevant for plotting:")
print(pp[, .(css_contig, pos=popseq_cM, sorted_chr, sorted_alphachr, popseq_chr, popseq_alphachr)])

# Save the guide map as a CSV file
write.csv(pp, "GDB136_pseudopopseq.csv", row.names = TRUE)