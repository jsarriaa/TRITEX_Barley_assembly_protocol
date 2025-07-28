library(data.table)
library(parallel)
library(igraph)
library(stringi)
library(zoo)

cat("Disabling multi-threading and optimizations of recent data.table versions\n")
cat("setDTthreads(1)\n")
cat("options(data.table.optimize=1)\n")
setDTthreads(1)
options(data.table.optimize=1)

# Define mapping table between numeric and proper chromosome names (21 -> 7D in Triticum aestivum)
chrNames<-function(agp=F, species="wheat") {
 if(species == "wheat"){
  data.table(alphachr=apply(expand.grid(1:7, c("A", "B", "D"), stringsAsFactors=F), 1, function(x) paste(x, collapse="")), chr=1:21)->z
 } else if (species == "barley"){
  data.table(alphachr=apply(expand.grid(1:7, "H", stringsAsFactors=F), 1, function(x) paste(x, collapse="")), chr=1:7)->z
 }
 else if (species == "rye"){
  data.table(alphachr=apply(expand.grid(1:7, "R", stringsAsFactors=F), 1, function(x) paste(x, collapse="")), chr=1:7)->z
 }
 else if (species == "lolium"){
  data.table(alphachr=as.character(1:7), chr=1:7)->z
 }
 else if (species == "maize"){
  data.table(alphachr=as.character(1:10), chr=1:10)->z
 }
 else if (species == "sharonensis"){
  data.table(alphachr=apply(expand.grid(1:7, "S", stringsAsFactors=F), 1, function(x) paste(x, collapse="")), chr=1:7)->z
 }
 else if (species == "oats"){
  data.table(alphachr=sub(" ", "", apply(expand.grid(1:21, "M", stringsAsFactors=F), 1, function(x) paste(x, collapse=""))), chr=1:21)->z
 }
 else if (species == "oats_new"){
  data.table(alphachr=apply(expand.grid(1:7, c("A", "C", "D"), stringsAsFactors=F), 1, function(x) paste(x, collapse="")), chr=1:21)->z
 }
 else if (species == "avena_barbata"){
  data.table(alphachr=apply(expand.grid(1:7, c("A", "B"), stringsAsFactors=F), 1, function(x) paste(x, collapse="")), chr=1:14)->z
 }
 else if (species == "hordeum_bulbosum"){
  data.table(alphachr=sort(apply(expand.grid(1:7, c("H_"), c(1:2), stringsAsFactors=F), 1, function(x) paste(x, collapse=""))), chr=1:14)->z
 }
 else if (species == "faba_bean"){
  data.table(alphachr=apply(expand.grid(1:6, "", stringsAsFactors=F), 1, function(x) paste(x, collapse="")), chr=1:6)->z
 }
 else if (species == "Pginseng"){
  data.table(alphachr=as.character(1:24), chr=1:24)->z
 }
 else if (species == "pea"){
 data.table(alphachr=apply(expand.grid(1:7, "", stringsAsFactors=F), 1, function(x) paste(x, collapse="")), chr=1:7)->z
 }
 else if(species == "rapeseed"){
  expand.grid(1:10, c("A", "C"), stringsAsFactors=F) -> rs_chr
  rs_chr <- rs_chr[-nrow(rs_chr), ]
  rs_chr <- rs_chr[, c(2, 1)]
  rs_chr <- data.frame(chr = sprintf("%s%s", rs_chr$Var2, rs_chr$Var1))
  data.table(alphachr=apply(rs_chr, 1, function(x) paste(x, collapse="", sep="")), chr=1:19)->z
 }
if(agp){
  rbind(z, data.table(alphachr="Un", chr=0))[, agp_chr := paste0("chr", alphachr)]->z
 }
 z[]
}

# alias for backwards compatibility
wheatchr <- chrNames

# Calculate N50 (and similar statistics) from a vector scaffold lengths
n50<-function(l, p=0.5){
 l[order(l)] -> l
 l[head(which(cumsum(as.numeric(l)) >= (1 - p) * sum(as.numeric(l))), n=1)]
}

# Use POSPEQ, Hi-C and flow-sorting infromation to assign scaffolds to approximate chromosomal locations
anchor_scaffolds<-function(assembly, popseq, species=NULL,  
			   sorted_percentile=95,
			   popseq_percentile=90,
			   hic_percentile=98){

 if(is.null(species)){
  stop("Parameter 'species' is NULL. Please set 'species' to one of \"wheat\", \"barley\", \"Pginseng\", \"oats\", \"oats_new\", \"avena_barbata\", \"faba bean\", \"lolium\", \"sharonensis\", \"hordeum_bulbosum\" , \"maize\", \"rye\", \"pea\".")
 }

 if(!species %in% c("wheat", "barley", "avena_barbata", "hordeum_bulbosum", "faba_bean",  "Pginseng", "rye", "oats", "oats_new", "sharonensis", "lolium", "maize", "pea")){
  stop("Parameter 'species' is NULL. Please set 'species' to one of \"wheat\", \"barley\", \"oats\", \"oats_new\", \"faba_bean\", \"Pginseng\", \"lolium\", \"sharonensis\", \"hordeum_bulbosum\", \"avena_barbata\", \"maize\", \"pea\", \"rye\".")
 }

 setnames(chrNames(species=species), c("alphachr", "chr"), c("popseq_alphachr", "popseq_chr"))->wheatchr

 assembly$info -> fai
 assembly$cssaln -> cssaln
 if(is.null(assembly$fpairs) || nrow(assembly$fpairs) == 0){
  hic <- F
 } else {
  hic <- T
  assembly$fpairs -> fpairs
 }

 # Assignment of CARMA chromosome
 cssaln[!is.na(sorted_alphachr), .N, keyby=.(scaffold, sorted_alphachr)]->z
 z[order(-N)][, .(Ncss=sum(N), sorted_alphachr=sorted_alphachr[1], sorted_Ncss1=N[1],
			 sorted_alphachr2=sorted_alphachr[2], sorted_Ncss2=N[2]), keyby=scaffold]->z
 z[, sorted_pchr := sorted_Ncss1/Ncss]
 z[, sorted_p12 := sorted_Ncss2/sorted_Ncss1]

 # Assignment of CARMA chromosome arm
 cssaln[sorted_arm == "L", .(NL=.N), keyby=.(scaffold, sorted_alphachr)]->al
 cssaln[sorted_arm == "S", .(NS=.N), keyby=.(scaffold, sorted_alphachr)]->as
 al[z, on=c("scaffold", "sorted_alphachr")]->z
 as[z, on=c("scaffold", "sorted_alphachr")]->z
 z[is.na(NL), NL := 0]
 z[is.na(NS), NS := 0]
 z[, sorted_arm := ifelse(NS >=  NL, "S", "L")]
 z[sorted_alphachr %in% c("1H", "3B"), sorted_arm := NA]
 z[sorted_arm == "S", sorted_parm := NS/sorted_Ncss1]
 z[sorted_arm == "L", sorted_parm := NL/sorted_Ncss1]
 setnames(copy(wheatchr), c("sorted_alphachr", "sorted_chr"))[z, on="sorted_alphachr"]->z
 setnames(copy(wheatchr), c("sorted_alphachr2", "sorted_chr2"))[z, on="sorted_alphachr2"]->z
 z[fai, on="scaffold"]->info
 info[is.na(Ncss), Ncss := 0]
 info[is.na(NS), NS := 0]
 info[is.na(NL), NL := 0]
 info[is.na(sorted_Ncss1), sorted_Ncss1 := 0]
 info[is.na(sorted_Ncss2), sorted_Ncss2 := 0]

 # Assignment of POPSEQ genetic positions
 popseq[!is.na(popseq_alphachr), .(css_contig, popseq_alphachr, popseq_cM)][cssaln[, .(css_contig, scaffold)], on="css_contig", nomatch=0]->z
 z[, .(.N, popseq_cM=mean(popseq_cM), popseq_cM_sd=ifelse(length(popseq_cM) > 1, sd(popseq_cM), 0), popseq_cM_mad=mad(popseq_cM)), keyby=.(scaffold, popseq_alphachr)]->zz
 zz[, popseq_Ncss := sum(N), by=scaffold]->zz
 zz[order(-N)][, .(popseq_alphachr=popseq_alphachr[1], popseq_Ncss1=N[1],
			 popseq_alphachr2=popseq_alphachr[2], popseq_Ncss2=N[2]), keyby=scaffold]->x
 zz[, .(scaffold, popseq_alphachr, popseq_Ncss, popseq_cM, popseq_cM_sd,  popseq_cM_mad)][x, on=c("scaffold", "popseq_alphachr")]->zz
 zz[, popseq_pchr := popseq_Ncss1/popseq_Ncss]
 zz[, popseq_p12 := popseq_Ncss2/popseq_Ncss1]
 wheatchr[zz, on="popseq_alphachr"]->zz
 setnames(copy(wheatchr), paste0(names(wheatchr), 2))[zz, on="popseq_alphachr2"]->zz
 zz[info, on="scaffold"]->info
 info[is.na(popseq_Ncss), popseq_Ncss := 0]
 info[is.na(popseq_Ncss1), popseq_Ncss1 := 0]
 info[is.na(popseq_Ncss2), popseq_Ncss2 := 0]

 # Chromosome assignment based on Hi-C links
 if(hic){
  info[!(popseq_chr != sorted_chr)][, .(scaffold, chr=popseq_chr)]->info0
  setnames(copy(info0), names(info0), sub("$", "1", names(info0)))[fpairs, on="scaffold1"]->tcc_pos
  setnames(copy(info0), names(info0), sub("$", "2", names(info0)))[tcc_pos, on="scaffold2"]->tcc_pos
  tcc_pos[!is.na(chr1), .N, key=.(scaffold=scaffold2, hic_chr=chr1)]->z
  z[order(-N)][, .(Nhic=sum(N), hic_chr=hic_chr[1], hic_N1=N[1],
		   hic_chr2=hic_chr[2], hic_N2=N[2]), keyby=scaffold]->zz
  zz[, hic_pchr := hic_N1/Nhic]
  zz[, hic_p12 := hic_N2/hic_N1]
  zz[info, on="scaffold"]->info
  info[is.na(Nhic), Nhic := 0]
  info[is.na(hic_N1), hic_N1 := 0]
  info[is.na(hic_N2), hic_N2 := 0]
 }

 if(hic){
  measure <- c("popseq_chr", "hic_chr", "sorted_chr")
 } else {
  measure <- c("popseq_chr", "sorted_chr")
 }
 melt(info, id.vars="scaffold", measure.vars=measure, variable.factor=F, variable.name="map", na.rm=T, value.name="chr")->w
 w[, .N, key=.(scaffold, chr)]->w
 w[order(-N), .(Nchr_ass = sum(N), Nchr_ass_uniq = .N), keyby=scaffold]->w
 w[info, on="scaffold"]->info
 info[is.na(Nchr_ass), Nchr_ass := 0]
 info[is.na(Nchr_ass_uniq), Nchr_ass_uniq := 0]

 x<-info[Ncss >= 30, quantile(na.omit(sorted_p12), 0:100/100)][sorted_percentile+1]
 info[, bad_sorted := (sorted_p12 >= x & sorted_Ncss2 >= 2)]
 x<-info[popseq_Ncss >= 30, quantile(na.omit(popseq_p12), 0:100/100)][popseq_percentile+1]
 info[, bad_popseq := (popseq_p12 >= x & popseq_Ncss2 >= 2)]
 
 info[is.na(bad_sorted), bad_sorted := F]
 info[is.na(bad_popseq), bad_popseq:= F]

 if(hic){
  x<-info[Nhic >= 30, quantile(na.omit(hic_p12), 0:100/100)][hic_percentile+1]
  info[Nhic >= 30, bad_hic := hic_p12 >= x & hic_N2 >= 2]
  info[is.na(bad_hic), bad_hic := F]
 }

 melt(info, id.vars="scaffold", measure.vars=grep(value=T, "bad_", names(info)), variable.factor=F, variable.name="map", na.rm=T, value.name="bad")[bad == T]->w
 w[, .(Nbad=.N), key=scaffold]->w
 w[info, on="scaffold"]->info
 info[is.na(Nbad), Nbad := 0]

 assembly$info <- info
 assembly$popseq <- popseq
 if(hic){
  assembly$fpairs <- tcc_pos
 }
 assembly
}

# Create diagnostics plot for putative chimeras
plot_chimeras<-function(scaffolds, assembly, breaks=NULL, file, mindist=0, cores=20, species="wheat", mbscale=F, autobreaks=T, refname="Morex V2"){

 plot_popseq_carma_tcc<-function(scaffold, breaks=NULL, page, info, cssaln, tcc_pos, span, molcov, species){
  chrNames(species=species)->wheatchr

  i <- scaffold

  if(is.null(tcc_pos) || nrow(tcc_pos) == 0){
   hic <- F
  } else {
   hic <- T
  }

  if(is.null(molcov) || nrow(molcov) == 0){
   tenex <- F
  } else {
   tenex <- T
  }

  nrow <- 1
  if(hic){
   nrow <- nrow + 1
  }
  if(tenex){
   nrow <- nrow + 1
  }

  par(mfrow=c(nrow,3))

  par(oma=c(1,0,3,0))
  if(hic){
   par(mar=c(1,4,4,1))
  } else {
   par(mar=c(4,4,4,1))
  }
  par(cex=1)

  if(!is.null(breaks)){
   br <- breaks[i, on="scaffold", br]/1e6
  } else if (hic & autobreaks){
   span[i, on="scaffold"][order(r)][1, bin]/1e6->br
  } else {
   br <- NULL
  }

  if(hic){
   info[scaffold == i, paste0(", bad Hi-C: ", bad_hic)] -> badhic
   xlab<-""
  } else {
   badhic <- ""
   xlab<-"position in scaffold (Mb)"
  }

  cssaln[i, on="scaffold"]->z
  l=info[i, on="scaffold"][, length]/1e6

  if(species %in% c("rye", "lolium", "barley")){
   ymin <- 8
  } else {
   ymin <- 24
  }

  if(mbscale == F){
   ylab = "POPSEQ chromosome"
  } else {
   ylab = paste(refname, "chromosome")
  }

  plot(xlim=c(0,l), 0, ylim=c(ymin,1), bty='l', col=0, yaxt='n', pch=20, xlab=xlab, ylab=ylab)
  if(!is.null(br)){
   abline(v=br, lwd=3, col="blue")
  }
  if(mbscale == F){
   title("POPSEQ", line=1)
  } else {
   title(refname, line=1)
  }
  z[, points(pos/1e6, popseq_chr, pch=20)]
  axis(2, las=1, wheatchr$chr, wheatchr$alphachr)

  if(mbscale == F){
   info[scaffold == i, title(outer=T, line=0, paste("page: ", page, ", ", scaffold, sep="", " (", round(l,1), " Mb),\n",  
			   "bad POPSEQ: ", bad_popseq, ", bad CARMA: ", bad_sorted, badhic))]
  } else {
   info[scaffold == i, title(outer=T, line=0, paste("page: ", page, "\n", scaffold, sep="", " (", round(l,1), " Mb)"))]
  }

  w<-info[i, on="scaffold"]$popseq_chr
  if(is.na(w)){
   plot(xlab='', ylab='', type='n', 0, axes=F)
  } else {
   if(mbscale == F){
    ylab = "genetic position (cM)"
   } else {
    ylab = "pseudomolecule position (Mb)"
   }
   z[popseq_chr == w][, plot(pch=20, xlab=xlab, ylab=ylab, xlim=c(0,l), bty='l', las=1, pos/1e6, popseq_cM)]
   if(mbscale == F){
    title(main="Genetic positions of CSS contigs\nfrom major chromosome")
   } else {
    title(main="Physical positions of 1 kb single-copy tags\nfrom major chromosome")
   }
   if(hic){
    abline(v=br, lwd=3, col="blue")
   }
  } 

  if(!species %in% c('lolium', "oats", "oats_new")){
   if(mbscale == F){
    ylab="flow-sorting chromosome"
    ttt="flow-sorting CARMA"
   } else {
    ylab="chromsome arm"
    ttt="Hi-C based chromosome arm assignment"
   }
   plot(xlim=c(0,l), 0, bty='l', ylim=c(ymin,1), yaxt='n', pch=20, col=0, xlab=xlab, ylab=ylab)
   title(ttt, line=1)
   if(!is.null(br)){
    abline(v=br, lwd=3, col="blue")
   }
   if(species == "wheat"){
    z[, points(pos/1e6, sorted_chr, col=ifelse(sorted_alphachr == "3B" | sorted_arm == "S", 1, 2), pch=20)]
    legend(horiz=T, "bottomleft", pch=19, col=1:2, bg="white", legend=c("short arm / 3B", "long arm"))
   }
   if(species == "rye"){
    z[, points(pos/1e6, sorted_chr, col=1, pch=20)]
   }
   if(species == "barley"){
    z[, points(pos/1e6, sorted_chr, col=ifelse(sorted_alphachr == "1H" | sorted_arm == "S", 1, 2), pch=20)]
    if(mbscale == F){
     llx =  "short arm / 1H"
    } else {
     llx = "short arm"
    }
    legend(horiz=T, "bottomleft", pch=19, col=1:2, bg="white", legend=c(llx, "long arm"))
   }
   axis(2, las=1, wheatchr$chr, wheatchr$alphachr)
  } else {
   plot(0, type='n', xlab="", ylab="", axes=F)
  }

  if(hic){
   par(mar=c(4,4,4,1))
   plot(0, xlim=c(0,l), col=0, bty='l', ylab='Hi-C chromosome', yaxt='n', ylim=c(ymin,1), xlab="position in scaffold (Mb)")
   title("Interchromosomal Hi-C links", line=1)
   tcc_pos[i, on="scaffold1", points(pos1/1e6, col="#00000003", chr2, pch=20)]
   abline(v=br, lwd=3, col="blue")
   axis(2, las=1, wheatchr$chr, wheatchr$alphachr)

   span[i, on="scaffold"]->w
   if(nrow(w[!is.na(n)]) == 0){
    plot(xlab='', ylab='', type='n', 0, axes=F)
   } else {
    w[order(bin), plot(bty='l', bin/1e6, n, log='y', xlim=c(0,l), col=0, ylab='coverage', las=1, xlab="position in scaffold (Mb)")]
    title("Intrascaffold physical Hi-C coverage", line=1)
    abline(v=br, lwd=3, col="blue")
    w[order(bin), points(bin/1e6, n, xlim=c(0,l), type='l', lwd=3, col=1)]
   }

   if(nrow(w[!is.na(n)]) == 0){
    plot(xlab='', ylab='', type='n', 0, axes=F)
   } else {
    w[order(bin), plot(bty='l', bin/1e6, r, xlim=c(0,l), col=0, ylab='log2(observed/expected ratio)', las=1, xlab="position in scaffold (Mb)")]
    title("Hi-C expected vs. observed coverage", line=1)
    abline(v=br, lwd=3, col="blue")
    w[order(bin), points(bin/1e6, r, xlim=c(0,l), type='l', lwd=3, col=1)]
   }
  }

  if(tenex){
   molcov[i, on="scaffold"]->w
   if(nrow(w) == 0){
    plot(xlab='', ylab='', type='n', 0, axes=F)
   } else {
    w[order(bin), plot(bty='l', bin/1e6, n, log='y', xlim=c(0,l), col=0, ylab='coverage', las=1, xlab="position in scaffold (Mb)")]
    title("10X molecule coverage", line=1)
    if(!is.null(br)){
     abline(v=br, lwd=3, col="blue")
    }
    w[order(bin), points(bin/1e6, n, xlim=c(0,l), type='l', lwd=3, col=1)]
   }

   if(nrow(w) == 0){
    plot(xlab='', ylab='', type='n', 0, axes=F)
   } else {
    w[order(bin), plot(bty='l', bin/1e6, r, xlim=c(0,l), col=0, ylab='log2(observed/expected ratio)', las=1, xlab="position in scaffold (Mb)")]
    title("10X expected vs. observed coverage", line=1)
    if(!is.null(br)){
     abline(v=br, lwd=3, col="blue")
    }
    w[order(bin), points(bin/1e6, r, xlim=c(0,l), type='l', lwd=3, col=1)]
   }
  }
 }

 info<-assembly$info
 ff<-assembly$cov
 tcc_pos <- assembly$fpairs
 cssaln <- assembly$cssaln
 molcov <- assembly$molecule_cov

 width <- 2500
 height <- 1000
 res <- 150

 if(!is.null(tcc_pos) && nrow(tcc_pos) > 0){
  height <- height + 1000
 }
 if(!is.null(molcov) && nrow(molcov) > 0){
  height <- height + 1000
 }

 scaffolds <- scaffolds[!duplicated(scaffold), .(scaffold)]

 bad <- copy(scaffolds)
 out <- file
 bad[, f:=tempfile(fileext=".png"), by=scaffold]
 bad[, p:=tempfile(fileext=".pdf"), by=scaffold]
 mclapply(mc.cores=cores, 1:nrow(bad), function(i){
  file <- bad[i, f]
  pdf <- bad[i, p]
  s<-bad[i,scaffold]
  cat(paste0(i, " ", bad[i, scaffold], "\n"))
  png(file, height=height, res=res, width=width)
  plot_popseq_carma_tcc(s, breaks=breaks, page=i, info, cssaln, tcc_pos, ff[d >= mindist], molcov, species)
  dev.off()
  system(paste("convert", file, pdf))
  unlink(file)
 })
 system(paste0("gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=", out, " ", paste(bad$p, collapse=' ')))
 unlink(bad$p)
}

# Calculate the physical coverage with 10X molecules in sliding windows along the scaffolds
add_molecule_cov<-function(assembly, scaffolds=NULL, binsize=200, cores=1){

 info <- assembly$info

 if(is.null(assembly$molecules) || nrow(assembly$molecules) == 0){
  stop("The assembly object does not have a molecule table; aborting.")
 }

 if("mr_10x" %in% colnames(info)){
  stop("assembly$info already has mr_10x column; aborting.")
 }
 
 if(is.null(scaffolds)){
  scaffolds <- info$scaffold
  copy(assembly$molecules) -> mol
  null=T
 } else {
  info[scaffolds, on="scaffold"] -> info
  assembly$molecules[scaffolds, on="scaffold"] -> mol
  null=F
 }

 mol -> f
 f[, bin1 := start %/% binsize * binsize]
 f[, bin2 := end %/% binsize * binsize]
 f[bin2 - bin1 > 2 * binsize] -> f
 setkey(f, scaffold)
 f[, i := 1:.N]

 rbindlist(mclapply(mc.cores=cores, unique(f$scaffold), function(j){
  f[j][, .(scaffold, bin=seq(bin1+binsize, bin2-binsize, binsize)), key=i][, .(n=.N), key=.(scaffold, bin)]
 })) -> ff

 if(nrow(ff) > 0){
  info[, .(scaffold, length)][ff, on="scaffold"]->ff
  ff[, d := pmin(bin, (length-bin) %/% binsize * binsize)]
  ff[, nbin := .N, key="scaffold"]
  ff[, mn := mean(n), key=d]
  ff[, r := log2(n/mn)]

  ff[, .(mr_10x = min(r)), key=scaffold][info, on="scaffold"] -> info_mr
 } else {
  copy(info) -> info_mr
  info_mr[, mr_10x := NA]
 }

 if(null){
  assembly$info <- info_mr
  assembly$molecule_cov <- ff
  assembly$mol_binsize <- binsize

  assembly
 } else {
  list(info=info_mr, molecule_cov=ff)
 }
}

# Read the files as produced by run_10x_mapping.zsh
read_10x_molecules<-function(files, ncores=1){
 f <- files
 rbindlist(mclapply(mc.cores=ncores, names(f), function(i){
  fread(cmd=paste('gzip -cd', f[i]), head=F, col.names=c("scaffold", "start", "end", "barcode", "npairs"))->z
  z[, sample := i]
 })) -> mol
 mol[, length := end - start]
 mol[, start := start + 1]
 mol[]
}

# Initialized an assembly object with genetic map/flowsorting information, 10X and Hi-C links
init_assembly<-function(fai, cssaln, fpairs=NULL, molecules=NULL, rename=NULL){

 copy(fai)->info
 info[, orig_start := 1]
 info[, orig_end := length]
 if(is.null(rename)) {
  info[, orig_scaffold := scaffold]
 } else {
  setnames(info, "scaffold", "orig_scaffold")
  rename[info, on="orig_scaffold"] -> info
 }

 copy(cssaln)->z
 if(is.null(rename)) {
  z[, orig_scaffold := scaffold] 
 } else {
  setnames(z, "scaffold", "orig_scaffold")
  rename[z, on="orig_scaffold"] -> z
 }
 z[, orig_pos := pos]
 z[, orig_scaffold_length := scaffold_length]

 if(!is.null(molecules)){
  copy(molecules)->y
  y[, orig_scaffold := scaffold]
  y[, orig_start := start]
  y[, orig_end := end]
  y[, scaffold := NULL]
  info[, .(orig_scaffold, scaffold)][y, on="orig_scaffold"]->y
 } else {
  y <- data.table()
 }

 if(!is.null(fpairs)){
  copy(fpairs)->tcc
  tcc[, orig_scaffold1 := scaffold1]
  tcc[, orig_pos1 := pos1]
  tcc[, orig_scaffold2 := scaffold2]
  tcc[, orig_pos2 := pos2]
  tcc[, scaffold1 := NULL]
  tcc[, scaffold2 := NULL]
  info[, .(orig_scaffold1=orig_scaffold, scaffold1=scaffold)][tcc, on="orig_scaffold1"]->tcc
  info[, .(orig_scaffold2=orig_scaffold, scaffold2=scaffold)][tcc, on="orig_scaffold2"]->tcc
 } else {
  tcc <- data.table()
 }

 list(info=info, cssaln=z, fpairs=tcc, molecules=y)
}

# Find breakpoints in chimeric scaffolds based on drop in 10X molecule coverage
find_10x_breaks<-function(assembly, scaffolds=NULL, interval = 5e4, minNbin = 20, dist = 5e3, ratio = -3){
 cov <- copy(assembly$molecule_cov)
 if(!is.null(scaffolds)){
  cov[scaffolds, on="scaffold"]->cov
 } 
 cov[, b := bin %/% interval * interval]
 cov[nbin >= minNbin & pmin(bin, length - bin) >= dist & r <= ratio]->e
 if(nrow(e) == 0){
  return(NULL)
 }
 e[order(r)][, idx := 1:.N, by=.(scaffold, b)][idx == 1]->e
 setnames(e, "bin", "br")[]
 e[]
}

# Calculate physical coverage with Hi-C links in sliding windows along the scaffolds
add_hic_cov<-function(assembly, scaffolds=NULL, binsize=1e3, binsize2=1e5, minNbin=50, innerDist=1e5, cores=1){

 print("Inside add_hic_cov")

 info<-assembly$info
 
 # print info columns
 print (colnames(info))

 
 if("mr" %in% colnames(info) | "mri" %in% colnames(info)){
  stop("assembly$info already has mr and/or mri columns; aborting.")
 }

 fpairs<-assembly$fpairs

 # print fpairs columns
 print (colnames(fpairs))

 if(is.null(scaffolds)){
  scaffolds <- info$scaffold
  null=T
 } else {
  info[scaffold %in% scaffolds]->info
  fpairs[scaffold1 %in% scaffolds]->fpairs
  null=F
 }
  
  

 fpairs[scaffold1 == scaffold2 & pos1 < pos2][, .(scaffold = scaffold1, bin1 = pos1 %/% binsize * binsize, bin2 =pos2 %/% binsize * binsize)]->f
 f[bin2 - bin1 > 2*binsize]->f
 f[, i := 1:.N]
 f[, b := paste0(scaffold, ":", bin1 %/% binsize2)]
 setkey(f, b)

 rbindlist(mclapply(mc.cores=cores, unique(f$b), function(j){
 f[j][, .(scaffold=scaffold, bin=seq(bin1+binsize, bin2-binsize, binsize)), key=i][, .(n=.N), key=.(scaffold, bin)]
 }))->ff

 if(nrow(ff) > 0){
  ff[, .(n=sum(n)), key=.(scaffold, bin)]->ff
  info[, .(scaffold, length)][ff, on="scaffold"]->ff
  ff[, d := pmin(bin, (length-bin) %/% binsize * binsize)]
  ff[, nbin := .N, key="scaffold"]
  ff[, mn := mean(n), key=d]
  ff[, r := log2(n/mn)]
  ff[nbin >= minNbin, .(mr=suppressWarnings(min(r))), key=scaffold][order(mr)]->z
  ff[nbin > minNbin & d >= innerDist, .(mri=suppressWarnings(min(r))), key=scaffold]->zi
  z[ff, on="scaffold"]->ff
  zi[ff, on="scaffold"]->ff
  z[info, on="scaffold"]->info_mr
  zi[info_mr, on="scaffold"]->info_mr
 } else {
  copy(info) -> info_mr
  info_mr[, c("mri", "mr") := list(NA, NA)]
 }

 if(null){
  assembly$info=info_mr
  assembly$cov=ff

  assembly$binsize <- binsize
  assembly$minNbin <- minNbin
  assembly$innerDist <- innerDist

  assembly
 } else {
  list(info=info_mr, cov=ff)
 }
}

# Break scaffolds at specified points and lift positional information to updated assembly
break_scaffolds<-function(breaks, assembly, prefix, slop, cores=1, species, 
			  regex1="(^.*[^-0-9])(([0-9]+)(-[0-9]+)?$)", 
			  regex2="(^.*[^-0-9])([0-9]+(-[0-9]+)?$)"){
 info <- assembly$info
 cov <- assembly$cov
 fpairs <- assembly$fpairs
 cssaln <- assembly$cssaln
 molecules <- assembly$molecules
 regex2_0 <- regex2

 br<-copy(breaks)
 info[, .(scaffold, orig_scaffold, old_scaffold=scaffold, orig_start, orig_end, length)]->fai

 cat("Split scaffolds\n")
 j <- 0
 while(nrow(br) > 0){
  j <- j + 1
  o <- nrow(fai)
  fai[br, on="scaffold"] -> br
  br[, orig_br := orig_start + br - 1]
  br[order(scaffold, br)] -> br
  br[duplicated(scaffold)] -> nbr
  br[!duplicated(scaffold)] -> br

  max(as.integer(sub(regex1, "\\3", fai$scaffold)))->maxidx
  br[, idx := 3*1:.N-2]
  br[, scaffold1 := paste0(prefix, maxidx+idx)]
  br[, start1 := 1]
  br[, end1 := pmax(0, br - slop - 1)]
  br[, scaffold2 := paste0(prefix, maxidx+idx+1)]
  br[, start2 := pmax(1, br - slop)]
  br[, end2 := pmin(br + slop - 1, length)]
  br[, scaffold3 := paste0(prefix, maxidx+idx+2)]
  br[, start3 := pmin(length + 1, br + slop)]
  br[, end3 := length]
  br[, length1 := 1 + end1 - start1]
  br[, length2 := 1 + end2 - start2]
  br[, length3 := 1 + end3 - start3]
  rbind(
   br[, .(scaffold=scaffold1, length=length1, orig_scaffold, orig_start=orig_start+start1-1, orig_end=orig_start+end1-1, old_scaffold)],
   br[, .(scaffold=scaffold2, length=length2, orig_scaffold, orig_start=orig_start+start2-1, orig_end=orig_start+end2-1, old_scaffold)],
   br[, .(scaffold=scaffold3, length=length3, orig_scaffold, orig_start=orig_start+start3-1, orig_end=orig_start+end3-1, old_scaffold)],
   fai[!scaffold %in% br$scaffold, .(orig_scaffold, length,  orig_start, orig_end, old_scaffold,
			      scaffold=paste0(prefix, sub(regex2, "\\2", scaffold)))]

  ) -> fai

  cat(paste0("Iteration ", j, " finished. "))
  fai[length > 0]->fai
  cat(paste0("The number of scaffolds increased from ", o, " to ", nrow(fai), ".\n"))
  fai[, .(scaffold, orig_scaffold, orig_start, orig_br=orig_start)][nbr[, .(orig_scaffold, orig_br)], on=c("orig_scaffold", "orig_br"), roll=T]->nbr
  nbr[, br := orig_br - orig_start + 1]
  nbr[, .(scaffold, br)] -> br

  regex1="(^.*[^-0-9])(([0-9]+)(-[0-9]+)?$)"
  regex2="(^.*[^-0-9])([0-9]+(-[0-9]+)?$)"
 }

 fai[, split := F]
 fai[old_scaffold %in% breaks$scaffold, split := T]
 fai[, old_scaffold := NULL]

 assembly_new<-list(info=fai)

 cat("Transpose cssaln\n")
 copy(cssaln) -> z
 z[, scaffold_length := NULL]
 z[, scaffold := NULL]
 fai[, .(scaffold, scaffold_length=length, orig_scaffold, orig_start, orig_pos=orig_start)][z, on=c("orig_scaffold", "orig_pos"), roll=T]->z
 z[, pos := orig_pos - orig_start + 1]
 z[, orig_start := NULL]
 assembly_new$cssaln <- z
 
 if("fpairs" %in% names(assembly) && nrow(fpairs) > 0){
  cat("Transpose fpairs\n")
  assembly$fpairs[, .(orig_scaffold1, orig_scaffold2, orig_pos1, orig_pos2)]->z
  fai[, .(scaffold1=scaffold, orig_scaffold1=orig_scaffold, orig_start1=orig_start, orig_pos1=orig_start)][z, on=c("orig_scaffold1", "orig_pos1"), roll=T]->z
  fai[, .(scaffold2=scaffold, orig_scaffold2=orig_scaffold, orig_start2=orig_start, orig_pos2=orig_start)][z, on=c("orig_scaffold2", "orig_pos2"), roll=T]->z
  z[, pos1 := orig_pos1 - orig_start1 + 1]
  z[, pos2 := orig_pos2 - orig_start2 + 1]
  z[, orig_start1 := NULL]
  z[, orig_start2 := NULL]
  assembly_new$fpairs <- z
 } else {
  assembly_new$fpairs <- data.table()
 }

 if("molecules" %in% names(assembly) && nrow(molecules) > 0){
  cat("Transpose molecules\n")
  copy(molecules) -> z
  z[, scaffold := NULL]
  fai[, .(scaffold, orig_scaffold, orig_start, s_length=orig_end - orig_start + 1, orig_pos=orig_start)][z, on=c("orig_scaffold", "orig_start"), roll=T]->z
  z[, start := orig_start - orig_pos + 1]
  z[, end := orig_end - orig_pos + 1]
  z[end <= s_length]->z
  z[, orig_pos := NULL]
  z[, s_length  := NULL]
  assembly_new$molecules <- z
 } else {
  assembly_new$molecules <- data.table()
 }

 assembly_new$breaks <- breaks

 cat("Anchor scaffolds\n")
 anchor_scaffolds(assembly_new, popseq=assembly$popseq, species=species)->assembly_new
 
 if("mr_10x" %in% names(assembly_new$info)){
  assembly_new$info[, mr_10x := NULL]
 }

 if("mr" %in% names(assembly_new$info)){
  assembly_new$info[, mr := NULL]
  assembly_new$info[, mri := NULL]
 }

 if("cov" %in% names(assembly) & nrow(fpairs) > 0){
  cat("Hi-C coverage\n")
  add_hic_cov(assembly_new, scaffolds=fai[split == T]$scaffold, binsize=assembly$binsize, minNbin=assembly$minNbin, innerDist=assembly$innerDist, cores=cores)->cov
  assembly$cov[!scaffold %in% breaks$scaffold]->x
  x[, scaffold:=paste0(prefix, sub(regex2_0, "\\2", scaffold))]
  if(nrow(cov$cov) > 0){
   rbind(x, cov$cov)->assembly_new$cov
  } else {
   x -> assembly_new$cov
  }
  info[!scaffold %in% breaks$scaffold]->x
  x[, scaffold:=paste0(prefix, sub(regex2_0, "\\2", scaffold))]
  x[, split := F]
  rbind(x[, names(cov$info), with=F], cov$info)->assembly_new$info
 } else {
  assembly_new$cov <- data.table()
 }

 if("molecule_cov" %in% names(assembly) & nrow(molecules) > 0){
  cat("10X molecule coverage\n")
  add_molecule_cov(assembly_new, scaffolds=fai[split == T]$scaffold, binsize=assembly$mol_binsize, cores=cores)->cov
  info[!breaks$scaffold, on="scaffold"]->x
  x[, scaffold := paste0(prefix, sub(regex2_0, "\\2", scaffold))]
  x[, split := F]
  rbind(x[, names(cov$info), with=F], cov$info)->assembly_new$info
  assembly_new$mol_binsize <- assembly$mol_binsize

  assembly$molecule_cov[!breaks$scaffold, on="scaffold"]->x
  x[, scaffold := paste0(prefix, sub(regex2_0, "\\2", scaffold))]
  if(nrow(cov$molecule_cov) > 0){
   rbind(x, cov$molecule_cov)->assembly_new$molecule_cov
  } else {
   x -> assembly_new$molecule_cov 
  }
 } else {
  assembly_new$molecule_cov <- data.table()
 }

 assembly_new$binsize <- assembly$binsize
 assembly_new$innerDist <- assembly$innerDist
 assembly_new$minNbin <- assembly$minNbin

 assembly_new
}

# Read BAM files with alignment of genetic markers (originally Wheat CSS contigs, the "marker sequences" of the POPSEQ map). Merge with POPSEQ data. Extract flow-sorting information from contig names/
# Specific to bread wheat and CSS POPSEQ map.
read_cssaln<-function(bam, popseq, fai, minqual=30, minlen=1000){
 cssaln <- fread(cmd=paste("samtools view -q ", minqual, "-F260", bam, "| cut -f 1,3,4"), col.names=c("css_contig", "scaffold", "pos")) 

 cssaln[, sorted_lib := sub("_[0-9]+$", "", sub("ta_contig_", "", css_contig))]
 cssaln[, sorted_alphachr := sub("L|S$", "", sorted_lib)]
 cssaln[, sorted_subgenome := sub("[1-7]", "", sorted_alphachr)]
 cssaln[, sorted_arm := sub("[1-7][A-D]", "", sorted_lib)]

 popseq[cssaln, on="css_contig"]->cssaln

 chrNames(species="wheat")->wheatchr
 setnames(copy(wheatchr), c("popseq_alphachr", "popseq_chr"))[cssaln, on="popseq_alphachr"] -> cssaln
 setnames(copy(wheatchr), c("sorted_alphachr", "sorted_chr"))[cssaln, on="sorted_alphachr"] -> cssaln
 fai[, .(scaffold, scaffold_length = length)][cssaln, on="scaffold"]->cssaln
 cssaln[css_contig_length >= minlen]
}

# Genetic function for reading minimap2 alignment of genetic maps and merging them genetic positional infromation. First used for barley. cv Morex, hence the name.
read_morexaln_minimap<-function(paf, popseq, minqual=30, minlen=500, prefix=T){
 fread(cmd=paste0("zgrep tp:A:P ", paf, " | awk -v l=", minlen," -v q=", minqual, " '$2 >= l && $12 >= q' | cut -f 1,6-8"),
              head=F, col.names=c("css_contig", "scaffold", "scaffold_length", "pos"))->z
 if(prefix){
  z[, css_contig := paste0("morex_", css_contig)]
 }
 popseq[z, on="css_contig"]
}

read_guide_map <- read_morexaln_minimap

# Generic function to read PAF files
read_paf<-function(file, primary_only=T, save=F){
 fread(head=F, cmd=paste("zcat", file, "| cut -f -13"), 
       col.names=c("query", "query_length", "query_start", "query_end",
		   "orientation", "reference", "reference_length",
		   "reference_start", "reference_end", "matches", "alnlen", "mapq", "type"))->z
 z[, c("query_start", "query_end", "reference_start", "reference_end") := list(query_start + 1, query_end + 1, reference_start + 1, reference_end + 1)]
 z[, orientation := as.integer(ifelse(orientation == "+", 1, -1))]
 z[, type := sub("tp:A:", "", type)]
 if(primary_only){
  z["P", on="type"]->z
 }
 if(save){
  saveRDS(z, file=sub("paf.gz$", "Rds", file))
 }
 z[]
}

# Aggregate alignments for each (reference, query) pair
summarize_paf<-function(paf){
 setnames(dcast(paf[, .(l=sum(alnlen)), key=.(query, reference, orientation)], query + reference ~ orientation, value.var="l", fill=0), c("-1", "1"), c("lenrev", "lenfw"))->fr
 paf[, .(query_start=min(query_start), query_end=max(query_end),
       reference_start=min(reference_start), reference_end=max(reference_end),
       matches=sum(matches), alnlen=sum(alnlen), naln=.N), 
    key=.(query, query_length, reference, reference_length)]->paf_summary
 fr[paf_summary, on=c("query", "reference")]->paf_summary
 setorder(paf_summary, -alnlen)
 paf_summary[, idx := paste(1:.N)][]
}

# Construct super-scaffold from 10X links. Use heuristics based on genetic map information to prune erroneous edges in the scaffold graph.
scaffold_10x <- function(assembly, prefix="super", min_npairs=5, max_dist=1e5, min_nmol=6, min_nsample=2, popseq_dist=5, max_dist_orientation=5, ncores=1, verbose=T, raw=F, unanchored=T){

 make_agp<-function(membership, gap_size=100){

  membership[, .(scaffold, length, super, bin, orientation)]->z

  setorder(z, super, bin, -length)
  z[, index := 2*1:.N-1]
  z[, gap := F]
  rbind(z, data.table(scaffold="gap", gap=T, super=z$super, bin=NA, length = gap_size,orientation = NA, index=z$index+1))->z
  z[order(index)][, head(.SD, .N-1), by=super]->z
  z[, n := .N, key=super]
  z[n > 1, super_start := cumsum(c(0, length[1:(.N-1)])) + 1, by = super]
  z[n == 1, super_start := 1]
  z[, super_end := cumsum(length), by = super]
  z[, n := NULL]

  z[, .(scaffold=scaffold, bed_start=0, bed_end=length, name=scaffold, score=1, strand=ifelse(is.na(orientation) | orientation == 1, "+", "-"), super=super)]->agp_bed

  list(agp=z, agp_bed=agp_bed)
 }

 make_super_scaffolds<-function(links, info, excluded=c(), ncores, prefix){
  info[, .(scaffold, chr=popseq_chr, cM=popseq_cM, length)]->info2
  excluded -> excluded_scaffolds
  info2[, excluded := scaffold %in% excluded_scaffolds] -> input
  setnames(input, "scaffold", "cluster")
  setnames(copy(links), c("scaffold1", "scaffold2"), c("cluster1", "cluster2"))->hl
  make_super(hl=hl, cluster_info=input, verbose=F, prefix=prefix, cores=ncores, paths=T, path_max=0, known_ends=F, maxiter=100)->s
  copy(s$membership) -> m
  setnames(m, "cluster", "scaffold")

  max(as.integer(sub(paste0(prefix, "_"), "", s$super_info$super))) -> maxidx
  rbind(m, info[!m$scaffold, on="scaffold"][, .(scaffold, bin=1, rank=0, backbone=T, chr=popseq_chr, cM=popseq_cM,
					length=length, excluded=scaffold %in% excluded_scaffolds, super=paste0(prefix, "_", maxidx + 1:.N))])->m

  m[, .(n=.N, nbin=max(bin), max_rank=max(rank), length=sum(length)), key=super]->res
  res[, .(super, super_size=n, super_nbin=nbin)][m, on="super"]->m
  list(membership=m, info=res)
 }

 if(verbose){
  cat("Finding links.\n")
 }
 assembly$info -> info
 assembly$molecules[npairs >= min_npairs] -> z
 info[, .(scaffold, scaffold_length=length)][z, on="scaffold"]->z
 z[end <= max_dist | scaffold_length - start <= max_dist]->z
 z[, nsc := length(unique(scaffold)), key=.(barcode, sample)]
 z[nsc >= 2]->z
 z[, .(scaffold1=scaffold, npairs1=npairs, pos1=as.integer((start+end)/2), sample, barcode)]->x 
 z[, .(scaffold2=scaffold, npairs2=npairs, pos2=as.integer((start+end)/2), sample, barcode)]->y
 y[x, on=.(sample, barcode), allow.cartesian=T][scaffold1 != scaffold2]->xy
 xy -> link_pos
 xy[, .(nmol=.N), key=.(scaffold1, scaffold2, sample)]->w
 w[nmol >= min_nmol]->ww
 ww[, .(nsample = length(unique(sample))), key=.(scaffold1, scaffold2)][nsample >= min_nsample]->ww2
 info[, .(scaffold1=scaffold, popseq_chr1=popseq_chr, length1=length, popseq_pchr1=popseq_pchr, popseq_cM1=popseq_cM)][ww2, on="scaffold1"]->ww2
 info[, .(scaffold2=scaffold, popseq_chr2=popseq_chr, length2=length, popseq_pchr2=popseq_pchr, popseq_cM2=popseq_cM)][ww2, on="scaffold2"]->ww2
 ww2[, same_chr := popseq_chr2 == popseq_chr1]
 ww2[, weight := -1 * log10((length1 + length2) / 1e9)]

 if(popseq_dist > 0){
  ww2[(popseq_chr2 == popseq_chr1 & abs(popseq_cM1 - popseq_cM2) <= popseq_dist)] -> links
 } else {
  ww2 -> links
 }

 ex <- c()

 run <- T

 if(verbose){
  cat("Finding initial super-scaffolds.\n")
 }
 # remove nodes until graph is free of branches of length > 1
 if(!raw){
  while(run){
   make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
   copy(out$membership) -> m
   copy(out$info) -> res
   m[unique(m[rank > 1][, .(super, bin)]), on=c("super", "bin")]->a
   a[rank == 0] -> add
   if(nrow(add) == 0){
    run <- F
   } else{
    c(ex, add$scaffold) -> ex
   }
  }
 } else {
  make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
  out$m -> m
  out$info -> res
 }

 if(!raw){
  if(verbose){
   cat("Tip removal.\n")
  }
  # remove short tips of rank 1 
  links[!scaffold2 %in% ex][, .(degree=.N), key=.(scaffold=scaffold1)]->b
  a <- b[m[m[rank == 1][, .(super=c(super, super, super), bin=c(bin, bin-1, bin+1))], on=c("super", "bin")], on="scaffold"]
  a[degree == 1 & length <= 1e4]$scaffold->add
  if(length(add) > 0){
   ex <- c(ex, add)
   make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
   out$m -> m
  }

  # remove short tips/bulges of rank 1 
  m[rank == 1 & length <= 1e4]$scaffold -> add
  ex <- c(ex, add)
  if(length(add) > 0){
   ex <- c(ex, add)
   make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
   out$m -> m
  }

  # resolve length-one-bifurcations at the ends of paths
  m[rank > 0][bin == 2 | super_nbin - 1 == bin ][, .(super, super_nbin, type = bin == 2, scaffold, length, bin0=bin)]->x
  unique(rbind(
  m[x[type == T, .(super, bin0, bin=1)], on=c("super", "bin")],
  m[x[type == F, .(super, bin0, bin=super_nbin)], on=c("super", "bin")]
  ))->a
  a[, .(super, bin0, scaffold2=scaffold, length2=length)][x, on=c("super", "bin0")][, ex := ifelse(length >= length2, scaffold2, scaffold)]$ex -> add

  if(length(add) > 0){
   ex <- c(ex, add)
   make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
   out$membership -> m
  }

  # remove short tips/bulges of rank 1 
  m[rank == 1 & length <= 1e4]$scaffold -> add
  if(length(add) > 0){
   ex <- c(ex, add)
   make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
   out$m -> m
  }
  
  # remove tips of rank 1 
  links[!scaffold2 %in% ex][, .(degree=.N), key=.(scaffold=scaffold1)]->b
  b[m[rank == 1], on="scaffold"][degree == 1]$scaffold -> add
  if(length(add) > 0){
   ex <- c(ex, add)
   make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
   out$m -> m
  }

  # remove remaining nodes of rank > 0
  m[rank > 0]$scaffold -> add
  if(length(add) > 0){
   ex <- c(ex, add)
   make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
   out$m -> m
  }

  if(popseq_dist > 0 & unanchored == T){
   if(verbose){
    cat("Including unanchored scaffolds.\n")
   }
   # use unanchored scaffolds to link super-scaffolds
   ww2[is.na(popseq_chr1), .(scaffold_link=scaffold1, link_length=length1, scaffold1=scaffold2)]->x
   ww2[is.na(popseq_chr1), .(scaffold_link=scaffold1, scaffold2=scaffold2)]->y
   x[y, on="scaffold_link", allow.cartesian=T][scaffold1 != scaffold2]->xy

   m[, .(scaffold1=scaffold, super1=super, chr1=chr, cM1=cM, size1=super_nbin, d1 = pmin(bin - 1, super_nbin - bin))][xy, on="scaffold1"]->xy
   m[, .(scaffold2=scaffold, super2=super, chr2=chr, cM2=cM, size2=super_nbin, d2 = pmin(bin - 1, super_nbin - bin))][xy, on="scaffold2"]->xy
   xy[super2 != super1 & d1 == 0 & d2 == 0 & size1 > 1 & size2 > 1 & chr1 == chr2]->xy
   xy[scaffold1 < scaffold2, .(nscl=.N), scaffold_link][xy, on="scaffold_link"]->xy
   xy[nscl == 1] -> xy
   xy[super1 < super2][, c("n", "g"):=list(.N, .GRP), by=.(super1, super2)][order(-link_length)][!duplicated(g)]->zz

   sel <- zz[, .(scaffold1=c(scaffold_link, scaffold_link, scaffold1, scaffold2),
	  scaffold2=c(scaffold1, scaffold2, scaffold_link, scaffold_link))]
   rbind(links, ww2[sel, on=c("scaffold1", "scaffold2")])->links2

   make_super_scaffolds(links=links2, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
   out$m -> m

   #resolve branches
   m[rank > 0][bin == 2 | super_nbin - 1 == bin ][, .(super, super_nbin, type = bin == 2, scaffold, length, bin0=bin)]->x
   unique(rbind(
   m[x[type == T, .(super, bin0, bin=1)], on=c("super", "bin")],
   m[x[type == F, .(super, bin0, bin=super_nbin)], on=c("super", "bin")]
   ))->a
   a[, .(super, bin0, scaffold2=scaffold, length2=length)][x, on=c("super", "bin0")][, ex := ifelse(length >= length2, scaffold2, scaffold)]$ex -> add

   if(length(add) > 0){
    ex <- c(ex, add)
    make_super_scaffolds(links=links2, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
    out$membership -> m
   }

   # remove remaining nodes of rank > 0
   m[rank > 0]$scaffold -> add
   if(length(add) > 0){
    ex <- c(ex, add)
    make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
   }
  }

  out$m -> m
  out$info -> res
 } 

 if(verbose){
  cat("Orienting scaffolds.\n")
 }
 #orient scaffolds
 m[super_nbin > 1, .(scaffold1=scaffold, bin1=bin, super1=super)][link_pos, on="scaffold1", nomatch=0]->a
 m[super_nbin > 1, .(scaffold2=scaffold, bin2=bin, super2=super)][a, on="scaffold2", nomatch=0]->a
 a[super1 == super2 & bin1 != bin2]->a
 a[, d := abs(bin2 - bin1)]->a
 a[d <= max_dist_orientation]->a
 a[, .(nxt=mean(pos1[bin2 > bin1]), prv=mean(pos1[bin2 < bin1])), key=.(scaffold=scaffold1)]->aa
 info[, .(scaffold, length)][aa, on="scaffold"]->aa
 aa[!is.nan(prv) & !is.nan(nxt), orientation := ifelse(prv <= nxt, 1, -1)]
 aa[is.nan(prv), orientation := ifelse(length - nxt <= nxt, 1, -1)]
 aa[is.nan(nxt), orientation := ifelse(length - prv <= prv, -1, 1)]

 aa[, .(orientation, scaffold)][m, on="scaffold"]->m
 m[, oriented := T]
 m[is.na(orientation), oriented := F]
 m[is.na(orientation), orientation := 1]
 setorder(m, super, bin, rank)

 m[, super_pos := 1 + cumsum(c(0, length[-.N])), by=super]

 if(verbose){
  cat("Anchoring super-scaffolds.\n")
 }
 # assign super scaffolds to genetic positions
 m[!is.na(chr), .(nchr=.N), key=.(chr, super)][, pchr := nchr/sum(nchr), by=super][order(-nchr)][!duplicated(super)]->y
 m[!is.na(cM)][y, on=c("chr", "super")]->yy
 y[res, on="super"]->res
 yy[, .(cM=mean(cM), min_cM=min(cM), max_cM=max(cM)), key=super][res, on='super']->res
 setorder(res, -length)

 make_agp(m, gap=100)->a

 list(membership=m, info=res, 
      agp=a$agp, agp_bed=a$agp_bed)
}

read_morex_aln <- function(assembly, super, paf, agp, fai0=NULL, minq=20){
 paf["P", on="type"]->z
 z[mapq >= minq]->z
 z[, seq := sub("-[0-9]+$", "", query)]
 if(!is.null(fai0)){
  fai0[, .(reference=scaffold0, scaffold)][z, on="reference"]->z
 } else {
  z[, scaffold := reference]
 }
 agp[z, on="seq"]->z
 z[, .(orig_scaffold=scaffold, orig_scaffold_start=reference_start, orig_scaffold_end=reference_end,
       orientation, agp_chr=chr, 
       agp_start = agp_start - 1 + query_start,
       agp_end = agp_start - 1 + query_end, seq, bac, hic_bin, cluster, mapq, alnlen)]->zz
 assembly$info[, .(scaffold, scaffold_length=length, orig_scaffold, orig_scaffold_start=orig_start, orig_start)]->b
 b[zz, on=c("orig_scaffold", "orig_scaffold_start"), roll=T]->zz
 zz[, scaffold_start :=  orig_scaffold_start - orig_start + 1]
 zz[, scaffold_end := orig_scaffold_end - orig_start + 1]
 super$membership[, .(scaffold, super, super_pos, super_orientation=orientation)][zz, on="scaffold"]->zz
 zz[super_orientation == 1, super_start := super_pos - 1 + scaffold_start]
 zz[super_orientation == -1, super_end := super_pos + (scaffold_length - scaffold_start)]
 zz[super_orientation == 1, super_end := super_pos - 1 + scaffold_end]
 zz[super_orientation == -1, super_start := super_pos + (scaffold_length - scaffold_end)]
 zz[super_orientation == -1, orientation := -1L * orientation]
 super$info[, .(super, super_chr=chr)][zz, on="super"]->zz
 zz[]
}

plot_morex_aln<-function(data, morex_aln, file, ncores=1, minlen=5e3){
 copy(data)[, idx := 1:.N] -> data

 plotfu<-function(data, aln, mem=NULL, minlen){
  cat(paste0(data$idx, " ", data$super, "\n"))
  aln[data$super, on="super"][agp_chr == super_chr][alnlen >= minlen] -> xx
  if(nrow(xx) > 0){
   quantile(xx$agp_start, 0:50/50)[c(2, 50)]/1e6 -> ylim
   xlim <- c(0, data$length)/1e6
   xx[, plot(0, xlim=xlim, type='n', ylim=ylim, las=1, xlab="10X super-scaffold position (Mb)",
	     ylab="Morex AGP position (Mb)", bty='l')]
   data[, title(main=paste0(sub("super", "super_scaffold", super), ", ", chr, "H, ", round(length/1e6, 1), " Mb"))]
   if(!is.null(mem)){
    abline(v=mem[data$super, on="super"]$super_pos/1e6, col="gray")
   }
   xx[, idx := 1:.N]
   xx[orientation == 1, lines(c(super_start/1e6, super_end/1e6), c(agp_start/1e6, agp_end/1e6)), by=idx]
   xx[orientation == -1, lines(c(super_start/1e6, super_end/1e6), c(agp_end/1e6, agp_start/1e6)), by=idx]
  } else {
   plot(0, type='n', axes=F, xlab="", ylab="")
   data[, title(main=paste0(sub("super", "super_scaffold", super), ", ", chr, "H, ", round(length/1e6, 1), " Mb",
			   "\nNo alignments to Morex pseudomolecules."))]
  }
 }
 parallel_plot(data=data, group="super", cores=ncores, aln=morex_aln, minlen=minlen,
		file=file, height=700, width=700, res=150, plot_function=plotfu)
}

# Read alignment of 10X super-scaffolds to another scaffolds. First used for comparison to NRGene assemblies, hence the name.
read_nrgene_aln <- function(assembly, super, paf, fai0=NULL, minq=20){
 paf["P", on="type"]->z
 z[mapq >= minq]->z
 if(!is.null(fai0)){
  fai0[, .(query=scaffold0, scaffold)][z, on="query"]->z
 } else {
  z[, scaffold := query]
 }
 z[, .(orig_scaffold=scaffold, orig_scaffold_start=query_start, orig_scaffold_end=query_end,
       nrgene=reference, nrgene_length=reference_length, orientation,
       nrgene_start=reference_start,  nrgene_end=reference_end, mapq, alnlen)]->zz
 assembly$info[, .(scaffold, scaffold_length=length, orig_scaffold, orig_scaffold_start=orig_start, orig_start)]->b
 b[zz, on=c("orig_scaffold", "orig_scaffold_start"), roll=T]->zz
 zz[, scaffold_start :=  orig_scaffold_start - orig_start + 1]
 zz[, scaffold_end := orig_scaffold_end - orig_start + 1]
 super$membership[, .(scaffold, super, super_pos, super_orientation=orientation)][zz, on="scaffold"]->zz
 zz[super_orientation == 1, super_start := super_pos - 1 + scaffold_start]
 zz[super_orientation == -1, super_end := super_pos + (scaffold_length - scaffold_start)]
 zz[super_orientation == 1, super_end := super_pos - 1 + scaffold_end]
 zz[super_orientation == -1, super_start := super_pos + (scaffold_length - scaffold_end)]
 zz[super_orientation == -1, orientation := -1L * orientation]
 super$info[, .(super, super_length=length, super_chr=chr)][zz, on="super"]->zz
 zz[]
}

# Plot alignment of 10X super-scaffolds to another assemblies. First used for comparison to NRGene assemblies, hence the name.
plot_nrgene_aln<-function(data, nrgene_aln, mem=NULL, file, ncores=1, minlen=2e4){
 copy(data)[, plot_idx := 1:.N] -> data

 plotfu<-function(data, aln, mem=NULL, minlen){
  cat(paste0(data$plot_idx, " ", data$super, "\n"))
  aln[data[, .(super, nrgene)], on=c("nrgene", "super")][alnlen >= minlen] -> xx
  if(nrow(xx) > 0){
   xlim <- c(0,0)
   ylim <- c(0,0)
   quantile(xx$super_start, 0:50/50)[2]/1e6 -> xlim[1]
   quantile(xx$super_end, 0:50/50)[50]/1e6 -> xlim[2]
   quantile(xx$nrgene_start, 0:50/50)[2]/1e6 -> ylim[1]
   quantile(xx$nrgene_end, 0:50/50)[50]/1e6 -> ylim[2]
   xx[, plot(0, xlim=xlim, type='n', ylim=ylim, las=1, xlab="10X super-scaffold position (Mb)",
	     ylab="NRGene scaffold position (Mb)", bty='l')]
   data[, title(main=paste0(sub("super", "super_scaffold", super), ", ", round(super_length/1e6, 1), " Mb\n",
		            sub("scaffold", "NRGene scaffold", nrgene), ", ", round(nrgene_length/1e6, 1), " Mb"))]
   if(!is.null(mem)){
    abline(v=mem[data$super, on="super"]$super_pos/1e6, col="gray")
   }
   xx[, idx := 1:.N]
   xx[orientation == 1, lines(c(super_start/1e6, super_end/1e6), c(nrgene_start/1e6, nrgene_end/1e6)), by=idx]
   xx[orientation == -1, lines(c(super_start/1e6, super_end/1e6), c(nrgene_end/1e6, nrgene_start/1e6)), by=idx]
  } else {
   plot(0, type='n', axes=F, xlab="", ylab="")
   data[, title(main=paste0(sub("super", "super_scaffold", super), ", ", round(super_length/1e6, 1), " Mb\n",
		            sub("scaffold", "NRGene scaffold", nrgene), ", ", round(nrgene_length/1e6, 1), " Mb"))]
  }
 }
 parallel_plot(data=data, group="idx", cores=ncores, aln=nrgene_aln, minlen=minlen, mem=mem,
		file=file, height=700, width=700, res=150, plot_function=plotfu)
}

# Iteratively break scaffolds using 10X physical coverage. Proceed until no more breakpoints are found.
break_10x<-function(assembly, prefix="scaffold_corrected", species="wheat", ratio=-3, interval=5e4, minNbin=20, dist=2e3, slop=1e3, 
		    intermediate=F, ncores=1, maxcycle=Inf){

 if(dist <= 2 * slop){
  dist <- 2 * slop + 1
  cat(paste0("Setting dist to ", 2 * slop + 1, "\n"))
 }
 find_10x_breaks(assembly=assembly, interval=interval, minNbin=minNbin, dist=dist, ratio=ratio) -> breaks
 i <- 1
 lbreaks <- list()
 lbreaks[[1]] <- breaks
 if(intermediate){
  assemblies <- list()
  assemblies[[1]] <- assembly
 }
 while(nrow(breaks) > 0 & i <= maxcycle){
  cat(paste0("Cycle ", i, ": ", breaks[, .N], " break points detected\n"))
  i <- i + 1
  break_scaffolds(breaks=breaks, assembly=assembly, species=species, prefix=paste0(prefix, "_"), slop=slop, cores=ncores) -> assembly
  if(intermediate){
   assemblies[[i]] <- assembly
  }
  find_10x_breaks(assembly=assembly, interval=interval, minNbin=minNbin, dist=dist, ratio=ratio) -> breaks
  if(is.null(breaks)){
   break
  }
  lbreaks[[i]] <- breaks
 }
 if(!intermediate){
  list(assembly=assembly, breaks=lbreaks)
 } else {
  list(assemblies=assemblies, breaks=lbreaks)
 }
}

# Initialize an assembly object for super-scaffolds created by scaffold_10x(). Lift positional information from scaffolds to super-scaffolds.
init_10x_assembly<-function(assembly, map_10x, molecules=F){
 super <- map_10x

 copy(assembly$cssaln)->z
 z[, orig_scaffold := NULL]
 z[, orig_scaffold_length := NULL]
 z[, orig_pos := NULL]
 z[, scaffold_length := NULL]
 super$agp[, .(scaffold, super, super_start, super_end, orientation)][z, on="scaffold"]->z
 z[orientation == 1, pos := super_start - 1 + pos]
 z[orientation == -1, pos := super_end - pos + 1]
 z[, scaffold := NULL]
 setnames(z, "super", "scaffold")
 z[, c("super_start", "super_end", "orientation") := list(NULL, NULL, NULL)]
 super$info[, .(scaffold=super, scaffold_length=length)][z, on="scaffold"]->z
 z->s_cssaln

 if(molecules){
  copy(assembly$molecules)->z
  z[, c("orig_scaffold", "orig_start") := list(NULL, NULL)]
  super$agp[, .(scaffold, super, super_start, super_end, orientation)][z, on="scaffold"]->z
  z[orientation == 1, start := super_start - 1 + start]
  z[orientation == 1, end := super_start - 1 + end]
  z[orientation == -1, start := super_end - end + 1]
  z[orientation == -1, end := super_end - start + 1]
  z[, scaffold := NULL]
  z[, c("super_start", "super_end", "orientation") := list(NULL, NULL, NULL)]
  setnames(z, "super", "scaffold")
  z -> s_molecules
 } else {
  s_molecules <- NULL
 }

 if(!is.null(assembly$fpairs) && nrow(assembly$fpairs) > 0){
  copy(assembly$fpairs)->z
  z[, c("orig_scaffold1", "orig_pos1") := list(NULL, NULL)]
  z[, c("orig_scaffold2", "orig_pos2") := list(NULL, NULL)]
  z[, c("chr1", "chr2") := list(NULL, NULL)]
  super$agp[, .(scaffold1=scaffold, super, super_start, super_end, orientation)][z, on="scaffold1"]->z
  z[orientation == 1, pos1 := super_start - 1 + pos1]
  z[orientation == -1, pos1 := super_end - pos1 + 1]
  z[, scaffold1 := NULL]
  setnames(z, "super", "scaffold1")
  z[, c("super_start", "super_end", "orientation") := list(NULL, NULL, NULL)]
  super$agp[, .(scaffold2=scaffold, super, super_start, super_end, orientation)][z, on="scaffold2"]->z
  z[orientation == 1, pos2 := super_start - 1 + pos2]
  z[orientation == -1, pos2 := super_end - pos2 + 1]
  z[, scaffold2 := NULL]
  setnames(z, "super", "scaffold2")
  z[, c("super_start", "super_end", "orientation") := list(NULL, NULL, NULL)]
 } else {
  z <- NULL
 }

 init_assembly(fai=super$agp_bed[, .(length=sum(bed_end - bed_start)), key=.(scaffold=super)], cssaln=s_cssaln, molecules=s_molecules, fpairs=z) 
}

# Convert original scaffolds positions to AGP positions
orig_scaffold_to_agp<-function(assembly, map_10x=NULL, hic_map){
 if(!is.null(map_10x)){
  super <- map_10x
  assembly$info[, .(orig_scaffold, orig_scaffold_start=orig_start, 
			      orig_scaffold_end=orig_end, scaffold)]->zz
  super$agp[, .(scaffold, scaffold_length=length, super, super_start, super_end,
		       super_orientation=orientation)][zz, on="scaffold"]->zz
  hic_map$agp[scaffold != "gap", .(agp_chr, agp_start0 = agp_start, agp_end0 = agp_end, super=scaffold, agp_orientation = ifelse(is.na(orientation), 1, orientation))]->a
  a[zz, on="super"]->zz
  zz[, orientation := agp_orientation * super_orientation]
  zz[agp_orientation == 1, agp_start := agp_start0 - 1 + super_start] 
  zz[agp_orientation == 1, agp_end := agp_start0 - 1 + super_end] 
  zz[agp_orientation == -1, agp_start := agp_end0 + 1 - super_end] 
  zz[agp_orientation == -1, agp_end := agp_end0 + 1 - super_start] 
  zz[, agp_start0 := NULL]
  zz[, agp_end0 := NULL]
  zz[]
 } else {
  assembly$info[, .(orig_scaffold, orig_scaffold_start=orig_start, 
			       orig_scaffold_end=orig_end, scaffold)]->zz
  hic_map$agp[scaffold != "gap", .(agp_chr, agp_start, agp_end, scaffold, agp_orientation = ifelse(is.na(orientation), 1, orientation))]->a
  a[zz, on="scaffold"]->zz
  zz[, orientation := agp_orientation]
  zz[]
 }
}

# Lift alignment coordinate: query to AGP; reference (=IBSC2017 assembly) BACs to pseudonmolecule
lift_morex_aln_agp<-function(paf, morex_agp, scaffold_to_agp){
 paf[, .(seq=sub("-[0-9]+$", "", query), seq_start=query_start, seq_end=query_end,
	       orig_scaffold=reference, orig_scaffold_start = reference_start, orig_scaffold_end = reference_end,
	       aln_orientation = orientation, mapq, alnlen)]->p
 morex_agp[, .(seq, morex_chr=agp_chr, morex_agp_start0=agp_start)][p, on="seq"]->p
 p[, morex_agp_start := morex_agp_start0 - 1 + seq_start]
 p[, morex_agp_end := morex_agp_start0 - 1 + seq_end]
 scaffold_to_agp[, .(orig_scaffold, orig_pos = orig_scaffold_start, orig_scaffold_start, agp_chr, agp_start0=agp_start, agp_end0=agp_end, orientation)][p, on=c("orig_scaffold", "orig_scaffold_start"), roll=T]->p
 p[, aln_orientation := orientation * aln_orientation]
 p[orientation == 1, agp_start := agp_start0 + (orig_scaffold_start - orig_pos)]
 p[orientation == 1, agp_end := agp_start0 + (orig_scaffold_end - orig_pos)]
 p[orientation == -1, agp_start := agp_end0 - (orig_scaffold_end - orig_pos)]
 p[orientation == -1, agp_end := agp_end0 - (orig_scaffold_start - orig_pos)]
 p[, agp_start0 := NULL]
 p[, agp_end0 := NULL]
 p[, orig_pos := NULL]
 p[]
}

# Wrapper function for Hi-C mapping
make_hic_info<-function(cluster_info, super_global, chrs){
 s<-super_global$super_info
 s[!duplicated(s$chr),]->s
 s[chr %in% chrs]->s

 super_global$membership[, .(cluster, super, bin, rank, backbone)]->tmp
 tmp[super %in% s$super]->tmp
 tmp[, super := NULL]
 setnames(tmp, c("cluster", "hic_bin", "hic_rank", "hic_backbone"))
 tmp[cluster_info, on="cluster"]->cluster_info
 cluster_info[order(chr, hic_bin, hic_rank, cluster)]
}

# Insert nodes not in MST backbone to where they fit best
insert_node<-function(df, el){
 while({
  setkey(el[cluster2 %in% df$cluster & cluster1 %in% unique(el[!cluster1 %in% df$cluster & cluster2 %in% df$cluster]$cluster1)], "cluster2")[setkey(setnames(df[, .(cluster, bin)], c("cluster2", "bin")), "cluster2"), allow.cartesian=T][!is.na(cluster1)][
   order(cluster1, bin)][, dist:={if(.N == 1) {as.integer(NA)} else { as.integer(c(bin[2:.N],NA)-bin)}}, by=cluster1]->y
   y[order(cluster1, bin)]->y
   length(which(y$dist == 1)->idx) > 0
  }) {

   setkeyv(setnames(el[, .(cluster1, cluster2, weight)], c("path_node1", "path_node2", "old_path")), c("path_node1", "path_node2"))[
    setkeyv(data.table(cluster=y[idx, cluster1], path_node1=y[idx, cluster2], path_node2=y[idx+1, cluster2], 
	      weight1=y[idx, weight], weight2=y[idx+1, weight], bin=y[idx, bin]), c("path_node1", "path_node2"))]->z
    head(z[ ,diff:= weight1 + weight2 - old_path][order(diff)],1)->z
   
   m=z$bin
   data.table(rbind(df[1:m], data.table(cluster=z$cluster, bin=m+1), df[(m+1):nrow(df),][, bin:=bin+1][]))->df
 }
 df
}

# Traveling salesman heuristic for Hi-C mapping construction
node_relocation<-function(df, el, maxiter=100, verbose=T){
 i<-0
 if(nrow(df) > 2){
  while({
   i<-i+1
   df[order(bin)]->df
   x<-data.table(old_node1=df[1:(nrow(df)-2)]$cluster, cluster=df[2:(nrow(df)-1)]$cluster, old_node2=df[3:nrow(df)]$cluster)
   setkeyv(el[, .(cluster1, cluster2, weight)], c("cluster1", "cluster2"))->ee
   setnames(setnames(ee, c("cluster1", "cluster2"), c("old_node1", "old_node2"))[setkeyv(x, c("old_node1", "old_node2"))], "weight", "new_edge3")->x
   setnames(setnames(ee, c("old_node1", "old_node2"), c("cluster", "old_node1"))[setkeyv(x, c("cluster", "old_node1"))], "weight", "old_edge1")->x
   setnames(setnames(ee, "old_node1", "old_node2")[setkeyv(x, c("cluster", "old_node2"))], "weight", "old_edge2")->x
   x[!is.na(new_edge3) ]->x

   setkey(el[, .(cluster1, cluster2)], "cluster2")[setkey(copy(df), "cluster")][order(cluster1, bin)]->t
   which(t[, dist:={if(.N == 1) {as.integer(NA)} else { as.integer(c(bin[2:.N],NA)-bin)}}, by=cluster1]$dist == 1)->idx

   data.table(cluster=t$cluster1[idx], new_node1=t$cluster2[idx], new_node2=t$cluster2[idx+1])->t
   setkeyv(el[, .(cluster1, cluster2, weight)], c("cluster1", "cluster2"))->ee
   setnames(setnames(ee, c("cluster1", "cluster2"), c("cluster", "new_node1"))[setkeyv(t, c("cluster", "new_node1"))], "weight", "new_edge1")->t
   setnames(setnames(ee, "new_node1", "new_node2")[setkeyv(t, c("cluster", "new_node2"))], "weight", "new_edge2")->t
   setnames(setnames(ee, c("cluster", "new_node2"), c("new_node1", "new_node2"))[setkeyv(t, c("new_node1", "new_node2"))], "weight", "old_edge3")->t
   setkey(x, "cluster")[setkey(t, "cluster")][!is.na(old_node1)]->m
   m[ ,diff := new_edge1+new_edge2+new_edge3 - old_edge1 - old_edge2 - old_edge3]
   nrow(m<-m[diff < 0][order(diff)]) > 0 & i <= maxiter
   }) {
   ne<-head(data.frame(m), 1)
   idx<-data.frame(cluster=df$cluster, idx=1:nrow(df))
   invisible(lapply(c("cluster", "new_node1", "new_node2", "old_node1", "old_node2"), function(i) {
    merge(ne, by.x=i, by.y="cluster", idx)->>ne
    colnames(ne)[which(colnames(ne) == "idx")]<<-paste(sep="_", "idx", i)
   }))
   df[with(ne, {
    min_new <- min(idx_new_node1, idx_new_node2)
    max_new <- max(idx_new_node1, idx_new_node2)
    min_old <- min(idx_old_node1, idx_old_node2)
    max_old <- max(idx_old_node1, idx_old_node2)
    if(min_old < min_new) {
     c(1:min_old, max_old:min_new, idx_cluster, max_new:nrow(df))
    } else {
     c(1:min_new, idx_cluster, max_new:min_old, max_old:nrow(df))
    }
   }),]->df
   df$bin<-1:nrow(df)
   df<-data.table(df)
  }
  if(verbose){
   cat(paste0("Node relocation steps: ", i-1, "\n"))
  }
 }
 df
}

# Traveling salesman heuristic for Hi-C mapping construction
kopt2<-function(df, el){
 setkeyv(el, c("cluster1", "cluster2"))->el
 while({
  df[, .(cluster1=cluster[1:(nrow(df)-1)], cluster2=cluster[2:nrow(df)])]->m
  el[setkeyv(m, c("cluster1", "cluster2"))]->m
  m[ ,.(cluster1, cluster2, weight)]->m
  setnames(m, "weight", "weight12")->m
  m[, dummy:=1]
  setkey(setnames(copy(df), c("cluster", "bin"), c("cluster1", "bin1")), "cluster1")[setkey(m, "cluster1")]->m
  setkey(setnames(copy(df), c("cluster", "bin"), c("cluster2", "bin2")), "cluster2")[setkey(m, "cluster2")]->m
  copy(m)->n
  setnames(n, c("cluster1", "cluster2", "bin1", "bin2", "weight12"), c("cluster3", "cluster4", "bin3", "bin4", "weight34"))
  setkey(m, "dummy")[setkey(n, "dummy"), allow.cartesian=T]->mn
  mn[, dummy:=NULL]
  mn[bin1 < bin3]->mn
  o<-el[, .(cluster1, cluster2, weight)]
  setkeyv(setnames(copy(o), c("cluster1", "cluster3", "weight13")), c("cluster1", "cluster3"))[setkeyv(mn, c("cluster1", "cluster3"))]->mn
  setkeyv(setnames(copy(o), c("cluster2", "cluster4", "weight24")), c("cluster2", "cluster4"))[setkeyv(mn, c("cluster2", "cluster4"))]->mn
  mn[, old:=weight12+weight34]
  mn[, new:=weight13+weight24]
  mn[, diff:=old-new]
  mn[order(-diff)]->mn
  nrow(mn[diff > 0]) > 0 
 }) {
  head(mn, 1)->x
  bin1<-df[cluster == x$cluster1]$bin
  bin2<-df[cluster == x$cluster2]$bin
  bin3<-df[cluster == x$cluster3]$bin
  bin4<-df[cluster == x$cluster4]$bin
  df[c(1:bin1, bin3:bin2, bin4:nrow(df))]->df
  df[, bin:=1:nrow(df)]->df
 }
 df
}

# Set up Hi-C graph structure, determine ends of chromosome from genetic map, run Hi-C ordering for each chromosome
make_super<-function(hl, cluster_info, prefix="super", cores=1, paths=T, path_max=0, known_ends=F, 
		     maxiter=100, verbose=T){

 hl[cluster1 %in% cluster_info[excluded == F]$cluster & cluster2 %in% cluster_info[excluded == F]$cluster]->hl

 hl[cluster1 < cluster2]->e
 graph.edgelist(as.matrix(e[, .(cluster1, cluster2)]), directed=F)->g
 E(g)$weight<-e$weight

 data.table(cluster=V(g)$name, super=paste(prefix, sep="_", clusters(g)$membership))->mem
 cluster_info[mem, on="cluster"]->mem

 mem[, .(super_size=.N, length=.N, chr=unique(na.omit(chr))[1], cM=mean(na.omit(cM))), keyby=super]->info
 mem[, .(cluster1=cluster, super)][hl, on="cluster1"]->e

 list(super_info=info, membership=mem, graph=g, edges=e)->s

 if(paths){
  if(path_max > 0){
   idx<-head(s$super_info[order(-length)], n=path_max)$super
  } else {
   idx<-s$super_info$super
  }
  rbindlist(mclapply(mc.cores=cores, idx, function(i) {
   start <- end <- NULL
   # Take terminal nodes from genetic map
   if(known_ends){
    s$mem[super == i & !is.na(cM)][order(cM)]$cluster->x
    start=x[1]
    end=tail(x,1)
   } 
   make_super_path(s, idx=i, start=start, end=end, maxiter=maxiter, verbose=verbose)->x
   if(verbose){
    cat(paste0("Chromosome ", head(s$mem[super == i]$chr, n=1), " finished.\n"))
   }
   x
  }))[s$membership, on="cluster"]->s$membership
 }

 s
}

# Order scaffolds using Hi-C links for one chromosome
make_super_path<-function(super, idx=NULL, start=NULL, end=NULL, maxiter=100, verbose=T){

 # Get backbone from minimum spanning tree (MST)
 submem<-super$mem[super == idx]
 super$edges[super == idx, .(cluster1, cluster2, weight)]->el
 minimum.spanning.tree(induced.subgraph(super$graph, submem$cluster))->mst

 E(mst)$weight <- 1

 if(is.null(start) | is.null(end)){
  V(mst)[get.diameter(mst)]$name->dia
 } else {
  V(mst)[get.shortest.paths(mst, from=start, to=end)$vpath[[1]]]$name->dia
 }

 data.table(cluster=dia, bin=1:length(dia))->df

 df<-insert_node(df, el)
 # Traveling salesman heuristics
 df<-kopt2(df, el)
 df<-node_relocation(df, el, maxiter=maxiter, verbose=verbose)

 data.frame(df)->df
 data.frame(el)->el
 data.frame(cluster=df$cluster, rank = 0)->ranks
 r=0

 while(length(n<-unique(subset(el, !cluster1 %in% df$cluster & cluster2 %in% df$cluster)$cluster1)) > 0) {
  r = r+1
  subset(el, cluster2 %in% df$cluster & cluster1 %in% n)->tmp
  tmp[!duplicated(tmp$cluster1),]->tmp
  rbind(ranks, data.frame(cluster=tmp$cluster1, rank=r))->ranks
  merge(tmp, df[c("cluster", "bin")], by.x="cluster2", by.y="cluster")->x
  rbind(df, data.frame(cluster=x$cluster1, bin=x$bin))->df
 }

 merge(df, submem)->df
 df$bin<-as.integer(df$bin)
 flip<-with(df, suppressWarnings(cor(bin, cM))) < 0
 if((!is.na(flip)) & flip) {
  with(df, max(bin) - bin + 1)->df$bin
 }
 ranks$rank<-as.numeric(ranks$rank)
 merge(df, ranks)->df
 df$backbone <- df$cluster %in% dia

 data.table(df)[, .(cluster, bin, rank, backbone)]
}

# Wrapper function for Hi-C mapping
make_hic_map<-function(hic_info, links, ncores=1, maxiter=100, known_ends=T){

 copy(links)->hl
 copy(hic_info)->info

 setnames(hl, c("scaffold1", "scaffold2"), c("cluster1", "cluster2"))
 setnames(info, "scaffold", "cluster")
 setkey(info, "cluster")
 chrs <- info[!is.na(chr), unique(chr)]

 make_hic_info(info, 
  super_global<-make_super(hl, cluster_info=info, cores=ncores, maxiter=maxiter,
			   known_ends=known_ends, path_max=length(chrs)), chrs=chrs)->res
 res[order(chr, hic_bin, cluster)][, .(scaffold=cluster, chr, cM, hic_bin, hic_backbone, hic_rank)][!is.na(hic_bin)] -> res
 setnames(res, "hic_bin", "hic_bin0")
 res[, hic_bin := 1:.N, by=chr]
 res[]
}

# Flip AGP orientation of specified scaffolds
correct_inversions<-function(hic_map, scaffolds, species){
 copy(hic_map$hic_map)->z
 z[scaffold %in% scaffolds, consensus_orientation := ifelse(is.na(consensus_orientation), -1, -1 * consensus_orientation)]
 make_agp(z, gap_size=hic_map$gap_size, species=species)->a

 list()->new
 new$hic_map <- z
 new$agp <- a$agp
 new$gap_size <- copy(hic_map$gap_size)
 new$chrlen <- copy(hic_map$chrlen)
 new$binsize <- copy(hic_map$binsize)
 new$hic_map_bin <- copy(hic_map$hic_map_bin)
 new$max_cM_dist <- copy(hic_map$max_cM_dist)
 new$min_nfrag_bin <- copy(hic_map$min_nfrag_bin)
 new$agp_bed <- a$agp_bed
 new$corrected_inversions <- scaffolds
 new
}

# Call Hi-C function for ordering, orient by ordering parts of scaffolds, create AGP table
hic_map<-function(info, assembly, frags, species, ncores=1, min_nfrag_scaffold=50, max_cM_dist = 20, 
		  binsize=5e5, min_nfrag_bin=30, gap_size=100, maxiter=100, orient=T, agp_only=F,
		  map=NULL, known_ends=T, orient_old=F, min_binsize=1e5, min_nbin=5){
 if(!agp_only){
  copy(info)->hic_info
  hic_info[, excluded := nfrag < min_nfrag_scaffold]

  print("running Hi-C mapping")
  
  assembly$fpairs[scaffold1 != scaffold2, .(nlinks=.N), key=.(scaffold1, scaffold2)]->hl
  hic_info[, .(scaffold1=scaffold, chr1=chr, cM1=cM)][hl, nomatch=0, on="scaffold1"]->hl
  hic_info[, .(scaffold2=scaffold, chr2=chr, cM2=cM)][hl, nomatch=0, on="scaffold2"]->hl
  hl[chr1 == chr2]->hl
  hl<-hl[abs(cM1-cM2) <= max_cM_dist | is.na(cM1) | is.na(cM2)]
  hl[, weight:=-log10(nlinks)]

  cat("Scaffold map construction started.\n")
  make_hic_map(hic_info=hic_info, links=hl, ncores=ncores, known_ends=known_ends)->hic_map
  cat("Scaffold map construction finished.\n")

  if(orient){
   if(orient_old){
    options(scipen = 1000)
    frags[, .(nfrag=.N), keyby=.(scaffold, pos = start %/% binsize * binsize)]->fragbin
    fragbin[, id := paste(sep=":", scaffold, pos)]
    fragbin<-hic_info[excluded == F, .(scaffold, chr, cM)][fragbin, on="scaffold", nomatch=0]

    assembly$fpairs[, .(nlinks=.N), keyby=.(scaffold1, pos1 = pos1 %/% binsize * binsize, scaffold2, pos2 = pos2 %/% binsize * binsize)]->binl
    binl[, id1 := paste(sep=":", scaffold1, pos1)]
    binl[, id2 := paste(sep=":", scaffold2, pos2)]
    binl[id1 != id2]->binl
    fragbin[, .(id1=id, chr1=chr, cM1=cM)][binl, on="id1"]->binl
    fragbin[, .(id2=id, chr2=chr, cM2=cM)][binl, on="id2"]->binl
    binl[, c("scaffold1", "scaffold2", "pos1", "pos2") := list(NULL, NULL, NULL, NULL)]
    setnames(binl, c("id1", "id2"), c("scaffold1", "scaffold2"))

    cat("Scaffold bin map construction started.\n")
    fragbin[, .(scaffold=id, nfrag, chr, cM)]->hic_info_bin
    hic_info_bin[, excluded:=nfrag < min_nfrag_bin]
    binl[chr1 == chr2 & (abs(cM1-cM2) <= max_cM_dist | is.na(cM1) | is.na(cM2))]->binl
    binl[, weight:=-log10(nlinks)]

    make_hic_map(hic_info=hic_info_bin, links=binl, ncores=ncores, maxiter=maxiter, known_ends=known_ends)->hic_map_bin
    cat("Scaffold bin map construction finished.\n")

    w<-hic_map_bin[!is.na(hic_bin), .(id=scaffold, scaffold=sub(":.*$", "", scaffold), pos=as.integer(sub("^.*:", "", scaffold)), chr, hic_bin)]
    w<-w[, .(gbin=mean(na.omit(hic_bin)),
	     hic_cor=as.numeric(suppressWarnings(cor(method='s', hic_bin, pos, use='p')))), keyby=scaffold][!is.na(hic_cor)]
    hic_map[!is.na(hic_bin) & scaffold %in% w$scaffold][order(chr, hic_bin)]->z0
    z0[,.(scaffold1=scaffold[1:(.N-2)], scaffold2=scaffold[2:(.N-1)], scaffold3=scaffold[3:(.N)]), by=chr]->z
    z0[, data.table(key="scaffold1", scaffold1=scaffold, hic_bin1=hic_bin)][setkey(z, "scaffold1")]->z
    z0[, data.table(key="scaffold2", scaffold2=scaffold, hic_bin2=hic_bin)][setkey(z, "scaffold2")]->z
    z0[, data.table(key="scaffold3", scaffold3=scaffold, hic_bin3=hic_bin)][setkey(z, "scaffold3")]->z
    w[, data.table(key="scaffold1", scaffold1=scaffold, gbin1=gbin)][setkey(z, "scaffold1")]->z
    w[, data.table(key="scaffold2", scaffold2=scaffold, gbin2=gbin)][setkey(z, "scaffold2")]->z
    w[, data.table(key="scaffold3", scaffold3=scaffold, gbin3=gbin)][setkey(z, "scaffold3")]->z
    z[, cc:= apply(z[, .(hic_bin1, hic_bin2, hic_bin3, gbin1, gbin2, gbin3)],1,function(x) {
		    suppressWarnings(cor(x[1:3], x[4:6]))
		    })]
    z[, data.table(key="scaffold", scaffold=scaffold2, cc=ifelse(cc > 0, 1, -1))]->ccor
    ccor[w]->m
    m[, hic_orientation:=ifelse(hic_cor > 0, 1 * cc, -1 * cc)]
    m[, .(scaffold, hic_cor, hic_invert=cc, hic_orientation)][hic_map, on="scaffold"]->hic_map_oriented

    setnames(hic_map_oriented, "chr", "consensus_chr")
    setnames(hic_map_oriented, "cM", "consensus_cM")
    hic_map_oriented[, consensus_orientation := hic_orientation]
   } else {
    options(scipen = 1000)
    assembly$info[, .(scaffold, binsize=pmax(min_binsize, length %/% min_nbin))][frags, on='scaffold']->f
    f[, .(nfrag=.N), keyby=.(scaffold, binsize, pos = start %/% binsize * binsize)]->fragbin
    fragbin[, id := paste(sep=":", scaffold, pos)]
    fragbin<- hic_map[, .(scaffold, chr, cM=hic_bin)][fragbin, on="scaffold", nomatch=0]

    unique(fragbin[, .(scaffold1=scaffold, binsize1=binsize)])[assembly$fpairs, on='scaffold1']->fp
    unique(fragbin[, .(scaffold2=scaffold, binsize2=binsize)])[fp, on='scaffold2']->fp
    fp[, .(nlinks=.N), keyby=.(scaffold1, pos1 = pos1 %/% binsize1 * binsize1, scaffold2, pos2 = pos2 %/% binsize2 * binsize2)]->binl
    binl[, id1 := paste(sep=":", scaffold1, pos1)]
    binl[, id2 := paste(sep=":", scaffold2, pos2)]
    binl[id1 != id2]->binl

    fragbin[, .(id1=id, chr1=chr, cM1=cM)][binl, on="id1"]->binl
    fragbin[, .(id2=id, chr2=chr, cM2=cM)][binl, on="id2"]->binl
    binl[, c("scaffold1", "scaffold2", "pos1", "pos2") := list(NULL, NULL, NULL, NULL)]
    setnames(binl, c("id1", "id2"), c("scaffold1", "scaffold2"))

    fragbin[, .(scaffold=id, nfrag, chr, cM)]->hic_info_bin
    hic_info_bin[, excluded:=nfrag < min_nfrag_bin]
    binl[chr1 == chr2 & (abs(cM1-cM2) <= 2 | is.na(cM1) | is.na(cM2))]->binl
    binl[, weight:=-log10(nlinks)]

    make_hic_map(hic_info=hic_info_bin, links=binl, ncores=ncores, maxiter=maxiter, known_ends=known_ends)->hic_map_bin

    w<-hic_map_bin[!is.na(hic_bin), .(id=scaffold, scaffold=sub(":.*$", "", scaffold), pos=as.integer(sub("^.*:", "", scaffold)), chr, hic_bin)]
    w<-w[, .(gbin=mean(na.omit(hic_bin)),
	     hic_cor=as.numeric(suppressWarnings(cor(method='s', hic_bin, pos, use='p')))), keyby=scaffold][!is.na(hic_cor)]
    hic_map[!is.na(hic_bin) & scaffold %in% w$scaffold][order(chr, hic_bin)]->z0
    z0[,.(scaffold1=scaffold[1:(.N-2)], scaffold2=scaffold[2:(.N-1)], scaffold3=scaffold[3:(.N)]), by=chr]->z
    z0[, data.table(key="scaffold1", scaffold1=scaffold, hic_bin1=hic_bin)][setkey(z, "scaffold1")]->z
    z0[, data.table(key="scaffold2", scaffold2=scaffold, hic_bin2=hic_bin)][setkey(z, "scaffold2")]->z
    z0[, data.table(key="scaffold3", scaffold3=scaffold, hic_bin3=hic_bin)][setkey(z, "scaffold3")]->z
    w[, data.table(key="scaffold1", scaffold1=scaffold, gbin1=gbin)][setkey(z, "scaffold1")]->z
    w[, data.table(key="scaffold2", scaffold2=scaffold, gbin2=gbin)][setkey(z, "scaffold2")]->z
    w[, data.table(key="scaffold3", scaffold3=scaffold, gbin3=gbin)][setkey(z, "scaffold3")]->z
    z[, cc:= apply(z[, .(hic_bin1, hic_bin2, hic_bin3, gbin1, gbin2, gbin3)],1,function(x) {
		    suppressWarnings(cor(x[1:3], x[4:6]))
		    })]
    z[, data.table(key="scaffold", scaffold=scaffold2, cc=ifelse(cc > 0, 1, -1))]->ccor
    ccor[w]->m
    m[, hic_orientation:=ifelse(hic_cor > 0, 1 * cc, -1 * cc)]
    m[, .(scaffold, hic_cor, hic_invert=cc, hic_orientation)][hic_map, on="scaffold"]->hic_map_oriented
    hic_map_oriented[is.na(hic_orientation), hic_orientation := ifelse(hic_cor > 0, 1, -1)]
    setnames(hic_map_oriented, "chr", "consensus_chr")
    setnames(hic_map_oriented, "cM", "consensus_cM")
    hic_map_oriented[, consensus_orientation := hic_orientation]
   }
  } else {
   hic_map_oriented<-copy(hic_map)
   setnames(hic_map_oriented, "chr", "consensus_chr")
   setnames(hic_map_oriented, "cM", "consensus_cM")
   hic_map_oriented[, consensus_orientation := as.numeric(NA)]
   hic_map_oriented[, hic_cor := as.numeric(NA)]
   hic_map_oriented[, hic_invert := as.numeric(NA)]
   hic_map_oriented[, hic_orientation := as.numeric(NA)]
   hic_map_bin <- NA
  }

  if("orientation" %in% names(hic_info)){
   hic_info[, .(scaffold, old_orientation=orientation)][hic_map_oriented, on="scaffold"]->hic_map_oriented
   hic_map_oriented[!is.na(old_orientation), consensus_orientation := old_orientation]
   hic_map_oriented[, old_orientation := NULL]
  }

  hic_map_oriented[assembly$info, on="scaffold"]->hic_map_oriented
 } else {
  hic_map_oriented <- map$hic_map
  hic_map_bin <- map$hic_map_bin
  min_nfrag_scaffold <- map$min_nfrag_scaffold
  binsize <- map$binsize
  max_cM_dist <- map$max_cM_dist
  min_nfrag_bin <- map$min_nfrag_bin
  gap_size <- map$gap_size
 }

 make_agp(hic_map_oriented, gap_size=gap_size, species=species)->a

 a$agp[, .(length=sum(scaffold_length)), key=agp_chr]->chrlen
 chrlen[, alphachr := sub("chr", "", agp_chr)]
 chrNames(species=species)[chrlen, on="alphachr"]->chrlen
 chrlen[, truechr := !grepl("Un", alphachr)]
 chrlen[order(!truechr, chr)]->chrlen
 chrlen[, offset := cumsum(c(0, length[1:(.N-1)]))]
 chrlen[, plot_offset := cumsum(c(0, length[1:(.N-1)]+1e8))]

 list(agp=a$agp, agp_bed=a$agp_bed, chrlen=chrlen, hic_map=hic_map_oriented, hic_map_bin=hic_map_bin)->res
 invisible(lapply(sort(c("min_nfrag_scaffold", "max_cM_dist", "binsize", "min_nfrag_bin", "gap_size")), function(i){
  res[[i]] <<- get(i)
 }))
 res
}

# Convert Hi-C map table into an AGP
make_agp<-function(hic_map_oriented, gap_size=100, species){

 hic_map_oriented[, .(scaffold, chr = consensus_chr,
	  popseq_cM=ifelse(consensus_chr == popseq_chr | is.na(consensus_chr), popseq_cM, NA),
	  scaffold_length = length, hic_bin, orientation=consensus_orientation)]->z
 chrNames(agp=T, species=species)[z, on="chr"]->z

 z[, agp_chr := "chrUn"]
 z[!is.na(hic_bin), agp_chr := sub("NA", "Un", paste0("chr", alphachr))]
 z[, alphachr := NULL]
 z[order(agp_chr, hic_bin, chr, popseq_cM, -scaffold_length)]->z
 z[, index := 2*1:.N-1]
 z[, gap := F]
 rbind(z, data.table(scaffold="gap", gap=T, chr=NA, popseq_cM=NA, scaffold_length = gap_size, hic_bin = NA, orientation = NA, agp_chr=z$agp_chr, index=z$index+1))->z
 z[order(index)][, head(.SD, .N-1), by=agp_chr]->z
 z[, agp_start := {if(nrow(.SD) > 1){
    cumsum(c(0, scaffold_length[1:(.N-1)]))+1
 } else {1}
 }, by = agp_chr]
 z[, agp_end := cumsum(scaffold_length), by = agp_chr]

 z[, .(scaffold=scaffold, bed_start=0, bed_end=scaffold_length, name=scaffold, score=1, strand=ifelse(is.na(orientation) | orientation == 1, "+", "-"), agp_chr=agp_chr)]->agp_bed

 list(agp=z, agp_bed=agp_bed)
}

# Read BED files with positions of restriction fragments on the input assembly, lift coordinates to updated assemblies
read_fragdata<-function(info, map_10x=NULL, assembly_10x=NULL, file){
 fragbed<-fread(file, head=F, col.names=c("orig_scaffold", "start", "end"))
 fragbed[, length := end - start]
 fragbed[, start := start + 1]
 info[, .(scaffold, start=orig_start, orig_start, orig_scaffold)][fragbed, on=c("orig_scaffold", "start"), roll=T]->fragbed
 fragbed[, start := start - orig_start + 1]
 fragbed[, end := end - orig_start + 1]
 fragbed[, orig_start := NULL]
 fragbed[, orig_scaffold := NULL]
 if(!is.null(assembly_10x)){
  map_10x$agp[gap == F, .(super, orientation, super_start, super_end, scaffold)][fragbed, on="scaffold"]->fragbed
  fragbed[orientation == 1, start := super_start - 1 + start]
  fragbed[orientation == 1, end := super_start - 1 + end]
  fragbed[orientation == -1, start := super_end - end + 1]
  fragbed[orientation == -1, end := super_end - start + 1]
  fragbed[, c("orientation", "super_start", "super_end", "scaffold") := list(NULL, NULL, NULL, NULL)]
  setnames(fragbed, "super", "orig_scaffold")

  assembly_10x$info[, .(scaffold, start=orig_start, orig_start, orig_scaffold)][fragbed, on=c("orig_scaffold", "start"), roll=T]->fragbed
  fragbed[, start := start - orig_start + 1]
  fragbed[, end := end - orig_start + 1]
  fragbed[, orig_start := NULL]
  fragbed[, orig_scaffold := NULL]

  fragbed[, .(nfrag = .N), keyby=scaffold][assembly_10x$info, on="scaffold"][is.na(nfrag), nfrag := 0]->z
 } else {
  fragbed[, .(nfrag = .N), keyby=scaffold][info, on="scaffold"][is.na(nfrag), nfrag := 0]->z
 }
 list(bed=fragbed[], info=z[])
}

# Add Hi-C link information to a hic_map object
add_psmol_fpairs<-function(assembly, hic_map, nucfile, map_10x=NULL, assembly_10x=NULL, cov=NULL){
 if(is.null(map_10x)){
  assembly$fpairs[, .(scaffold1, scaffold2, pos1, pos2)] -> z
 } else {
  assembly_10x$fpairs[, .(scaffold1, scaffold2, pos1, pos2)] -> z
 }
 hic_map$agp[agp_chr != "chrUn", .(chr, scaffold, orientation=orientation, agp_start, agp_end)]->a
 a[is.na(orientation), orientation := 1]
 setnames(copy(a), paste0(names(a), 1))[z, on="scaffold1"]->z
 setnames(copy(a), paste0(names(a), 2))[z, on="scaffold2"]->z
 z[orientation1 == 1, pos1 := agp_start1 - 1 + pos1]
 z[orientation1 == -1, pos1 := agp_end1 + 1 - pos1]
 z[orientation2 == 1, pos2 := agp_start2 - 1 + pos2]
 z[orientation2 == -1, pos2 := agp_end2 + 1 - pos2]
 z[!is.na(chr1) & !is.na(chr2), .(chr1, chr2, start1=pos1, start2=pos2)]->links

 n<-c("orig_scaffold", "orig_start", "orig_end", "frag_id", "nA", "nC", "nG", "nT", "nN", "length")
 nuc<-fread(nucfile, select=c(1:4,7:11,13), head=T, col.names=n)
 assembly$info[, .(scaffold, orig_scaffold, orig_start, off=orig_start)][nuc, on=c("orig_scaffold", "orig_start"), roll=T]->z
 z[, start := orig_start - off + 1]
 z[, end := orig_end - off + 1]
 if(!is.null(map_10x)){
  map_10x$agp[gap == F, .(super, orientation, super_start, super_end, scaffold)][z, on="scaffold"]->z
  z[orientation == 1, start := super_start - 1 + start]
  z[orientation == 1, end := super_start - 1 + end]
  z[orientation == -1, start := super_end - end + 1]
  z[orientation == -1, end := super_end - start + 1]
  z[, c("orientation", "super_start", "super_end", "scaffold") := list(NULL, NULL, NULL, NULL)]
  setnames(z, "super", "orig_scaffold")
  assembly_10x$info[, .(scaffold, start=orig_start, orig_start, orig_scaffold)][z, on=c("orig_scaffold", "start"), roll=T]->z
  z[, start := start - orig_start + 1]
  z[, end := end - orig_start + 1]
  z[, orig_start := NULL]
  z[, orig_scaffold := NULL]
 }
 a[z, on="scaffold", nomatch=0]->z
 z[orientation == 1, start := agp_start - 1 + start]
 z[orientation == -1, start := agp_end + 1 - start]
 z[orientation == 1, end := agp_start - 1 + end]
 z[orientation == -1, end := agp_end + 1 - end]
 z[, c("chr", "start", "end", n[4:length(n)]), with=F]->z
 z[, cov := 1]
 z->frags
 
 copy(hic_map) -> res
 res$links <- links
 res$frags <- frags
 res
}

# Create a heatmap plot of Hi-C contact matrix
contact_matrix<-function(hic_map, links, file, species, chrs=NULL, boundaries=T, grid=NULL, ncol=100, trafo=NULL, v=NULL){
 colorRampPalette(c("white", "red"))(ncol)->whitered

 if(is.null(chrs)){
  chrs <- unique(links$chr1) 
 }

 binsize <- min(links[dist > 0]$dist)

 links[chr1 %in% chrs, .(chr=chr1, bin1, bin2, l=log10(nlinks_norm))]->z
 if(is.null(trafo)){
  z[, col := whitered[cut(l, ncol, labels=F)]]
 } else {
  z[, col := whitered[cut(trafo(l), ncol, labels=F)]]
 }

 chrNames(agp=T, species=species)[z, on="chr"]->z

 pdf(file)
 lapply(chrs, function(i){
  z[chr == i, plot(0, las=1, type='n', bty='l', xlim=range(bin1/1e6), ylim=range(bin2/1e6), xlab="position (Mb)", ylab="position (Mb)", col=0, main=chrNames(species=species)[chr == i, alphachr])]
  if(boundaries){
   hic_map$agp[gap == T & agp_chr == chrNames(agp=T, species=species)[chr == i, agp_chr], abline(lwd=1, col='gray', v=(agp_start+agp_end)/2e6)]
   hic_map$agp[gap == T & agp_chr == chrNames(agp=T, species=species)[chr == i, agp_chr], abline(lwd=1, col='gray', h=(agp_start+agp_end)/2e6)]
  }
  if(!is.null(grid)){
   max(hic_map$agp[gap == T & agp_chr == chrNames(agp=T, species=species)[chr == i, agp_chr]]$agp_end)->end
   abline(v=seq(0, end, grid)/1e6, col="blue", lty=2)
  }
  if(!is.null(v)){
   abline(v=v/1e6, col="blue", lty=2)
  }
  z[chr == i, rect((bin1-binsize)/1e6, (bin2-binsize)/1e6, bin1/1e6, bin2/1e6, col=col, border=NA)]
 })
 dev.off()
}

# helper function to supply scaffolds order on the short of 5AS (a single bin in the Chapman et al. POPSEQ map of bread wheat)
correct_5A<-function(hic_info, bam, map, assembly){
 fread(cmd=paste("samtools view -q 30 -F260", bam, "| cut -f 1,3,4"))->i90sam
 setnames(i90sam, c("i90_marker", "scaffold", "pos"))

 i90map<-fread(map, select=c(1,3,4))
 setnames(i90map, c("i90_marker", "i90_alphachr", "i90_cM"))
 i90map[, n := .N, by=i90_marker]
 i90map[n == 1][, n := NULL][]->i90map
 setnames(chrNames(species="wheat"), c("i90_alphachr", "i90_chr"))[i90map, on="i90_alphachr"]->i90map
 i90map[i90sam, on="i90_marker"][i90_chr == 5, .(orig_scaffold=scaffold, orig_pos=pos, cM=i90_cM)]->z
 copy(z)->i90k_5A

 assembly$info[, .(scaffold, orig_scaffold, orig_start, orig_pos=orig_start)][z, on=c("orig_scaffold", "orig_pos"), roll=T]->z
 z[, .(scaffold, pos=orig_pos-orig_start+1, cM)]->i90k_5A

 z[, .(chr=5, icM=median(cM)), key="scaffold"][hic_info, on=c("scaffold", "chr")]->hic_info
 mx <- hic_info[chr == 5, max(na.omit(icM))]
 hic_info[chr == 5, cM := mx - icM]->hic_info
 hic_info[, icM := NULL][]
}

# Plot the biases computed by find_inversions() along the genome
plot_inversions<-function(inversions, hic_map, bad=c(), species, file, height=5, width=20){
 yy<-inversions$summary
 w<-inversions$ratio
 yy[scaffold %in% bad]->yy2

 pdf(file, height=height, width=width)
 lapply(sort(unique(w$chr)), function(i){
  w[chr == i, plot(bin/1e6, r, pch=20, las=1, bty='l', ylab="r", xlab="genomic position (Mb)", main=alphachr[i])]
  hic_map$agp[gap == T & agp_chr == chrNames(agp=T, species=species)[chr == i]$agp_chr, abline(lwd=1, col='gray', v=(agp_start+agp_end)/2e6)]
  if(nrow(yy2) > 0 & i %in% yy2$chr){
   yy2[chr == i, rect(col="#FF000011", agp_start/1e6, -1000, agp_end/1e6, 1000), by=scaffold]
  }
 })
 dev.off()
}

# Find inverted scaffolds based on directionality biases
find_inversions<-function(hic_map, links, species, chrs=NULL, cores=1, winsize=15, maxdist=1e8, threshold=40, factor=100){
 
 if(is.null(chrs)){
  unique(links$chr1) -> chrs
 }

 links[chr1 %in% chrs]->zz
 zz[, l := log(factor*nlinks_norm/sum(nlinks_norm)*1e6)]->zz
 hic_map$chrlen[, .(chr1=chr, length)][zz, on="chr1"]->zz
 zz[dist <= bin1 & dist <= bin2 & (length - bin1) >= dist & (length - bin2) >= dist]->zz
 zz[l >= 0 & dist <= maxdist, .(r = sum(sign(bin1 - bin2) * l)), key=.(chr=chr1, bin=bin1)]->w
 chrNames(agp=T, species=species)[w, on="chr"]->w

 data.table(w, rbindlist(mclapply(mc.cores=cores, chrs, function(i){
  ww<-w[chr == i][order(bin)]
  data.table(
   cc=ww[, rollapply(r, width=winsize, FUN=function(x) cor(x, 1:length(x)), align="left", fill=NA)],
   lm=ww[, rollapply(r, width=winsize, FUN=function(x) lm(data=.(a=x, b=1:winsize/winsize), a~b)$coefficient[2], align="left",fill=NA)]
  )
 })))->w

 hic_map$agp[gap == F, .(agp_chr, bin=agp_start, scaffold)]->x
 copy(w)->y
 y[, bin := bin + 1]
 x[y, on=c("agp_chr", "bin"), roll=T]->y
 y[agp_chr != "chrUn"]->y
 y[, .(s=sd(r)), key=scaffold][!is.na(s)][order(-s)]->yy
 hic_map$hic_map[, .(scaffold, chr=consensus_chr, hic_bin)][yy, on="scaffold"]->yy
 hic_map$agp[, .(scaffold, agp_start, agp_end)][yy, on="scaffold"]->yy
 list(ratio=w, summary=yy)
}
 
# Flip the orientation of groups of scaffolds
correct_multi_inversions<-function(hic_map, ranges, species){
 copy(ranges)->y
 chrNames(agp=T, species=species)[, .(consensus_chr=chr, agp_chr)][hic_map$hic_map, on="consensus_chr"]->m
 y[, i:=1:.N]

 y[, .(b=seq(start,end)), by=.(agp_chr, i)][, .N, key=.(agp_chr, b)][N > 1]->dups
 if(nrow(dups) > 0){
  stop(paste("Overlapping ranges specified: bin(s)",
        paste(dups[, paste(sep=":", agp_chr, b)], collapse=", "), "are contained in more than one inversion."))
 }

 y[, .(agp_chr=agp_chr[1], hic_bin = start:end, new_bin = end:start, new=T), by=i][, i := NULL][m, on=c("agp_chr", "hic_bin")][is.na(new), new := F]->m
 m[new == T, hic_bin := new_bin]
 m[new == T, consensus_orientation := consensus_orientation * -1]
 m[new == T & is.na(consensus_orientation), consensus_orientation := -1]
 m[, c("new", "new_bin", "agp_chr") := list(NULL, NULL, NULL)]
 make_agp(m, gap_size=hic_map$gap_size, species=species)->a

 list()->new
 new$hic_map <- m
 new$agp <- a$agp
 new$gap_size <- copy(hic_map$gap_size)
 new$chrlen <- copy(hic_map$chrlen)
 new$binsize <- copy(hic_map$binsize)
 new$hic_map_bin <- copy(hic_map$hic_map_bin)
 new$max_cM_dist <- copy(hic_map$max_cM_dist)
 new$min_nfrag_bin <- copy(hic_map$min_nfrag_bin)
 new$agp_bed <- a$agp_bed
 new$corrected_multi_inversions <- copy(ranges)
 new
}

lift_nucfile<-function(assembly, map_10x, nucfile, outfile){
 n <- c("orig_scaffold", "orig_start", "orig_end", "frag_id", "nA", "nC", "nG", "nT", "nN", "length")
 nuc <- fread(nucfile, select=c(1:4,7:11,13), head=T, col.names=n)
 assembly$info[, .(scaffold, orig_scaffold, orig_start, off=orig_start)][nuc, on=c("orig_scaffold", "orig_start"), roll=T]->z
 z[, start := orig_start - off + 1]
 z[, end := orig_end - off + 1]
 map_10x$agp[gap == F, .(super, orientation, super_start, super_end, scaffold)][z, on="scaffold"]->z
 z[orientation == 1, start := super_start - 1 + start]
 z[orientation == 1, end := super_start - 1 + end]
 z[orientation == -1, start := super_end - end + 1]
 z[orientation == -1, end := super_end - start + 1]
 z[, .(super, start - 1, end, frag_id, ".", ".", nA, nC, nG, nT, nN, ".", length)]->zz
 fwrite(zz, file=outfile, col.names=T, row.names=F, sep="\t", quote=F)
}

# Export Hi-C map information for use in the R Shiny app
inspector_export<-function(hic_map, assembly, inversions, species, file, trafo=NULL){
 hic_map$hic_1Mb$norm -> x
 if(!is.null(trafo)){
  x[, nlinks_norm := trafo(nlinks_norm)]
 }
 chrNames(agp=T, species=species)[, .(chr1=chr, agp_chr)][x, on="chr1"]->x
 x[, c("chr1", "chr2", "id2", "id1", "dist") := list(NULL, NULL, NULL, NULL, NULL)]
 
 inversions$ratio[, .(agp_chr, bin, r)] -> y

 copy(hic_map$agp)[, chr := NULL][]->a
 assembly$info[, .(scaffold, popseq_chr=popseq_alphachr)][a, on="scaffold"]->a
 saveRDS(version=2, list(agp=a, binsize=1e6, links=x, ratio=y), file=file)
}

# Plot heatmap of interchromosomal Hi-C contact matrix
interchromosomal_matrix_plot<-function(hic_map, binsize=1e6, file, species, resolution=72,
				       height=3000, width=3000, col=NULL, ncol=100, cex=2,
				       oma=c(5,7,5,1)){
 colorRampPalette(c("white", "red"))(ncol)->whitered
 copy(hic_map$hic_1Mb$mat)->z
 z[, nl := log10(nlinks)]
 if(is.null(col)){
  z[, col := whitered[cut(nl, ncol, labels=F)]]
 } else{
  cc <- col
  z[, col := cc] 
 }

 nrow(chrNames(species=species)) -> n

 z[, chr1 := as.character(chr1)]
 z[, chr2 := as.character(chr2)]

 copy(hic_map$chrlen)[, chr := paste(chr)][paste(1:n), on="chr"]$length -> ll

 png(file, res=resolution, height=height, width=width)
 par(mar=c(0.1, 0.1, 0.1, 0.1))
 layout(matrix(ncol=n, nrow=n, 1:(n*n), byrow=F), widths=ll, heights=ll)
 par(oma=oma)
 lapply(1:n, function(i){
  lapply(1:n, function(j){
   z[paste(i), on="chr1"][paste(j), on="chr2"] -> zz
   zz[, plot(0, las=1, xaxt='n', yaxt='n', type='n', xlim=range(bin1/1e6), ylim=range(bin2/1e6), xlab="", ylab="", col=0, main="")]
   zz[, rect((bin1-binsize)/1e6, (bin2-binsize)/1e6, bin1/1e6, bin2/1e6, col=col, border=NA)]
   if(j == 1){
    title(xpd=NA, chrNames(species=species)[chr == i]$alphachr, line=1, cex.main=cex, font=2)
   }
   if(i == 1){
    mtext(chrNames(species=species)[chr == j]$alphachr, side=2, line=1.5, cex=cex*0.6, font=2, las=1, xpd=NA)
   }
  })
 })
 dev.off()
}

# Bin Hi-C links in windows of a fixed size
bin_hic_step<-function(hic, frags, binsize, step=NULL, chrlen, chrs=NULL, cores=1){

 if(is.null(chrs)){
  chrs <- unique(chrlen$chr)
 }

 if(is.null(step)){
  step <- binsize
 }

# Assign chromosome names for alphachr values starting with GDB136_
# Assign chromosome numbers for alphachr values starting with GDB136_
 chrlen[grepl("^GDB136_", alphachr), chr := as.integer(sub("^GDB136_(\\d+).*$", "\\1", alphachr))]
 chrlen[grepl("^GDB136_", alphachr), alphachr := as.character(chr)]
 chrlen[alphachr == "Un", chr := 0]

 # From chrlen, take alphacr and build chrs 
 chrs <- unique(chrlen[truechr == TRUE, chr]) # added by JOAN, exclude 0

 options(scipen=20)
 print("loading chrlen after chr change")
 print(chrs)
 print(chrlen)
 chrlen[chr %in% chrs, .(bin=seq(0, length, step)), key=chr][, end := bin + binsize - 1][, start := bin] -> bins
 print("loading bin")
 bins[, id := paste0(chr, ":", bin)]
 print("loaded bins:")
 print(bins)
 
 print("loading hic")
 hic[chr1 %in% chrs & chr2 %in% chrs] -> z
 print ("loading hic done")
 print(z)

 # print step and binsize
 print(paste("step:", step, "binsize:", binsize))

 # Transform all chr1 and chr2 to integers
      z[, chr1 := as.integer(chr1)]
      z[, chr2 := as.integer(chr2)]
      z[, start1 := as.integer(start1)]
      z[, start2 := as.integer(start2)]

 if(step == binsize){
  z[, bin1 := as.integer(start1 %/% binsize * binsize)]
  z[, bin2 := as.integer(start2 %/% binsize * binsize)]
  z[, .(nlinks=.N), keyby=.(chr1,bin1,chr2,bin2)]->z
  z[, id1 := paste(sep=":", chr1, bin1)]
  z[, id2 := paste(sep=":", chr2, bin2)]
  z->mat

  frags[chr %in% chrs]->f
  f[, bin := as.integer(start %/% binsize * binsize)]
  f[, .(nfrags=.N, eff_length=sum(length), cov=weighted.mean(cov, length), 
	gc=(sum(nC)+sum(nG))/(sum(nA)+sum(nC)+sum(nG)+sum(nT))), keyby=.(chr, bin)]->f
  f[, id := paste(sep=":", chr, bin)]
 } else {
  copy(bins)->b
  b[, idx := 1:.N]
  b[, .(win=seq(bin, bin+binsize-1, step)),.(chr, bin)]->b

  z[chr1 == chr2]->z
  z[, win1:=as.integer(start1 %/% step * step)]
  z[, win2:=as.integer(start2 %/% step * step)]

  rbindlist(mclapply(mc.cores=cores, mc.preschedule=F, chrs, function(i){
   z[chr1 == i]->y
   y[, .(nlinks=.N), keyby=.(chr1,win1,chr2,win2)]->y
   setnames(copy(b), paste0(names(b), 1))[y, on=c("chr1", "win1"), allow.cartesian=T]->y
   setnames(copy(b), paste0(names(b), 2))[y, on=c("chr2", "win2"), allow.cartesian=T]->y
   y[, .(nlinks=.N), keyby=.(chr1,bin1,chr2,bin2)]
  }))->y
  y[, id1:=stri_c(sep=":", chr1, bin1)]
  y[, id2:=stri_c(sep=":", chr2, bin2)]
  y[, dist := abs(bin1 - bin2)]
  y->mat
  
  copy(frags)->f
  f[, win := as.integer(start %/% step * step)]
  b[f, on=c("chr", "win"), nomatch=0, allow.cartesian=T]->f
  f[, .(nfrags=.N, eff_length=sum(length), cov=weighted.mean(cov, length), 
	gc=(sum(nC)+sum(nG))/(sum(nA)+sum(nC)+sum(nG)+sum(nT))), keyby=.(chr, bin)]->f
  f[, id :=paste(sep=":", chr, bin)]
 }

 list(bins=f[], mat=mat[])
}

# Inter-chromosomal normalization, taken from Hi-C norm https://academic.oup.com/bioinformatics/article/28/23/3131/192582
normalize_mat_trans<-function(u, v_a, v_b){
#change matrix into vector
 u_vec<-c(as.matrix(u))
#get cov matrix
 len_m<-as.matrix(log(v_a[,eff_length]%o%v_b[,eff_length]))
 gcc_m<-as.matrix(log(v_a[,gc]%o%v_b[,gc]))
 map_m<-as.matrix(log(v_a[,cov]%o%v_b[,cov]))
#centralize cov matrix of enz, gcc
 len_m<-(len_m-mean(c(len_m)))/sd(c(len_m))
 gcc_m<-(gcc_m-mean(c(gcc_m)))/sd(c(gcc_m))
#change matrix into vector
 len_vec<-c(len_m)
 gcc_vec<-c(gcc_m)
 map_vec<-c(map_m)
#fit Poisson regression: u~len+gcc+offset(map)
 fit<-glm(u_vec~len_vec+gcc_vec+offset(map_vec),family="poisson")
#summary(fit)
 coeff <- round(fit$coeff,4)
 res <- round(u/exp(coeff[1]+coeff[2]*len_m+coeff[3]*gcc_m+map_m), 4)
 data.table(id1=v_a$id, res)->res
 melt(res, id.var="id1", value.name="nlinks_norm", variable.name="id2")->res
 res[nlinks_norm > 0]
}

# Intra-chromosomal normalization, taken from HiCNorm https://academic.oup.com/bioinformatics/article/28/23/3131/192582
normalize_mat_cis<-function(u, v){
 u_vec<-u[upper.tri(u,diag=F)]

#get cov matrix
 len_m<-as.matrix(log(v[, eff_length]%o%v[, eff_length]))
 gcc_m<-as.matrix(log(v[, gc]%o%v[, gc]))
 map_m<-as.matrix(log(v[, cov]%o%v[, cov]))

#centralize cov matrix of enz, gcc
 len_m<-(len_m-mean(c(len_m)))/sd(c(len_m))
 gcc_m<-(gcc_m-mean(c(gcc_m)))/sd(c(gcc_m))

#change matrix into vector
 len_vec<-len_m[upper.tri(len_m,diag=F)]
 gcc_vec<-gcc_m[upper.tri(gcc_m,diag=F)]
 map_vec<-map_m[upper.tri(map_m,diag=F)]

#fit Poisson regression: u~len+gcc+offset(map)
 fit<-glm(u_vec~len_vec+gcc_vec+offset(map_vec),family="poisson")

#summary(fit)
 coeff <- round(fit$coeff,4)
 res <- round(u/exp(coeff[1]+coeff[2]*len_m+coeff[3]*gcc_m+map_m), 4)
 data.table(id1=colnames(res), res)->res
 melt(res, id.var="id1", value.name="nlinks_norm", variable.name="id2")->res
 res[nlinks_norm > 0]
}

# format intrachromosomal Hi-C matrix for HiCNorm 
normalize_cis<-function(binhic, ncores=1, chrs=NULL, percentile=0, omit_smallest=0){
 if(is.null(chrs)){
  chrs <- unique(binhic$bins$chr)
 }

 mf <- binhic$bins[chr %in% chrs]
 ab <- binhic$mat

 if(percentile > 0){
  mf[eff_length >= quantile(eff_length, 0:100/100)[percentile + 1]]->mf
 }

 if(omit_smallest > 0){
  setorder(mf, eff_length)
  mf[, omit_idx := 1:.N, by=chr]
  mf[omit_idx > omit_smallest][, omit_idx := NULL] -> mf
 }

 rbindlist(mclapply(mc.cores=ncores, chrs, function(i) {
  ab[chr1 == chr2 & chr1 == i & id1 %in% mf$id & id2 %in% mf$id]->abf
  dcast.data.table(abf, id1 ~ id2, value.var="nlinks", fill=as.integer(0))->mat

  u<-as.matrix(mat[, setdiff(names(mat), "id1"), with=F])
  setkey(mf, id)[colnames(u)]->v
  normalize_mat_cis(u,v)
 }))->nhic
 mf[, data.table(key="id1", id1=id, chr1=chr, bin1=bin)][setkey(nhic, id1)]->nhic
 mf[, data.table(key="id2", id2=id, chr2=chr, bin2=bin)][setkey(nhic, id2)]->nhic
 nhic[, dist := abs(bin2 - bin1)][]
}

# format interchromosomal Hi-C matrix for HiCNorm 
normalize_trans<-function(binhic, ncores=1, chrs=NULL, percentile=0){
 if(is.null(chrs)){
  chrs <- unique(binhic$bins$chr)
 }

 mf <- binhic$bins[chr %in% chrs]
 ab <- binhic$mat

 if(percentile > 0){
  mf[eff_length >= quantile(eff_length, 0:100/100)[percentile + 1]]->mf
 }

 rbindlist(mclapply(mc.cores=ncores, chrs, function(i) rbindlist(lapply(setdiff(1:21, i), function(j) {
  ab[chr1 == i & chr2 == j & id1 %in% mf$id & id2 %in% mf$id]->abf
  dcast.data.table(abf, id1 ~ id2, value.var="nlinks", fill=as.integer(0))->mat
  setkey(mf, id)[mat$id1]->v_a
  u<-as.matrix(mat[, setdiff(names(mat), "id1"), with=F])
  setkey(mf, id)[colnames(u)]->v_b
  normalize_mat_trans(u, v_a, v_b)
 }))))->nhic
 mf[, data.table(key="id1", id1=id, chr1=chr, bin1=bin)][setkey(nhic, id1)]->nhic
 mf[, data.table(key="id2", id2=id, chr2=chr, bin2=bin)][setkey(nhic, id2)]->nhic
}

# Like scaffold_10x, but using a Hi-C map of the scaffolds, instead of a genetic map to prune erroneous edges
scaffold_10x_hic <- function(assembly, prefix="super", min_npairs=5, max_dist=1e5, min_nmol=6, min_nsample=2, popseq_dist=5, max_dist_orientation=5, ncores=1, verbose=T, raw=F, unanchored=T, hic_map=NULL){

 make_agp<-function(membership, gap_size=100){

  membership[, .(scaffold, length, super, bin, orientation)]->z

  setorder(z, super, bin, -length)
  z[, index := 2*1:.N-1]
  z[, gap := F]
  rbind(z, data.table(scaffold="gap", gap=T, super=z$super, bin=NA, length = gap_size,orientation = NA, index=z$index+1))->z
  z[order(index)][, head(.SD, .N-1), by=super]->z
  z[, n := .N, key=super]
  z[n > 1, super_start := cumsum(c(0, length[1:(.N-1)])) + 1, by = super]
  z[n == 1, super_start := 1]
  z[, super_end := cumsum(length), by = super]
  z[, n := NULL]

  z[, .(scaffold=scaffold, bed_start=0, bed_end=length, name=scaffold, score=1, strand=ifelse(is.na(orientation) | orientation == 1, "+", "-"), super=super)]->agp_bed

  list(agp=z, agp_bed=agp_bed)
 }

 make_super_scaffolds<-function(links, info, excluded=c(), ncores, prefix){
  info[, .(scaffold, chr=popseq_chr, cM=popseq_cM, length)]->info2
  excluded -> excluded_scaffolds
  info2[, excluded := scaffold %in% excluded_scaffolds] -> input
  setnames(input, "scaffold", "cluster")
  setnames(copy(links), c("scaffold1", "scaffold2"), c("cluster1", "cluster2"))->hl
  make_super(hl=hl, cluster_info=input, verbose=F, prefix=prefix, cores=ncores, paths=T, path_max=0, known_ends=F, maxiter=100)->s
  copy(s$membership) -> m
  setnames(m, "cluster", "scaffold")

  max(as.integer(sub(paste0(prefix, "_"), "", s$super_info$super))) -> maxidx
  rbind(m, info[!m$scaffold, on="scaffold"][, .(scaffold, bin=1, rank=0, backbone=T, chr=popseq_chr, cM=popseq_cM,
					length=length, excluded=scaffold %in% excluded_scaffolds, super=paste0(prefix, "_", maxidx + 1:.N))])->m

  m[, .(n=.N, nbin=max(bin), max_rank=max(rank), length=sum(length)), key=super]->res
  res[, .(super, super_size=n, super_nbin=nbin)][m, on="super"]->m
  list(membership=m, info=res)
 }

 if(verbose){
  cat("Finding links.\n")
 }

 copy(assembly$info) -> info
 if(!is.null(hic_map)){
  info[, c("popseq_chr", "popseq_cM") := list(NULL, NULL)]
  hic_map$agp[, .(popseq_chr = chr, popseq_cM = hic_bin, scaffold)][info, on="scaffold"]->info
 } 

 assembly$molecules[npairs >= min_npairs] -> z
 info[, .(scaffold, scaffold_length=length)][z, on="scaffold"]->z
 z[end <= max_dist | scaffold_length - start <= max_dist]->z
 z[, nsc := length(unique(scaffold)), key=.(barcode, sample)]
 z[nsc >= 2]->z
 z[, .(scaffold1=scaffold, npairs1=npairs, pos1=as.integer((start+end)/2), sample, barcode)]->x 
 z[, .(scaffold2=scaffold, npairs2=npairs, pos2=as.integer((start+end)/2), sample, barcode)]->y
 y[x, on=.(sample, barcode), allow.cartesian=T][scaffold1 != scaffold2]->xy
 xy -> link_pos
 xy[, .(nmol=.N), key=.(scaffold1, scaffold2, sample)]->w
 w[nmol >= min_nmol]->ww
 ww[, .(nsample = length(unique(sample))), key=.(scaffold1, scaffold2)][nsample >= min_nsample]->ww2
 info[, .(scaffold1=scaffold, popseq_chr1=popseq_chr, length1=length, popseq_pchr1=popseq_pchr, popseq_cM1=popseq_cM)][ww2, on="scaffold1"]->ww2
 info[, .(scaffold2=scaffold, popseq_chr2=popseq_chr, length2=length, popseq_pchr2=popseq_pchr, popseq_cM2=popseq_cM)][ww2, on="scaffold2"]->ww2
 ww2[, same_chr := popseq_chr2 == popseq_chr1]
 ww2[, weight := -1 * log10((length1 + length2) / 1e9)]

 if(popseq_dist > 0){
  ww2[(popseq_chr2 == popseq_chr1 & abs(popseq_cM1 - popseq_cM2) <= popseq_dist)] -> links
 } else {
  ww2 -> links
 }

 ex <- c()

 run <- T

 if(verbose){
  cat("Finding initial super-scaffolds.\n")
 }
 # remove nodes until graph is free of branches of length > 1
 if(!raw){
  while(run){
   make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
   copy(out$membership) -> m
   copy(out$info) -> res
   m[unique(m[rank > 1][, .(super, bin)]), on=c("super", "bin")]->a
   a[rank == 0] -> add
   if(nrow(add) == 0){
    run <- F
   } else{
    c(ex, add$scaffold) -> ex
   }
  }
 } else {
  make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
  out$m -> m
  out$info -> res
 }

 if(!raw){
  if(verbose){
   cat("Tip removal.\n")
  }
  # remove short tips of rank 1 
  links[!scaffold2 %in% ex][, .(degree=.N), key=.(scaffold=scaffold1)]->b
  a <- b[m[m[rank == 1][, .(super=c(super, super, super), bin=c(bin, bin-1, bin+1))], on=c("super", "bin")], on="scaffold"]
  a[degree == 1 & length <= 1e4]$scaffold->add
  if(length(add) > 0){
   ex <- c(ex, add)
   make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
   out$m -> m
  }

  # remove short tips/bulges of rank 1 
  m[rank == 1 & length <= 1e4]$scaffold -> add
  ex <- c(ex, add)
  if(length(add) > 0){
   ex <- c(ex, add)
   make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
   out$m -> m
  }

  # resolve length-one-bifurcations at the ends of paths
  m[rank > 0][bin == 2 | super_nbin - 1 == bin ][, .(super, super_nbin, type = bin == 2, scaffold, length, bin0=bin)]->x
  unique(rbind(
  m[x[type == T, .(super, bin0, bin=1)], on=c("super", "bin")],
  m[x[type == F, .(super, bin0, bin=super_nbin)], on=c("super", "bin")]
  ))->a
  a[, .(super, bin0, scaffold2=scaffold, length2=length)][x, on=c("super", "bin0")][, ex := ifelse(length >= length2, scaffold2, scaffold)]$ex -> add

  if(length(add) > 0){
   ex <- c(ex, add)
   make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
   out$membership -> m
  }

  # remove short tips/bulges of rank 1 
  m[rank == 1 & length <= 1e4]$scaffold -> add
  if(length(add) > 0){
   ex <- c(ex, add)
   make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
   out$m -> m
  }
  
  # remove tips of rank 1 
  links[!scaffold2 %in% ex][, .(degree=.N), key=.(scaffold=scaffold1)]->b
  b[m[rank == 1], on="scaffold"][degree == 1]$scaffold -> add
  if(length(add) > 0){
   ex <- c(ex, add)
   make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
   out$m -> m
  }

  # remove remaining nodes of rank > 0
  m[rank > 0]$scaffold -> add
  if(length(add) > 0){
   ex <- c(ex, add)
   make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
   out$m -> m
  }

  if(popseq_dist > 0 & unanchored == T){
   if(verbose){
    cat("Including unanchored scaffolds.\n")
   }
   # use unanchored scaffolds to link super-scaffolds
   ww2[is.na(popseq_chr1), .(scaffold_link=scaffold1, link_length=length1, scaffold1=scaffold2)]->x
   ww2[is.na(popseq_chr1), .(scaffold_link=scaffold1, scaffold2=scaffold2)]->y
   x[y, on="scaffold_link", allow.cartesian=T][scaffold1 != scaffold2]->xy

   m[, .(scaffold1=scaffold, super1=super, chr1=chr, cM1=cM, size1=super_nbin, d1 = pmin(bin - 1, super_nbin - bin))][xy, on="scaffold1"]->xy
   m[, .(scaffold2=scaffold, super2=super, chr2=chr, cM2=cM, size2=super_nbin, d2 = pmin(bin - 1, super_nbin - bin))][xy, on="scaffold2"]->xy
   xy[super2 != super1 & d1 == 0 & d2 == 0 & size1 > 1 & size2 > 1 & chr1 == chr2]->xy
   xy[scaffold1 < scaffold2, .(nscl=.N), scaffold_link][xy, on="scaffold_link"]->xy
   xy[nscl == 1] -> xy
   xy[super1 < super2][, c("n", "g"):=list(.N, .GRP), by=.(super1, super2)][order(-link_length)][!duplicated(g)]->zz

   sel <- zz[, .(scaffold1=c(scaffold_link, scaffold_link, scaffold1, scaffold2),
	  scaffold2=c(scaffold1, scaffold2, scaffold_link, scaffold_link))]
   rbind(links, ww2[sel, on=c("scaffold1", "scaffold2")])->links2

   make_super_scaffolds(links=links2, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
   out$m -> m

   #resolve branches
   m[rank > 0][bin == 2 | super_nbin - 1 == bin ][, .(super, super_nbin, type = bin == 2, scaffold, length, bin0=bin)]->x
   unique(rbind(
   m[x[type == T, .(super, bin0, bin=1)], on=c("super", "bin")],
   m[x[type == F, .(super, bin0, bin=super_nbin)], on=c("super", "bin")]
   ))->a
   a[, .(super, bin0, scaffold2=scaffold, length2=length)][x, on=c("super", "bin0")][, ex := ifelse(length >= length2, scaffold2, scaffold)]$ex -> add

   if(length(add) > 0){
    ex <- c(ex, add)
    make_super_scaffolds(links=links2, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
    out$membership -> m
   }

   # remove remaining nodes of rank > 0
   m[rank > 0]$scaffold -> add
   if(length(add) > 0){
    ex <- c(ex, add)
    make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
   }
  }

  out$m -> m
  out$info -> res
 } 

 if(verbose){
  cat("Orienting scaffolds.\n")
 }
 #orient scaffolds
 m[super_nbin > 1, .(scaffold1=scaffold, bin1=bin, super1=super)][link_pos, on="scaffold1", nomatch=0]->a
 m[super_nbin > 1, .(scaffold2=scaffold, bin2=bin, super2=super)][a, on="scaffold2", nomatch=0]->a
 a[super1 == super2 & bin1 != bin2]->a
 a[, d := abs(bin2 - bin1)]->a
 a[d <= max_dist_orientation]->a
 a[, .(nxt=mean(pos1[bin2 > bin1]), prv=mean(pos1[bin2 < bin1])), key=.(scaffold=scaffold1)]->aa
 info[, .(scaffold, length)][aa, on="scaffold"]->aa
 aa[!is.nan(prv) & !is.nan(nxt), orientation := ifelse(prv <= nxt, 1, -1)]
 aa[is.nan(prv), orientation := ifelse(length - nxt <= nxt, 1, -1)]
 aa[is.nan(nxt), orientation := ifelse(length - prv <= prv, -1, 1)]

 aa[, .(orientation, scaffold)][m, on="scaffold"]->m
 m[, oriented := T]
 m[is.na(orientation), oriented := F]
 m[is.na(orientation), orientation := 1]
 setorder(m, super, bin, rank)

 m[, super_pos := 1 + cumsum(c(0, length[-.N])), by=super]

 if(verbose){
  cat("Anchoring super-scaffolds.\n")
 }
 # assign super scaffolds to genetic positions
 m[!is.na(chr), .(nchr=.N), key=.(chr, super)][, pchr := nchr/sum(nchr), by=super][order(-nchr)][!duplicated(super)]->y
 m[!is.na(cM)][y, on=c("chr", "super")]->yy
 y[res, on="super"]->res
 yy[, .(cM=mean(cM), min_cM=min(cM), max_cM=max(cM)), key=super][res, on='super']->res
 setorder(res, -length)

 make_agp(m, gap=100)->a

 list(membership=m, info=res, 
      agp=a$agp, agp_bed=a$agp_bed)
}

#function to plot contact matrix, Hi-C physical coverage and directionality bias, one page for each chromosome
big_hic_plot<-function(hic_map, cov=NULL, inv, species, file, chrs=NULL, cores=1, breaks=NULL){
 hmap <- hic_map

 print("inside big_hic_plot")

 ncol=100
 colorRampPalette(c("white", "red"))(ncol)->whitered
 links <- hmap$hic_1Mb$norm
 if(is.null(chrs)){
  chrs <- unique(links$chr1)
 }
 links[chr1 %in% chrs, .(chr=chr1, bin1, bin2, l=log10(nlinks_norm))]->z
 chrNames(agp=T, species=species)[z, on="chr"]->z
 binsize <- min(links[dist > 0]$dist)
 z[, col := whitered[cut(l, ncol, labels=F)]]
 nchr <- length(chrs)

 tempdir() -> d
 tmp <- paste0(d, "/plot")

 print ("checkpoint 1   tempdir")
 print("temp dir is:")
 print(d)
 mclapply(mc.cores=cores, chrs, function(i){
  f <- paste0(tmp, "_", chrNames(agp=T, species=species)[chr == i, agp_chr], ".png")
  png(f, res=500, height=5/3*2500, width=5/3*1500)
  layout(matrix(1:(3*1), ncol=1), heights=rep(c(3,1,1), 1))
   
   z[chr == i, plot(0, las=1, type='n', bty='l', xlim=range(bin1/1e6), ylim=range(bin2/1e6), xlab="position (Mb)", ylab="position (Mb)", col=0, main=chrNames(species=species)[chr == i, alphachr])]
    hmap$agp[gap == T & agp_chr == chrNames(agp=T, species=species)[chr == i, agp_chr], abline(lwd=1, col='gray', v=(agp_start+agp_end)/2e6)]
    hmap$agp[gap == T & agp_chr == chrNames(agp=T, species=species)[chr == i, agp_chr], abline(lwd=1, col='gray', h=(agp_start+agp_end)/2e6)]
   z[chr == i, rect((bin1-binsize)/1e6, (bin2-binsize)/1e6, bin1/1e6, bin2/1e6, col=col, border=NA)]

   if(!is.null(cov)){
    cov[chr == i][, plot(bin/1e6, type='n', las=1, xlab="genomic position (Mb)",
			ylab="r", bty='l', r)]
    hmap$agp[chr == i, abline(col="gray", v=c(0,agp_end/1e6))]
    if(!is.null(breaks)){
     abline(col="blue", lty=2, v=breaks[chr == i, bin/1e6])
    }
    abline(h=-2, lty=2, col='red')
    cov[chr == i][, points(bin/1e6, r, pch=".")]
    title(main="physical Hi-C coverage")
   } else {
    plot(axes=F, 0, type='n', xlab="", ylab="") 
   }

   inv$ratio[chr == i][, plot(bin/1e6, type='n', las=1, 
			      xlab="genomic position (Mb)",
			      ylab="r", bty='l', r)]
   hmap$agp[chr == i, abline(col="gray", v=c(0,agp_end/1e6))]
   inv$ratio[chr == i][, points(bin/1e6, pch=".", r)]
   title(main="directionality bias")
  dev.off()
  print(paste("created", f))
  # system(paste("convert", f, sub(".png$", ".pdf", f)))    doesnt fucking work due to permisions :)
 system(paste("img2pdf", f, "-o", sub(".png$", ".pdf", f)))
 })

 print("checkpoint 2")

 system(paste0("gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=", file, " ", 
	       paste(paste0(tmp, "_", 
			    setdiff(chrNames(agp=T, species=species)[chr %in% chrs]$agp_chr, "chrUn"), ".pdf"), collapse=' ')))
 # unlink(tmp, force=T) # commented to avoid deleting the temporary files
}

#calculate Hi-C physical coverage for pseudomolecules
hic_cov_psmol <- function(hic_map, binsize=1e3, binsize2=1e5, maxdist=1e6, 
			  cores=1){
 fpairs <- copy(hic_map$links)
 info <- hic_map$chrlen
 setnames(fpairs, c("start1", "start2"), c("pos1", "pos2"))

 print("inside hic_cov_psmol")

 fpairs[chr1 == chr2 & pos1 < pos2][, .(chr = chr1, bin1 = pos1 %/% binsize * binsize, bin2 =pos2 %/% binsize * binsize)]->f
 f[bin2 - bin1 > 2*binsize & bin2 - bin1 <= maxdist]->f
 f[, i := 1:.N]
 f[, b := paste0(chr, ":", bin1 %/% binsize2)]
 setkey(f, b)

 print("inside hic_cov_psmol 2")

 rbindlist(mclapply(mc.cores=cores, unique(f$b), function(j){
  f[j][, .(chr, bin=seq(bin1+binsize, bin2-binsize, binsize)), key=i][, .(n=.N), key=.(chr, bin)]
 }))->ff

 print("created ff") 
 print(ff)
 print("ff columns")
 print(colnames(ff))
 print(head(ff, 1))

 ff[, .(n=sum(n)), key=.(chr, bin)]->ff
 print("Checkpoint 1")
 print("Check chr and length from info")
 print(head(info, 1))

 # Ensure both chr are integer
 info[, chr := as.integer(chr)]
 ff[, chr := as.integer(chr)]

 info[, .(chr, length)][ff, on="chr"]->ff
      print("Checkpoint 2")
 ff[, d := pmin(bin, (length-bin) %/% binsize * binsize)]
 print("Checkpoint 3")
 ff[, nbin := .N, key="chr"]
 print("Checkpoint 4")
 ff[, mn := mean(n), key=d]
      print("Checkpoint 5")
 ff[, r := log2(n/mn)][]

 print("returning ff")
 print(ff)
}

#compare two Hi-C maps
diff_hic_map <- function(rds, new, species){
 old <- readRDS(rds)
 new$hic_map[, .(scaffold, new_chr=consensus_chr, new_bin=hic_bin, new_orientation=consensus_orientation)] -> x
 old$hic_map[, .(scaffold, length, old_chr=consensus_chr, old_bin=hic_bin, old_orientation=consensus_orientation)] -> y
 x[y, on="scaffold"] -> xy
 xy[, diff_chr :=  ((new_chr != old_chr & !is.na(new_chr) & !is.na(old_chr)) | (is.na(new_chr) & !is.na(old_chr)) | (!is.na(new_chr) & is.na(old_chr)))]
 xy[, diff_bin :=  ((new_bin != old_bin & !is.na(new_bin) & !is.na(old_bin)) | (is.na(new_bin) & !is.na(old_bin)) | (!is.na(new_bin) & is.na(old_bin)))]
 xy[, diff_orientation :=  ((new_orientation != old_orientation & !is.na(new_orientation) & !is.na(old_orientation)) | (is.na(new_orientation) & !is.na(old_orientation)) | (!is.na(new_orientation) & is.na(old_orientation)))]
 xy[diff_chr | diff_bin | diff_orientation] -> xy
 setnames(xy, "new_chr", "new_chrnum")
 setnames(xy, "old_chr", "old_chrnum")
 chrNames(species=species)[, .(old_chrnum=chr, old_chr=alphachr)][xy, on="old_chrnum"] -> xy
 chrNames(species=species)[, .(new_chrnum=chr, new_chr=alphachr)][xy, on="new_chrnum"] -> xy
 setorder(xy, new_chr, new_bin)
 setcolorder(xy, c("scaffold", "length", "diff_chr", "new_chr", "old_chr", "new_chrnum", "old_chrnum", "diff_bin", "new_bin", "old_bin", "diff_orientation", "new_orientation", "old_orientation"))
 xy
}

#write Hi-C map as Excel file for manual editing
write_hic_map <- function(rds, file, species=NULL){
 nmap_o <- readRDS(rds)
 nmap <- nmap_o$hic_map[, .(scaffold, length, old_chrnum=consensus_chr, old_bin=hic_bin, old_orientation=consensus_orientation)]
 chrNames(species=species)[, .(old_chrnum=chr, old_chr=alphachr)][nmap, on="old_chrnum"] -> nmap
 nmap[is.na(old_chr), old_chr := "Un"]
 nmap[is.na(old_chrnum), old_chr := 0]
 setorder(nmap, old_chrnum, old_bin)
 nmap[, new_orientation := old_orientation]
 setcolorder(nmap, c("scaffold", "length", "old_chr", "old_chrnum", "old_bin", "old_orientation", "new_orientation"))
 with(unique(nmap[, .(old_chrnum, old_chr)]), {x=old_chr; names(x)=paste(old_chrnum, old_chr); x}) -> n
 lapply(n, function(i) nmap[i, on="old_chr"]) -> ll
 openxlsx::write.xlsx(ll, file)
}

#read edited Hi-C 
read_hic_map <- function(rds, file){
 setnames(data.table(do.call(rbind, strsplit(readxl::excel_sheets(file), " ")))[, 1], "chrnum")[, idx := 1:.N][]  -> n
 n[chrnum == "NA", chrnum := NA]
 n[, chrnum := as.integer(n$chrnum)]
 readRDS(rds) -> nmap
 nmap$hic_map[, consensus_chr := NULL]
 nmap$hic_map[, consensus_orientation := NULL]
 nmap$hic_map[, hic_bin := NULL]
 m <- rbindlist(lapply(1:nrow(n), function(i) {
  data.table(openxlsx::read.xlsx(file, sheet=i))[, chrnum := n[i, chrnum]][, nbin := 1:.N] -> x
  x[, .(scaffold, consensus_chr = chrnum, hic_bin = nbin, consensus_orientation = new_orientation)]
 }))
 m[is.na(consensus_chr), hic_bin := NA]
 m[nmap$hic_map, on="scaffold"] -> nmap$hic_map
 nmap
}
 
#wrapper function to generate all Hi-C matrix plots with one command   
hic_plots<-function(rds, assembly, noplot=F, prefix=NULL, nuc, cores=1, species, cov=T){
 print("YOURE USIGNG A DEBUG VERSION")
 cat(paste0("Reading Hi-C map RDS file: ", rds, "\n"))
 readRDS(rds) -> hic_map

 if(is.null(prefix)){
 prefix <- sub(".Rds", "", sub("^[0-9]+", "", rds))
 }
 cat(paste0("Output prefix: ", prefix, "\n"))

 cat("Convert Hi-C pairs to pseudomolecule coordinates.\n")
 add_psmol_fpairs(assembly=assembly, hic_map=hic_map, nucfile=nuc)->hic_map

 cat("Bin and normalize Hi-C link counts.\n")
 bin_hic_step(hic=hic_map$links, frags=hic_map$frags, binsize=1e6, chrlen=hic_map$chrlen, cores=cores)->hic_map$hic_1Mb
 normalize_cis(hic_map$hic_1Mb, ncores=cores, percentile=0, omit_smallest=1)->hic_map$hic_1Mb$norm

 if(noplot){
  return(hic_map) 
 }

 if(cov){
  cat("Calculate physical Hi-C coverage.\n")
  print(paste0("Using binsize ", format(1e3, scientific=T), " bp and binsize2 ", format(1e6, scientific=T), " bp for Hi-C coverage calculation."))
  hic_cov_psmol(hic_map=hic_map, binsize=1e3, binsize2=1e6, maxdist=1e6, cores=40) -> cov

  cat("Find Hi-C coverage drops.\n")
  print(names(cov))
  print(cov)
  print("validation:")
  print(cov[, .(chr, bin, r)])

  setorder(cov, r)[r < -2][, .(r=r[1]), key=.(chr, bin %/% 1e5 * 1e5)]->br0
  print("order setted")

  print(typeof(hic_map$agp$chr))
  print(typeof(br0$chr))
 
  hic_map$agp[, chr := as.character(chr)]
  br0[, chr := as.character(chr)]

  hic_map$agp[, .(chr, scaffold, scaffold_length=agp_end-agp_start+1, agp_start, agp_end, orientation, bin=agp_start)][br0, on=c("chr", "bin"), roll=T]->br 
  print("defined br")
  br[orientation == 1, br := bin - agp_start]
  br[orientation == -1, br := agp_end - bin]
  print ("checkpoint")

  print(typeof(hic_map$chrlen$chr))
  print(typeof(br$chr))

  hic_map$chrlen[, chr := as.character(chr)]
  br[, chr := as.character(chr)]

  hic_map$chrlen[,.(chr, length)][br, on="chr"]->br
  print("defined again br")
  br[, d := pmin(bin, length - bin)]
  br[, ds := pmin(br, scaffold_length - br)]
  copy(br) -> br1
  br1[ds >= 1e3 & d >= 1e5]->br
  print ("defined br again")

  print(typeof(hic_map$agp$chr))
  print(typeof(br$chr))
  
  hic_map$br <- br
  print("finished defining br")
 } else {
  cov <- NULL
 }

 cat("Calculate directionality bias.\n")

 print(typeof(hic_map$chrlen$chr))
 print(typeof(hic_map$hic_1Mb$norm$chr1))

 # transform to integer also chr1
 hic_map$chrlen[, chr := as.integer(chr)]
 hic_map$hic_1Mb$norm[, chr1 := as.integer(chr1)]

 find_inversions(hic_map=hic_map, links=hic_map$hic_1Mb$norm, species=species, cores=cores)->inv
 print("found inversions")
 print(inv)

 hic_map$inv <- inv

 print("inversions added to hic_map. Now plotting.")

 f <- paste0(prefix, "_intrachromosomal.pdf")
 cat(paste0("Plot intrachromosomal contact matrix: ", f, "\n"))

 big_hic_plot(hic_map=hic_map, cov=cov, inv=inv, species=species, file=f, breaks=br, cores=cores)

 print("pltted")

 # print("plotted in: ", f)

 f <- paste0(prefix, "_interchromosomal.png")
 cat(paste0("Plot interchromosomal contact matrix: ", f, "\n"))
 interchromosomal_matrix_plot(hic_map=hic_map, file=f, species=species)

 f <- paste0(prefix, "_export.Rds")
 cat(paste0("Export Hi-C map inspector file: ", f, "\n"))
 inspector_export(hic_map=hic_map, assembly=assembly, inv=inv, species=species, file=paste0(prefix, "_export.Rds"))

 hic_map
}

#Function that calls shell code to compile pseudomolecules
compile_psmol <- function(hic_map, assembly, fasta, output, agp_info, cores=1, chrUn=T,
 bedtools='/opt/Bio/bedtools/2.29.0/bin/bedtools',
 samtools='/opt/Bio/samtools/1.11/bin/samtools',
 tabix='/opt/Bio/bcftools/1.9/bin/tabix',
 bgzip='/opt/Bio/bcftools/1.9/bin/bgzip',
 parallel='/opt/Bio/parallel/20150222/bin/parallel',
 tritex='/filer-dg/agruppen/seq_shared/mascher/code_repositories/tritexassembly.bitbucket.io/'
 ){

 runcmd <- function(cmd){
  cat(paste0(cmd, "\n"))
  if(system(cmd) > 0) { suppressWarnings(sink()); stop()}
  cat("\n")
 }

 n50 <- paste0(tritex, "shell/n50")
 write <- paste0(tritex, "shell/write_psmol.zsh")
 check <- paste0(tritex, "shell/check_psmol.zsh")
 
 if(dir.exists(output)){
  stop(paste0("Output directory ", output, " exists."))
 }

 dir.create(output)

 deparse(substitute(assembly)) -> na
 deparse(substitute(hic_map)) -> nh
 write_agp_files(hic_map, assembly, output, name_a=na, name_h=nh) -> agp_info

 psmol <- paste0(output, "/", output, ".fasta")
 psmolbase <- sub(".fasta$", "", psmol)

 agpall <- agp_info$agp_bed
 bedall <- agp_info$psmol_bed
 agp <- sub(".bed", "_noUn.bed", agpall)
 bed <- sub(".bed", "_noUn.bed", bedall)

 sub(".bed$", "", agp_info$assembly_bed) -> base
 paste0(base, ".fasta") -> assembly
 paste0(base, "+gap.fasta") -> psmolinput

 paste("grep -v chrUn", agpall, ">", agp) -> cmd
 runcmd(cmd)

 paste("grep -v chrUn", bedall, ">", bed) -> cmd
 runcmd(cmd)

 cmd <-paste(bedtools, "getfasta -name -fi", fasta, " -bed", agp_info$assembly_bed, "-fo /dev/stdout | sed 's/::/ /' > ", assembly)
 runcmd(cmd)

 paste(samtools, "faidx", assembly) -> cmd
 runcmd(cmd)

 paste(n50, paste0(assembly, ".fai"), ">", paste0(base, ".n50")) -> cmd
 runcmd(cmd)

 paste("awk 'BEGIN {print \">gap\"; for(i=0;i<1000;i++) printf \"N\"; print \"\"}' | cat -", assembly, ">", psmolinput) -> cmd
 runcmd(cmd)

 paste(samtools, "faidx", psmolinput) -> cmd
 runcmd(cmd)

 cmd <- paste(write, "--bed", bed, "--pathbedtools", bedtools, "--tabix", tabix, "--parallel", parallel, "--bgzip", bgzip,
	       "--assembly", psmolinput, "--cores", cores, ">", psmol, "2>", paste0(psmolbase, ".err"))  
 runcmd(cmd)

 paste(samtools, "faidx", psmol) -> cmd
 runcmd(cmd)

 cmd <- paste(check,  "--bedtools", bedtools, "--tabix", tabix, "--samtools", samtools, "--bgzip", bgzip,
	      "--assembly", assembly, "--agp", agp, "--psmol", psmol, ">", 
	      paste0(psmolbase, "_check.out"), "2>",  
	      paste0(psmolbase, "_check.err"))
 runcmd(cmd)

 paste("cat", paste0(psmolbase, "_check.out"), "2>", paste0(psmolbase, "_check.err"), "| wc -l") -> cmd
 fread(cmd=cmd, head=F)$V1 -> wc

 if(wc > 0){
  stop(paste("Pseudomolecule validation failed. Check", paste0(psmolbase, "_check.out.")))
 }

 paste(n50, paste0(psmol, ".fai"), ">", paste0(psmolbase, ".n50")) -> cmd
 runcmd(cmd)

 paste0("Deleting ", psmolinput, ".\n")
 file.remove(psmolinput)

 if(chrUn){
  paste("grep chrUn", agpall, "| awk '$4 != \"gap\"' | cut -f 4 > ", paste0(psmolbase, "_chrUn.txt")) -> cmd
  runcmd(cmd)

  paste(samtools, "faidx -r", paste0(psmolbase, "_chrUn.txt"), assembly, ">", paste0(base, "_unanchored_contigs.fasta")) -> cmd
  runcmd(cmd)

  paste(samtools, "faidx", paste0(base, "_unanchored_contigs.fasta")) -> cmd
  runcmd(cmd)

  paste(n50, paste0(base, "_unanchored_contigs.fasta.fai"), ">", paste0(base, "_unanchored_contigs.n50")) -> cmd
  runcmd(cmd)

  cmd <- paste('cat', psmol, paste0(base, "_unanchored_contigs.fasta"), ">", paste0(psmolbase, "+unanchored_contigs.fasta"))
  runcmd(cmd)

  paste(samtools, "faidx", paste0(psmolbase, "+unanchored_contigs.fasta")) -> cmd
  runcmd(cmd)

  paste(n50, paste0(psmolbase, "+unanchored_contigs.fasta.fai"), ">", paste0(psmolbase, "+unanchored_contigs.n50")) -> cmd
  runcmd(cmd)
 }

 cat(paste0("Pseudomolecule FASTA: ", psmol, ".\n"))
 cat(paste0("Corrected assembly FASTA: ", assembly, ".\n"))
 cat("Done.\n")
}

#Function to write AGP/BED file corrected assemblies and pseudomolecules
write_agp_files <- function(hic_map, assembly, dir=".", nopsmol=F, name_a, name_h){
 prefix <- paste0(dir, "/", format(Sys.Date(), "%y%m%d"), "_")

 name_a -> na
 name_h -> nh
 paste0(prefix, na, "_", nh, "_psmol.bed") -> pbed
 paste0(prefix, na, "_", nh, "_AGP.bed") -> gbed
 paste0(prefix, na, ".bed") -> abed

 assembly$info[, .(orig_scaffold, orig_start = orig_start - 1, orig_end, scaffold, length)] -> oldnew
 write.table(sep="\t", col.names=F, row.names=F, quote=F, oldnew[, .(orig_scaffold, orig_start, orig_end, scaffold)],
	     file=abed)

 cat(paste0("Written corrected assembly AGP file: ", abed, ".\n\n"))
 
 if(nopsmol){
  return(list(assembly_bed=abed))
 }

 a <- hic_map

 bed <- a$agp[, .(scaffold, start=0, end=scaffold_length, agp_chr, agp_start, orientation=ifelse(orientation == 1, "+", "-"))] 
 bed[is.na(orientation), orientation := "+"]
 options(scipen=20)
 write.table(bed, col.names=F, row.names=F, quote=F, sep="\t", file=pbed)

 cat(paste0("Written pseudomolecule BED file: ", pbed, ".\n\n"))

 agpbed<-a$agp[, .(agp_chr, agp_start = agp_start-1, agp_end, scaffold, index, orientation=ifelse(orientation == 1, "+", "-"))] 
 agpbed[is.na(orientation), orientation := "+"]
 options(scipen=20)
 write.table(agpbed, col.names=F, quote=F, row.names=F, sep="\t", file=gbed)

 cat(paste0("Written pseudomolecule AGP file: ", gbed, ".\n\n"))

 list(assembly_bed=abed, psmol_bed=pbed, agp_bed=gbed)
}
