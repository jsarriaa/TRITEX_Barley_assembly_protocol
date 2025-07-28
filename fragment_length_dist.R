#!/opt/Bio/R/3.4.2/bin/Rscript --vanilla 

.libPaths("/filer-dg/agruppen/seq_shared/mascher/Rlib/current")
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(modeest))

args <- commandArgs(trailingOnly = TRUE)
out <- args[1]
cores <- args[2]
samples <- sort(args[-2:-1])

files<-paste0(samples, "/", samples, "_length_dist.txt")
names(files)<-sub("Sample_", "", samples)
mclapply(mc.cores=cores, files, function(f) {
 fread(f)->f
 setnames(f, c("count", "length"))
})->dist

files<-paste0(samples, "/", samples, "_length_dist_PE.txt")
names(files)<-sub("Sample_", "", samples)
mclapply(mc.cores=cores, files, function(f) {
 fread(f)->f
 setnames(f, c("count", "length"))
})->dist_pe

mlv_pe<-mclapply(mc.cores=cores, dist_pe, function(i) mlv(method='naive', do.call(c, apply(i,1 ,function(x) rep(x[2], x[1]))))$M)
mlv_re<-mclapply(mc.cores=cores, dist, function(i) mlv(method='naive', do.call(c, apply(i,1 ,function(x) rep(x[2], x[1]))))$M)

pdf(out, paper="a4r", width=30)
invisible(lapply(names(dist), function(i){
 par(mfrow=c(1,2))
 par(mar=c(5,6,3,1))
 dist[[i]][, plot(xlim=c(0,600), col=0, length, count/1e3, las=1, ylab="", main=paste(i, "(restriction site associated)"), xlab="fragment length")]
 abline(v=seq(0,600,50), lwd=2, col="gray")
 dist[[i]][, points(pch=20, length, count/1e3)]
 title(ylab="no. of fragments (x 1000)")
 abline(lwd=2, col=2, v=mlv_re[[i]])
 dist_pe[[i]][, plot(xlim=c(0,600), col=0, length, count/1e3, las=1, ylab="", main=paste(i, "(paired-end contamination)"), xlab="fragment length")]
 abline(v=seq(0,600,50), lwd=2, col="gray")
 dist_pe[[i]][, points(pch=20, length, count/1e3)]
 title(ylab="no. of fragments (x 1000)")
 abline(lwd=2, col=2, v=mlv_pe[[i]])
}))
invisible(dev.off())
