library(igraph)

read_cov_links<-function(prefix, save=F){

 f <- paste0(prefix, '_cov.tsv')
 fread(f, col.names=c("seq", "len", "kc", "km"), head=F) -> cov
 cov[, seq := paste(seq)]

 f <- paste0(prefix, '_links.tsv')
 fread(f, head=F, col.names=c("seq1", "seq2"))->e
 e[, c("seq1", "seq2") := list(as.character(seq1), as.character(seq2))]
 graph.edgelist(as.matrix(e), directed=F) -> g
 g + vertex(cov[!seq %in% c(e$seq1, e$seq2)]$seq) -> g
 data.table(seq=V(g)$name, degree=degree(g))[cov, on="seq"] -> cov

 res <- list(cov=cov, graph=graph, prefix=prefix)

 if(save){
  saveRDS(res, paste0(prefix, '_graph_cov.Rds'))
 }

 res
}

plot_assembly<-function(assembly, file=NULL){
 if(is.null(file)){
  file <- paste0(assembly$prefix, ".pdf")
 }

 cov <- assembly$cov

 cov[, .(km=mean(km)), key=len]->z
 copy(cov)[, len2 := len][, .(l=sum(as.numeric(len2))), key=len]->len
 len[, cs := cumsum(l)]
 len[, s := sum(as.numeric(cov$len)) - c(0, cs[-.N])]
 cov[, .(d=mean(degree)), key=.(len=10**(log10(len) %/% 0.02 * 0.02))]->y

 pdf(file, height=5, width=7)
 layout(matrix(1:3, nrow=3), height=c(4,5,5))
 par(mar=c(0,4.1,1,1))
 z[, plot(pch=20, col=1, bty='l', xaxt='n', yaxt='n', len, km, xlab="", ylab="coverage", log='xy')]
 axis(2, c(1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000), las=1)
 len[, plot(lwd=3, col=1, type='l', bty='l', xaxt='n', yaxt='n', bty='l', len, s/1e9, xlab="", ylab="assembly size", log='x')]
 axis(2, 0:5, paste(0:5, "Gb"), las=1)
 par(mar=c(4.1,4.1,1,1))
 y[, plot(pch=20, len, log='x', d, bty='l', xlab="", las=1, xaxt='n', ylab='degree')]
 axis(1, c(500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 2e4),
         c("500 bp", "1 kb", "2 kb", "5 kb", "10 kb", "20 kb", "50 kb", "100 kb", "200 kb"))
 title(xlab="unitig length")
 dev.off()
}
