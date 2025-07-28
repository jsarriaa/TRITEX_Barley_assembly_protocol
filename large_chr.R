bin_hic_step<-function(hic, frags, binsize, step=NULL, chrlen, chrs=NULL, cores=1){
 if(is.null(chrs)){
  chrs <- unique(chrlen$chr)
 }

 if(is.null(step)){
  step <- binsize
 }

 options(scipen=20)
 chrlen[chr %in% chrs, .(bin=seq(0, length, step)), key=chr][, end := bin + binsize - 1][, start := bin] -> bins
 bins[, id := paste0(chr, ":", bin)]
 
 hic[chr1 %in% chrs & chr2 %in% chrs] -> z

 if(step == binsize){
  z[, bin1 := as.numeric(start1 %/% binsize * binsize)]
  z[, bin2 := as.numeric(start2 %/% binsize * binsize)]
  z[, .(nlinks=.N), keyby=.(chr1,bin1,chr2,bin2)]->z
  z[, id1 := paste(sep=":", chr1, bin1)]
  z[, id2 := paste(sep=":", chr2, bin2)]
  z->mat

  frags[chr %in% chrs]->f
  f[, bin := as.numeric(start %/% binsize * binsize)]
  f[, .(nfrags=.N, eff_length=sum(length), cov=weighted.mean(cov, length), 
	gc=(sum(nC)+sum(nG))/(sum(nA)+sum(nC)+sum(nG)+sum(nT))), keyby=.(chr, bin)]->f
  f[, id := paste(sep=":", chr, bin)]
 } else {
  copy(bins)->b
  b[, idx := 1:.N]
  b[, .(win=seq(bin, bin+binsize-1, step)),.(chr, bin)]->b

  z[chr1 == chr2]->z
  z[, win1:=as.numeric(start1 %/% step * step)]
  z[, win2:=as.numeric(start2 %/% step * step)]

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
  f[, win := as.numeric(start %/% step * step)]
  b[f, on=c("chr", "win"), nomatch=0, allow.cartesian=T]->f
  f[, .(nfrags=.N, eff_length=sum(length), cov=weighted.mean(cov, length), 
	gc=(sum(nC)+sum(nG))/(sum(nA)+sum(nC)+sum(nG)+sum(nT))), keyby=.(chr, bin)]->f
  f[, id :=paste(sep=":", chr, bin)]
 }

 list(bins=f[], mat=mat[])
}

add_psmol_fpairs<-function(assembly, hic_map, nucfile, map_10x=NULL, assembly_10x=NULL, cov=NULL){
 if(is.null(map_10x)){
  assembly$fpairs[, .(scaffold1, scaffold2, pos1, pos2)] -> z
 } else {
  assembly_10x$fpairs[, .(scaffold1, scaffold2, pos1, pos2)] -> z
 }
 z[, pos1 := as.numeric(pos1)]
 z[, pos2 := as.numeric(pos2)]
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
