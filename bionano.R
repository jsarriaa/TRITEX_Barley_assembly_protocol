library(bit64)
library(zoo)

readXMap_2seq<-function(file, refmap, qrymap){
 fread(paste("grep -v '^#'", file), header=F)->xmap
 setnames(xmap, fread(header=F, paste("grep -m 1 '^#h'", file, "| tr ' \\t' '\\n' | tail -n +2"))$V1)
 xmap[,  tstrsplit(unlist(strsplit(gsub("\\(", "", Alignment), ")")), ","), by="XmapEntryID"]->aln
 setnames(aln, c("V1", "V2"), c("site_id_ref", "site_id_query"))
 aln[, site_id_query := as.integer(site_id_query)]
 aln[, site_id_ref := as.integer(site_id_ref)]
 setkey(aln, "XmapEntryID")
 setkey(xmap[, .(XmapEntryID, QryContigID, RefContigID, QryLen, RefLen, Confidence, Orientation)], XmapEntryID)[ setkey(aln, "XmapEntryID")]->aln
 setkey(qrymap[, .(QryContigID=CMapID, scaffold_qry=scaffold, site_id_query = site_id, position_query = position)], QryContigID, site_id_query)[setkey(aln, QryContigID, site_id_query)]->aln
 setkey(refmap[, .(RefContigID=CMapID, scaffold_ref=scaffold, site_id_ref = site_id, position_ref = position)], RefContigID, site_id_ref)[setkey(aln, RefContigID, site_id_ref)]->aln

 unique(refmap[, .(RefContigID=CMapID, scaffold_ref=scaffold)])[xmap, on="RefContigID"]->xmap
 unique(qrymap[, .(QryContigID=CMapID, scaffold_qry=scaffold)])[xmap, on="QryContigID"]->xmap
 xmap[, Alignment := NULL]
 xmap[, LabelChannel := NULL]
 xmap[, AlnLenRef := abs(RefEndPos - RefStartPos)]
 xmap[, AlnLenQry := abs(QryEndPos - QryStartPos)]
 
 list(xmap=xmap[], aln=aln[], refmap=refmap, qrymap=qrymap)
}

readXMap_seq<-function(file, refmap, qrymap){
 fread(cmd=paste("grep -v '^#'", file), header=F)->xmap
 setnames(xmap, fread(header=F, cmd=paste("grep -m 1 '^#h'", file, "| tr ' \\t' '\\n' | tail -n +2"))$V1)
 xmap[,  tstrsplit(unlist(strsplit(gsub("\\(", "", Alignment), ")")), ","), by="XmapEntryID"]->aln
 setnames(aln, c("V1", "V2"), c("site_id_ref", "site_id_query"))
 aln[, site_id_query := as.integer(site_id_query)]
 aln[, site_id_ref := as.integer(site_id_ref)]
 setkey(aln, "XmapEntryID")
 setkey(xmap[, .(XmapEntryID, QryContigID, RefContigID, QryLen, RefLen, Confidence, Orientation)], XmapEntryID)[ setkey(aln, "XmapEntryID")]->aln
 ref_prefix <- sub("_[0-9]+$", "_", head(n=1, refmap$opt_contig))
 qry_prefix <- sub("_[0-9]+$", "_", head(n=1, qrymap$opt_contig))
 aln[, QryContigID := sub("^", qry_prefix, QryContigID)]
 xmap[, QryContigID := sub("^", qry_prefix, QryContigID)]
 #aln[, RefContigID := sub("^", ref_prefix, RefContigID)]
 setkey(qrymap[, .(QryContigID=opt_contig, site_id_query = site_id, position_query = position)], QryContigID, site_id_query)[setkey(aln, QryContigID, site_id_query)]->aln
 setkey(refmap[, .(RefContigID=CMapID, scaffold, site_id_ref = site_id, position_ref = position)], RefContigID, site_id_ref)[setkey(aln, RefContigID, site_id_ref)]->aln

 unique(refmap[, .(RefContigID=CMapID, scaffold)])[xmap, on="RefContigID"]->xmap
 xmap[, Alignment := NULL]
 xmap[, LabelChannel := NULL]
 xmap[, AlnLenRef := abs(RefEndPos - RefStartPos)]
 xmap[, AlnLenQry := abs(QryEndPos - QryStartPos)]
 
 list(xmap=xmap[], aln=aln[], refmap=refmap, qrymap=qrymap)

}

plot_ref_aln<-function(aln, good, ref){
 k<-ref
 good[scaffold == k]->g
 unique(g$XmapEntryID)->zz
 g[, .(a=QryLen[1]+min(position_ref), RefLen = RefLen[1]), keyby=QryContigID][, max(c(a,RefLen))]->l
 par(mfrow=c(13, 1))
 par(oma=c(0,0,3,0))
 par(mar=c(1,0,0,0))
 n<<-0
 lapply(zz, function(i){
  n <<- n+1
  good[XmapEntryID == i]->z
  ori <- z$Orientation[1]
  qrylen <- z$QryLen[1]
  z[, index := .I]
  intr<-z[, c(min(site_id_ref), max(site_id_ref))]
  aln$refmap[scaffold == z$scaffold[1] & site_id %between% intr][, .(scaffold, site_id, position)]->mapr
  mapr[!z[, .(scaffold, site_id=site_id_ref)], on=c("scaffold", "site_id")]->mapr

  intq<-z[, c(min(site_id_query), max(site_id_query))]
  aln$qrymap[opt_contig == z$QryContigID[1] & site_id %between% intq][, .(opt_contig, site_id, position)]->mapq
  mapq[!z[, .(opt_contig=QryContigID, site_id=site_id_query)], on=c("opt_contig", "site_id")]->mapq

  y=0.1
  z[, plot(ylim=c(-3*y,1+y), col=0, 0, xlim=c(0, l), axes=F, xlab="", ylab="", xaxt='n', yaxt='n')]
  z[, min(position_ref)]->offset
  if(nrow(mapr) > 0){
   mapr[, lines(rep(position, 2), 1 + c(-y, +y), col=2), by=position]
  }
  if(nrow(mapq) > 0){
   if(ori == "+"){
    mapq[, lines(offset + rep(position, 2), 0 + c(-y, +y), col=2), by=position]
   } else {
    mapq[, lines(offset + qrylen - (rep(position, 2)), 0 + c(-y, +y), col=2), by=position]
   }
  }
  z[1][, rect(0, 1-y, RefLen, 1+y)]
  z[1][, rect(offset, 0-y, offset + QryLen, 0+y)]
  if(ori == "+"){
   z[, lines(c(offset + position_query, position_ref), c(0+y, 1-y)), by=index]
  } else {
   z[, lines(c(offset + qrylen - position_query, position_ref), c(0+y, 1-y)), by=index]
  }
  z[1][, text(xpd=NA, l/2, 0-3*y, paste(sep="", QryContigID, " (", round(QryLen/1e6, 1)," Mb), confidence = ", Confidence))]
  #if(n %% 13 == 1){
  # title(g[1][, paste(sep="", scaffold, " (", round(RefLen/1e6, 1)," Mb)")], outer=T)
  #}
 })
 n <- n+1
 if(n <= 13)
  lapply(n:13, function(i) plot(col=0, 0, axes=F, xlab="", ylab=""))
 title(g[1][, paste(sep="", scaffold, " (", round(RefLen/1e6, 1)," Mb)")], outer=T)
}


plot_qry_aln<-function(aln, good, qry, rows=1, ld=0.1, pages=T){
 k<-qry
 good[k, on="QryContigID"]->g
 g[order(position_query)][!duplicated(XmapEntryID)]$XmapEntryID->zz
# g[, .(a=RefLen[1]+min(position_query), QryLen = QryLen[1]), keyby=scaffold][, max(c(a,QryLen))]->l
 par(mfrow=c(rows, 1))
 par(oma=c(0,0,3,0))
 par(mar=c(3,0,0,0))
 n<<-0
 np <- ceiling(length(zz) / rows)
 pp <<- 0
 lapply(zz, function(i){
  n <<- n+1
  good[XmapEntryID == i]->z
  ori <- z$Orientation[1]
  qrylen <- z$QryLen[1]
  reflen <- z$RefLen[1]
  z[, index := .I]
  intr<-z[, c(min(site_id_ref), max(site_id_ref))]
  aln$refmap[scaffold == z$scaffold[1] & site_id %between% intr][, .(scaffold, site_id, position)]->mapr
  mapr[!z[, .(scaffold, site_id=site_id_ref)], on=c("scaffold", "site_id")]->mapr

  intq<-z[, c(min(site_id_query), max(site_id_query))]
  aln$qrymap[opt_contig == z$QryContigID[1] & site_id %between% intq][, .(opt_contig, site_id, position)]->mapq
  mapq[!z[, .(opt_contig=QryContigID, site_id=site_id_query)], on=c("opt_contig", "site_id")]->mapq

  y=0.1
  v <- qrylen * 0.1
  xl <- c(0-v, qrylen + v)
  z[, plot(ylim=c(-3*y,1+y), col=0, 0, xlim=xl, axes=F, xlab="", ylab="", xaxt='n', yaxt='n')]
  if(ori == "+"){
   z[, min(position_query) - min(position_ref) ]->offset
  } else {
   z[, min(position_query) - min(RefLen - position_ref) ]->offset
  }
  if(nrow(mapq) > 0){
   mapq[, lines(rep(position, 2), 1 + c(-y, +y), col=2), by=position]
  }
  if(nrow(mapr) > 0){
   if(ori == "+"){
    mapr[, lines(offset + rep(position, 2), 0 + c(-y, +y), col=2), by=position]
   } else {
    mapr[, lines(offset + reflen - (rep(position, 2)), 0 + c(-y, +y), col=2), by=position]
   }
  }
  z[1][, rect(0, 1-y, QryLen, 1+y)]
  z[1][, rect(offset, 0-y, offset + RefLen, 0+y)]
  if(ori == "+"){
   z[, lines(c(offset + position_ref, position_query), c(0+y, 1-y)), by=index]
  } else {
   z[, lines(c(offset + reflen - position_ref, position_query), c(0+y, 1-y)), by=index]
  }
  if(ori == "+"){
   minp <- offset + z[, min(position_ref)]
   maxp <- offset + z[, max(position_ref)]
  } else {
   minp <- offset + z[, min(reflen - position_ref)]
   maxp <- offset + z[, max(reflen - position_ref)]
  }
  #if(ori == "+"){
  # z[, text(minp, 0 - y - ld, paste0(round(min(position_ref/1e6),1), " Mb"))]
  # z[, text(maxp, 0 - y - ld, paste0(round(max(position_ref/1e6),1), " Mb"))]
  #}
  #if(ori == "-"){
  # z[, text(maxp, 0 - y - ld, paste0(round(min(position_ref/1e6),1), " Mb"))]
  # z[, text(minp, 0 - y - ld, paste0(round(max(position_ref/1e6),1), " Mb"))]
 # }

  #z[1][, text(xpd=NA, qrylen/2, 0-1.2*ld-(3)*y, paste(sep="", scaffold, " (", round(RefLen/1e6, 1)," Mb), confidence = ", Confidence, ", aln len = ", round((z[, max(position_ref)] - min(z[, position_ref]))/1e6, 1), 
  #  " Mb\nunaligned sites (ref): ", nrow(mapr), " / ", length(unique(z$site_id_ref)),
  #     ", unaligned sites (qry): ", nrow(mapq), " / ", length(unique(z$site_id_query))))]
  z[1][, text(xpd=NA, qrylen/2, 0-1.2*ld-(3)*y, paste(sep="", scaffold, " (", round(RefLen/1e6, 1), " Mb), ",
						      "start = ",  paste0(round(min(z$position_ref/1e6),1), " Mb"),
						      ", end = ",  paste0(round(max(z$position_ref/1e6),1), " Mb"),
						      ", orientation = ", ori, 
						      ", confidence = ", Confidence, ", aln len ref = ", round((z[, max(position_ref)] - min(z[, position_ref]))/1e6, 1), 
    " Mb\n",
    "query start = ", paste0(round(z[, min(position_query/1e6)],1), " Mb"),
    ", query end = ", paste0(round(z[, max(position_query/1e6)],1), " Mb"),
    ", aln len qry = ", round((z[, max(position_query)] - min(z[, position_query]))/1e6, 1),
    " Mb\nunaligned sites (ref): ", nrow(mapr), " / ", intr[2]-intr[1]+1,
       ", unaligned sites (qry): ", nrow(mapq), " / ", intq[2]-intq[1]+1))] 
  if(n %% rows == 1){
   pp <<- pp + 1
   if(pages){
    title(g[1][, paste(sep="", QryContigID, " (", round(QryLen/1e6, 1)," Mb), page ", pp, "/", np)], outer=T)
   } else {
    title(g[1][, paste(sep="", QryContigID, " (", round(QryLen/1e6, 1)," Mb)")], outer=T)
   }
  }
 })
 n <- n+1
 if(n <= rows)
  lapply(n:rows, function(i) plot(col=0, 0, axes=F, xlab="", ylab=""))
}

plot_ref_aln_v2<-function(aln, good, ref, rows=1){
 k<-ref
 good[k, on="scaffold"]->g
 unique(g$XmapEntryID)->zz
 par(mfrow=c(rows, 1))
 par(oma=c(0,3,3,3))
 par(mar=c(1,0,0,0))
 n<<-0
 np <- ceiling(length(zz) / rows)
 pp <<- 0
 lapply(zz, function(i){
  n <<- n+1
  good[XmapEntryID == i]->z
  ori <- z$Orientation[1]
  qrylen <- z$QryLen[1]
  reflen <- z$RefLen[1]
  z[, index := .I]
  intr<-z[, c(min(site_id_ref), max(site_id_ref))]
  aln$refmap[scaffold == z$scaffold[1] & site_id %between% intr][, .(scaffold, site_id, position)]->mapr
  mapr[!z[, .(scaffold, site_id=site_id_ref)], on=c("scaffold", "site_id")]->mapr

  intq<-z[, c(min(site_id_query), max(site_id_query))]
  aln$qrymap[opt_contig == z$QryContigID[1] & site_id %between% intq][, .(opt_contig, site_id, position)]->mapq
  mapq[!z[, .(opt_contig=QryContigID, site_id=site_id_query)], on=c("opt_contig", "site_id")]->mapq

  y=0.1
  plen <- z[, max(position_ref) - min(position_ref)]
  v <- plen * 0.1
  xl <- c(0-v, plen + v)
  z[, plot(ylim=c(-3*y,1+y), col=0, 0, xlim=xl, axes=F, xlab="", ylab="", xaxt='n', yaxt='n')]

  z[,  0 - min(position_ref) ]->offsetr

  if(ori == "+"){
   z[, 0 - min(position_query) ]->offset
  } else {
   z[, 0 - min(QryLen - position_query) ]->offset
  }
  if(nrow(mapr) > 0){
   mapr[, lines(rep(position + offsetr, 2), 1 + c(-y, +y), col=2), by=position]
  }
  if(nrow(mapq) > 0){
   if(ori == "+"){
    mapq[, lines(offset + rep(position, 2), 0 + c(-y, +y), col=2), by=position]
   } else {
    mapq[, lines(offset + qrylen - (rep(position, 2)), 0 + c(-y, +y), col=2), by=position]
   }
  }
  z[1][, rect(0 + offsetr, 1-y, RefLen + offsetr, 1+y)]
  z[1][, rect(offset, 0-y, offset + QryLen, 0+y)]
  if(ori == "+"){
   z[, lines(c(offset + position_query, offsetr + position_ref), c(0+y, 1-y)), by=index]
  } else {
   z[, lines(c(offset + qrylen - position_query, offsetr + position_ref), c(0+y, 1-y)), by=index]
  }
  if(ori == "+"){
   minp <- offset + z[, min(position_query)]
   maxp <- offset + z[, max(position_query)]
  } else {
   minp <- offset + z[, min(qrylen - position_query)]
   maxp <- offset + z[, max(qrylen - position_query)]
  }
  if(ori == "+"){
   z[, text(minp, 0 - y - 0.05, paste0(round(min(position_query/1e6),1), " Mb"))]
   z[, text(maxp, 0 - y - 0.05, paste0(round(max(position_query/1e6),1), " Mb"))]
  }
  if(ori == "-"){
   z[, text(maxp, 0 - y - 0.05, paste0(round(min(position_query/1e6),1), " Mb"))]
   z[, text(minp, 0 - y - 0.05, paste0(round(max(position_query/1e6),1), " Mb"))]
  }

  z[, abs(max(position_query) - min(position_query))] -> ll
  z[1][, text(xpd=NA, ll/2, 0-3*y, paste(sep="", QryContigID, " (", round(QryLen/1e6, 1)," Mb), confidence = ", Confidence))]

  minp <- offsetr + z[, min(position_ref)]
  maxp <- offsetr + z[, max(position_ref)]
  z[, text(xpd=NA, minp, 1 + y + 0.05, paste0(round(min(position_ref/1e6),1), " Mb"))]
  z[, text(xpd=NA, maxp, 1 + y + 0.05, paste0(round(max(position_ref/1e6),1), " Mb"))]

  if(n %% rows == 0){
   pp <<- pp + 1
   title(g[1][, paste(sep="", scaffold, " (", round(RefLen/1e6, 1)," Mb), page ", pp, "/", np)], outer=T)
  }
 })
 n <- n+1
 if(n <= rows)
  lapply(n:rows, function(i) plot(col=0, 0, axes=F, xlab="", ylab=""))
}

readCMap<-function(file, keyfile=NULL, key=NULL, prefix=NULL, rename=T){
 fread(cmd=paste("grep -v '^#'", file), header=F, select=1:9)->m
 setnames(m, c("CMapID", "scaffold_opt_length", "num_sites", "site_id", "label_channel", "position", "sd", "coverage", "occurrence"))
 if(!is.null(keyfile)){
  fread(keyfile, select=1:2)->key
  setnames(key, c("CMapID", "scaffold"))
  key[m, on="CMapID"]->m
 } else if(!is.null(key)){
  key[m, on="CMapID", nomatch=0]->m
 } else if(rename){
  setnames(m, "CMapID", "opt_contig")
  setnames(m, "scaffold_opt_length", "contig_length")
  m[, opt_contig := sub("^", paste0(prefix, "_"), opt_contig)]
 } else {
  m[, scaffold := sub("^",  paste0(prefix, "_"), CMapID)]
 }
 list(map=m, info=unique(m[, 1:4, with=F]))
}

read_optical_alignment<-function(xmap, key, prefix="optical_contig"){
 tempfile()->tmp
 system(paste("grep -v '^#'", xmap, ">", tmp))
 aln<-fread(tmp, select=1:10, sep="\t")
 unlink(tmp)
 setnames(aln, c("XmapEntryID", "QryContigID", "opt_contig", "fragment_start", "fragment_end", "opt_contig_start", "opt_contig_end", "orientation", "confidence", "hit_enum"))
 aln[, opt_contig:=sub("^", paste(sep="", prefix, "_"), opt_contig)]
 key<-fread(key, select=1:2)
 setnames(key, c("QryContigID", "fragment"))
 key[, bac:=toupper(sub("_[^:]+$", "", sub(":[0-9]+-[0-9]+$", "", sub("^[^:]+:", "", fragment))))]
 setkey(key, "QryContigID")
 setkey(aln, "QryContigID")
 key[aln][, c("QryContigID", "XmapEntryID"):=list(NULL, NULL)]
}

addMorexAln<-function(file, opt_map, morex_map){
 opt_map$morex_aln<-list()
 read_xmap(file, morex_map$map, opt_map$map)->opt_map$morex_aln$xmap
 opt_map$morex_aln$xmap[, .(queryStart=min(position_query),  
          queryEnd=max(position_query),  
          refStart=min(position_ref),
          refEnd=max(position_ref),
	  orientation=Orientation[1],
	  confidence=Confidence[1],
          queryLen=QryLen[1],
          refLen=RefLen[1]),
	  keyby=.(opt_morex_contig=RefContigID, opt_query_contig=QryContigID, XmapEntryID)]->opt_map$morex_aln$summary
 opt_map$morex_aln$summary[, alnLen := refEnd - refStart]
 opt_map$morex_aln$summary[, refAlnFraction := alnLen/refLen]
 opt_map$morex_aln$summary[, queryAlnFraction := alnLen/queryLen]
 opt_map$morex_aln$summary[, .(refStart = min(refStart),
          refEnd = min(refEnd),
          queryStart = min(queryStart),
          queryEnd = min(queryEnd),
          queryLen=queryLen[1],
          refLen=refLen[1],
	  alnLen=sum(alnLen),
	  nentries = .N,
	  orientation = as.character(ifelse(length(unique(orientation)) == 1, orientation[1], NA))),
	  keyby=.(opt_morex_contig, opt_query_contig)]->opt_map$morex_aln$info
 opt_map$morex_aln$info[, refAlnFraction := alnLen/refLen]
 opt_map$morex_aln$info[, queryAlnFraction := alnLen/queryLen]
 opt_map
}


addAGPAlignment<-function(opt_map, agp_data, xmap, key, prefix, min_conf=25){
 read_optical_alignment(xmap, key, prefix)->opt_map$bac_aln
 make_opt_aln_agp(opt_map$bac_aln, agp_data, opt_map$info, min_conf)->opt_map$agp_aln
 opt_map
}

readOpticalMap<-function(file, prefix, xmap, molinfo, header=0, cores=20){
 print("reading cmap")
 read_optical_map(file, prefix)->opt_map
 opt_map$map[, .(site_id=site_id[-.N], position=position[-.N], dist=position[-1] - position[-.N]), keyby=opt_contig]->site_dist
 setkey(site_dist[, .(opt_contig, site_id, dist)], opt_contig, site_id)[setkey(opt_map$map, opt_contig, site_id)]->opt_map$map

 cat("reading molecule info\n")
 readMolInfo(molinfo)->molinfo
 opt_map$molinfo <- molinfo

 cat("reading molecule alignment\n")
 opt_map$molaln<-readMolAln(cmd=xmap, refmap=opt_map$info, molinfo=molinfo, header=header)
 molaln<-opt_map$molaln
 molaln[, .(maxMolLen = max(QryLen)), keyby=.(opt_contig = RefContigID)][opt_map$info]->opt_map$info

 cat("calculating coverage statistics\n")
 rbindlist(mclapply(mc.cores=cores, unique(molaln$RefContigID), function(ref) {
  z<-molaln[RefContigID == ref]
  reflen<-unique(z$RefLen) %/% 1e3
  z[, .(QryContigID, RefStartPos = RefStartPos %/% 1e3, RefEndPos = RefEndPos %/% 1e3, AlnLen = AlnLen %/% 1e3)]->z
  z[order(-AlnLen)]->z
  data.table(key="pos", pos=1:(reflen+1), cov=0)->cov
  invisible(z[, cov[RefStartPos:RefEndPos + 1, cov := cov + 1],  by=QryContigID])
  cov[, .(opt_contig = ref, cov_mean=mean(cov), cov_sd=sd(cov), cov_max=max(cov), cov_min=min(cov[pos >= 100 & pos <= max(pos) - 100]), cov_median=median(cov), cov_mad=mad(cov))]
 }))->covstat
 setkey(covstat, opt_contig)[setkey(opt_map$info, opt_contig)]->opt_map$info

 cat("calculating repeat statistics\n")
 setkey(opt_map$info[, .(opt_contig, num_sites)], opt_contig)[setkey(site_dist, opt_contig)]->site_dist
 site_dist[, c("distMean", "distSD", "distMedian", "distMAD") := list(
    rollapply(fill=NA, dist, min(num_sites[1], 10), mean),
    rollapply(fill=NA, dist, min(num_sites[1], 10), sd),
    rollapply(fill=NA, dist, min(num_sites[1], 10), median),
    rollapply(fill=NA, dist, min(num_sites[1], 10), mad)),
   by=opt_contig]->molecule_div
 molecule_div[, repeatIndex :=  1 - distMAD/distMedian]
 setkey(molecule_div[, .(opt_contig, site_id, repeatIndex, distMean, distSD, distMedian, distMAD)], opt_contig, site_id)[setkey(opt_map$map, opt_contig, site_id)]->opt_map$map
 opt_map$map[, repetitive := rollapply(fill=NA, repeatIndex, min(num_sites[1], 10), function(x) all(x > 0.85)), by=opt_contig]
 opt_map$map[, .(repetitive = any(repetitive)), keyby=opt_contig][setkey(opt_map$info, opt_contig)]->opt_map$info
 opt_map
}

make_opt_aln_fpc<-function(opt_aln, fpc_info, fpc_table, info, min_conf){
 mem <- fpc_table
 cl <- fpc_info
 setkey(cl, fpc)

 x<-opt_aln
 setkey(x, bac)
 mem[!is.na(fpc), .(bac, fpc, cb_start, cb_end)]->y
 setkey(y, bac)

 y[x][!is.na(fpc)][order(-confidence)][!duplicated(fragment)]->x
 x[, .(confidence=sum(confidence)), keyby=.(opt_contig, fpc)][confidence >= min_conf]->opt_clu

 setkey(x, fpc, opt_contig)
 setkey(opt_clu, fpc, opt_contig)

 x[opt_clu[, .(fpc, opt_contig)]]->aln

 opt_clu[aln[, .(nbac=length(unique(bac)),
	    opt_contig_start=min(opt_contig_start), opt_contig_end=max(opt_contig_end),
	    cb_start=min(cb_start), cb_end=max(cb_end),
	    cor=suppressWarnings(cor(method='s', (opt_contig_start+opt_contig_end)/2, (cb_start+cb_end)/2))),
     keyby=.(fpc, opt_contig)]]->m
 m[, opt_aln_length:= opt_contig_end - opt_contig_start]
 m[, fpc_aln_length:= cb_end - cb_start]
 z<-info[, .(opt_contig, contig_length)]
 setnames(z, "contig_length", "opt_contig_length")
 setkey(z, "opt_contig")
 setkey(m, "opt_contig")
 z[m]->m
 m<-cl[, .(fpc, chr, cM, fpc_length = cb_length)][setkey(m, fpc)]
 list(info=m, aln=aln)
}

readInsilico<-function(file, key){
 fread(file)->m
 setnames(m, c("id", "contig_length", "num_sites", "site_id", "label_channel", "position", "sd", "coverage", "occurrence"))
 fread(key, select=1:2)->k
 setnames(k, c("id", "opt_contig"))
 setkey(k, id)[setkey(m, id)]->m
 list(map=m[, id := NULL], info=unique(m[, .(opt_contig, contig_length, num_sites)]))
}

plot_fpc<-function(input, i){
 input$fpc_table->z
 z[fpc == i]->z
 reflen<-max(ceiling(z$cb_end*1.24))
 z[, .(bac, RefStartPos = ceiling(1.24*cb_start), RefEndPos = ceiling(1.24*cb_end), AlnLen = ceiling(1.24*cb_length))]->z
 mem<-agp_data$oriented_info$mem
 setkey(mem[, .(cluster, bac, bac_bin=bin)], bac)[setkey(z, bac)]->z
 z[, lwd := ifelse(is.na(cluster), 2, 4)]
 z[, lib:=sub("[0-9]+[A-Z][0-9]{2}$", "", sub("HVVMRX", "", bac))]
 z[order(-lwd, -AlnLen)]->z
 data.table(key="pos", pos=1:(reflen+1), cov=0)->cov
 invisible(z[, cov[RefStartPos:RefEndPos + 1, cov := cov + 1],  by=bac])
 matrix(T, nrow=reflen+1, ncol=2*max(cov$cov))->occ
 setkey(setnames(z[, {
     head(which(apply(occ[RefStartPos:RefEndPos + 1,], 2, all)), n=1)->x
     occ[RefStartPos:RefEndPos + 1, x] <<- F
     x
 }, by=bac],  "V1", "y"), bac)[setkey(z, bac)]->z
 cov->baccov

 z->bacmap

 agp_data$agp[!is.na(cluster), .(posMin = min(agp_pos), posMax = max(agp_pos)), keyby=.(cluster, bac_bin)]->w
 setkey(w, cluster, bac_bin)[setkey(z[!is.na(cluster)], cluster, bac_bin)]->w
 w[order(RefStartPos)]->w
 w[, y := .GRP, by=cluster]
 suppressWarnings(agp_data$oriented_info$hic_map[, .(cluster, chr, cM, hic_bin)][w[, .(n=.N, agp_start=min(na.omit(posMin)), agp_end=max(na.omit(posMax)), y=unique(y)), keyby=cluster]])->clinfo

 layout(matrix(1:3, ncol=1), height=c(1,1,max(clinfo$y)*0.2))
 par(mar=c(0,5,3,1))
 baccov[, plot(type='l', lwd=3, pos, cov, axes=F, xlab='', las=1, ylab='coverage')]
 baccov[, abline(col='red', h=median(cov), lwd=3)]
 axis(2, las=1)
 title(i)

 par(mar=c(0,5,0,1))
 colors=c("orange", "green", "blue", "cyan", "violet", "black", "red")
 names(colors)<-c("83KHA", "ALLEA", "ALLHA", "ALLHB", "ALLHC", "ALLMA", "ALLRA")

 bacmap[, plot(xlim=c(0, reflen), ylim=c(max(30, max(y)), 0), col=0, 0, axes=F, xlab='', ylab='')]
 #bacmap[, lines(col=colors[lib], lwd=lwd, c(RefStartPos, RefEndPos), rep(y,2)), by=bac]
 bacmap[, lines(col=ifelse(is.na(cluster), 1,2), lwd=2, c(RefStartPos, RefEndPos), rep(y,2)), by=bac]
 #mtext(side=2, adj=0, names(colors), las=1, line=3, at=seq(3,max(30, max(bacmap$y)), length.out=7), xpd=NA, col=colors)

 par(mar=c(4.5,5,0,1))
 clinfo[, plot(xlim=c(0, reflen), ylim=c(max(y)+0.1, 0), col=0, 0, axes=F, xlab='', ylab='')]
 w[, lines(c(RefStartPos, RefEndPos), rep(y, 2), lwd=5, col="darkgray"), by=bac]
 clinfo[, 
  text(xpd=NA, 0, y - 0.5, adj=0, 
       paste(sep="", n, ifelse(n > 1, " BACs on ", " BAC on "), cluster, ifelse(is.na(chr), "", paste(sep="", ", ", chr, "H, ", round(cM,1), " cM, HiC: ", hic_bin, 
	      ", AGP length: ", 
	      ifelse(agp_start == Inf | agp_end == -Inf, "short", paste(sep="", round((agp_end - agp_start) / 1e6, 2), " Mb (", round(agp_start/1e6, 2), " - ", round(agp_end/1e6, 2), " Mb)"))))))
 , by=cluster]
 axis(1)
 title(xlab=paste("position in ", i, "(kb)"))

}

read_optical_alignment<-function(xmap, key, prefix="optical_contig"){
 tempfile()->tmp
 system(paste("grep -v '^#'", xmap, ">", tmp))
 aln<-fread(tmp, select=1:10, sep="\t")
 unlink(tmp)
 setnames(aln, c("XmapEntryID", "QryContigID", "opt_contig", "fragment_start", "fragment_end", "opt_contig_start", "opt_contig_end", "orientation", "confidence", "hit_enum"))
 aln[, opt_contig:=sub("^", paste(sep="", prefix, "_"), opt_contig)]
 key<-fread(key, select=1:2)
 setnames(key, c("QryContigID", "fragment"))
 key[, bac:=toupper(sub("_[^:]+$", "", sub(":[0-9]+-[0-9]+$", "", sub("^[^:]+:", "", fragment))))]
 setkey(key, "QryContigID")
 setkey(aln, "QryContigID")
 key[aln][, c("QryContigID", "XmapEntryID"):=list(NULL, NULL)]
}

make_opt_aln_agp<-function(opt_aln, agp_data, info, min_conf){
 mem <- agp_data$oriented_info$mem
 cl <- agp_data$oriented_info$hic_map
 setkey(cl, "cluster")

 agp_data$agp[!is.na(cluster), .(posMin = min(agp_pos), posMax = max(agp_pos)), keyby=.(cluster, bin=bac_bin)]->w
 w[setkey(mem, cluster, bin)]->mem

 x<-opt_aln
 setkey(x, "bac")
 mem[!is.na(cluster), .(bac, posMin, posMax, cluster, bin)]->y
 setkey(y, "bac")

 y[x][!is.na(cluster)][order(-confidence)][!duplicated(fragment)]->x
 x[, .(confidence=sum(confidence)), keyby=.(opt_contig, cluster)][confidence >= min_conf]->opt_clu

 setkey(x, cluster, opt_contig)
 setkey(opt_clu, cluster, opt_contig)

 x[opt_clu[, .(cluster, opt_contig)]]->aln

 opt_clu[aln[, .(nbac=length(unique(bac)),
	    opt_contig_start=min(opt_contig_start), opt_contig_end=max(opt_contig_end),
	    agp_start=min(na.omit(posMin)), agp_end=max(na.omit(posMax)),
	    bac_bin_start=min(na.omit(bin)), bac_bin_end=max(na.omit(bin)),
	    cor=suppressWarnings(cor(method='s', (opt_contig_start+opt_contig_end)/2, bin))),
     keyby=.(cluster, opt_contig)]]->m
 m[, opt_aln_length:= opt_contig_end - opt_contig_start]
 m[, agp_aln_length:= agp_end - agp_start]
 z<-info[, .(opt_contig, contig_length)]
 setnames(z, "contig_length", "opt_contig_length")
 setkey(z, "opt_contig")
 setkey(m, "opt_contig")
 z[m]->m
 m<-cl[, .(cluster, chr, cM, hic_bin)][setkey(m, cluster)]
 list(info=m, aln=aln)
}

read_molecule_cmap<-function(file, prefix){
 fread(file)->m
 setnames(m, c("opt_contig", "contig_length", "num_sites", "site_id", "label_channel", "position", "sd", "coverage", "occurrence"))
 m[, opt_contig:= sub("^", paste("opt_", prefix ,"_contig_", sep=""), opt_contig)]
 list(map=m, info=unique(m[, 1:3, with=F]))
}

plot_molaln<-function(opt_map, ref, log='', v=NULL, loess=F){

 molaln<-opt_map$molaln
 refmap<-opt_map
 opt_aln_agp<-opt_map$agp_aln

 refmap$map[opt_contig == ref]->site_dist

 z<-molaln[RefContigID == ref]
 reflen<-unique(z$RefLen) %/% 1e3
 z[, .(QryContigID, RefStartPos = RefStartPos %/% 1e3, RefEndPos = RefEndPos %/% 1e3, AlnLen = AlnLen %/% 1e3)]->z
 z[order(-AlnLen)]->z
 data.table(key="pos", pos=1:(reflen+1), cov=0)->cov
 invisible(z[, cov[RefStartPos:RefEndPos + 1, cov := cov + 1],  by=QryContigID])
 matrix(T, nrow=reflen+2, ncol=2*max(cov$cov))->occ
 setkey(setnames(z[, {
     head(which(apply(occ[RefStartPos:RefEndPos + 1,], 2, all)), n=1)->x
     occ[RefStartPos:RefEndPos + 1, x] <<- F
     x
 }, by=QryContigID],  "V1", "y"), QryContigID)[setkey(z, QryContigID)]->z
 z[, color := "black"]
 z[AlnLen >= 500, color := "blue"]
 z[AlnLen >= 1000, color := "red"]
 z[AlnLen >= 2000, color := "green"]
 z->molmap
 cov->molcov

 opt_aln_agp$aln[opt_contig == ref]->z
 if(nrow(z) > 0) {
  z[, .(fragment, cluster, RefStartPos = opt_contig_start %/% 1e3, RefEndPos = opt_contig_end %/% 1e3, AlnLen = (opt_contig_end - opt_contig_start) %/% 1e3)]->z
  z[order(-AlnLen)]->z
  reflen2<-max(z$RefEndPos)
  data.table(key="pos", pos=1:(reflen2+1), cov=0)->cov
  invisible(z[, cov[RefStartPos:RefEndPos + 1, cov := cov + 1],  by=fragment])
  matrix(T, nrow=reflen2+1, ncol=2*max(cov$cov))->occ
  setkey(setnames(z[, {
      head(which(apply(occ[RefStartPos:RefEndPos + 1,], 2, all)), n=1)->x
      occ[RefStartPos:RefEndPos + 1, x] <<- F
      x
  }, by=fragment],  "V1", "y"), fragment)[setkey(z, fragment)]->z
  z[, color := "black"]
  z->bacmap

  opt_aln_agp$info[opt_contig == ref]->z
  z[, .(cluster, chr, cM, hic_bin, opt_contig_length, agp_aln_length, agp_start, agp_end, RefStartPos = opt_contig_start %/% 1e3, RefEndPos = opt_contig_end %/% 1e3, AlnLen = (opt_aln_length) %/% 1e3)]->z
  z[order(RefStartPos)][, y := 1:.N]->clmap
  setkey(clmap[, .(cluster, yy=y)], cluster)[setkey(bacmap, cluster)]->bacmap
 } else {
  clmap <- NULL
  bacmap <- NULL
 }

 #if(!is.null(bacmap)){
 # layout(matrix(1:6), height=c(2,1.5,1.5,2,0.25*max(clmap$y),1))
 #} else {
 # layout(matrix(1:4), height=c(2,1.5,1.5,2))
 #}
 if(!is.null(bacmap)){
  layout(matrix(1:5), height=c(2,1.5,2,0.25*max(clmap$y),1))
 } else {
  layout(matrix(1:3), height=c(2,1.5,2))
 }

 par(mar=c(1,5,3,1))
 site_dist[, plot(log='y', position, rep(1, .N), col=0, ylim=c(1, 100), ylab="", axes=F, xlab="")]
 site_dist[, nextpos := c(position[2:.N], NA)]
 if(any(site_dist[!is.na(repetitive)]$repetitive)){
  site_dist[repetitive == T, rect(position, 0.1, nextpos, 1000, border=NA, col="#FF000033")]
 }
 site_dist[, points(pch=20, position, pmin(pmax(dist/1e3, 1), 100))]
 if(loess)
  site_dist[, lines(lwd=3, col='red', lowess(position, f=0.05, pmin(pmax(dist/1e3, 1), 1000)))]
 title(ylab="label distance (kb)")
 if(!is.null(bacmap)){
  title(ref)
 } else {
  title(paste(ref, "(no aligned clusters)"))
 }
 axis(2, las=1, xpd=NA)
 if(!is.null(v))
  abline(v=v, lwd=3, col="cyan")

 #par(mar=c(1,5,0,1))
 #site_dist[, plot(pch=20, position, ylab="repeat score", repeatIndex, axes=F, xlab="")]
 #axis(2, las=1)
 #if(!is.null(v))
 # abline(v=v, lwd=3, col="cyan")

 par(mar=c(0,5,0,1))
 molcov[, plot(type='l', lwd=3, pos, cov, axes=F, xlab='', log=ifelse(log == T, 'y', ''), las=1, ylab='coverage')]
 axis(2, las=1)
 if(!is.null(v))
  abline(v=v/1e3, lwd=3, col="cyan")

 if(is.null(bacmap)){
  par(mar=c(4,5,0,1))
 } else {
  par(mar=c(0,5,0,1))
 }
 molmap[, plot(xlim=c(0, reflen), ylim=c(max(y), 0), col=0, 0, axes=F, xlab='', ylab='')]
 molmap[, lines(col=color, c(RefStartPos, RefEndPos), rep(y,2)), by=QryContigID]
 title(ylab="aligned molecules")
 if(!is.null(v))
  abline(v=v/1e3, lwd=3, col="cyan")
 if(is.null(bacmap)){
  axis(1)
 }

 if(!is.null(bacmap)) {
  par(mar=c(0,5,0,1))
  clmap[, plot(xlim=c(0, reflen2), ylim=c(max(y)+0.1, 0), col=0, 0, axes=F, xlab='', ylab='')]
  title(xlab=paste("Position in", ref, "(kb)"))
  bacmap[, lines(c(RefStartPos, RefEndPos), rep(yy, 2), lwd=5, col="darkgray"), by=fragment]
  clmap[, 
   text(xpd=NA, 0, y - 0.5, adj=0, 
	paste(sep="", cluster, ", ", chr, "H, ", round(cM,1), " cM, HiC: ", hic_bin, 
	      ", aligned AGP length: ", round(agp_aln_length / 1e6, 1), " Mb (", round(agp_start/1e6, 1), " - ", round(agp_end/1e6, 1), " Mb)"))
  , by=cluster]

  par(mar=c(4,5,0,1))
  bacmap[, plot(xlim=c(0, reflen2), ylim=c(max(y)+1, 0), col=0, 0, axes=F, xlab='', ylab='')]
  title(xlab=paste("Position in", ref, "(kb)"))
  bacmap[, lines(col=color, lwd=5, c(RefStartPos, RefEndPos), rep(y,2)), by=fragment]
  title(ylab="aligned BACs", xpd=NA)
  axis(1)
  if(!is.null(v))
   abline(v=v/1e3, lwd=3, col="cyan")
 }
}

readMolInfo<-function(file){
 molinfo<-fread(file, sep="\t")
 molinfo[, V3 := as.character(V3)]
 setnames(molinfo, c("genotype", "labelchannel", "moleculeid", "length", "avgintensity", "snr", "nlabels", 
		"origmoleculeid", "scannumber", "scandirection", "chipid", "flowcell", "runid", 
		"globalscannumber"))
 molinfo
}

readMolAln<-function(cmd, refmap, molinfo, header=0){
 fread(cmd, sep="\t")->molaln
 molaln[, V2 := as.character(V2)]
 if(header == 0){
  setnames(molaln, c("XmapEntryID", "QryContigID", "RefContigID", 
		    "QryStartPos", "QryEndPos", "RefStartPos", "RefEndPos",
		    "Orientation",  "Confidence", "HitEnum"))
  setkey(molinfo[, .(QryContigID = moleculeid, QryLen = length)], QryContigID)[setkey(molaln, QryContigID)]->molaln
 } else {
  setnames(molaln, c("XmapEntryID","QryContigID", "RefContigID", "QryStartPos", "QryEndPos", 
    "RefStartPos", "RefEndPos", "Orientation", "Confidence", "HitEnum", "QryLen", "RefLen", "LabelChannel", "Alignment"))
 }
 molaln[, AlnLen := abs(QryEndPos - QryStartPos)]
 molaln[, AlnFraction := AlnLen / QryLen]
 ref_prefix <- sub("_[0-9]+$", "_", head(n=1, refmap$opt_contig))
 molaln[, RefContigID := sub("^", ref_prefix, RefContigID)]
 if(header == 0){
  setkey(refmap[, .(RefContigID = opt_contig, RefLen = contig_length)], RefContigID)[setkey(molaln, RefContigID)]->molaln
 }
 molaln
}

n50<-function(l){
 l[order(l)] -> l
 l[head(which(cumsum(l) >= 0.5 * sum(l)), n=1)]
}

read_smap<-function(file, refmap, qrymap){
 fread(cmd=paste("grep -v '^#'", file), header=F)->smap
 setnames(smap, fread(header=F, paste("grep -m 1 '^#h'", file, "| tr ' \\t' '\\n' | tail -n +2"))$V1)->smap
 setnames(smap, c("RefcontigID1", "RefcontigID2"), c("RefContigID1", "RefContigID2"))

 ref_prefix <- sub("_[0-9]+$", "_", head(n=1, refmap$opt_contig))
 qry_prefix <- sub("_[0-9]+$", "_", head(n=1, qrymap$opt_contig))
 smap[, RefContigID1:= sub("^", ref_prefix, RefContigID1)]
 smap[, RefContigID2:= sub("^", ref_prefix, RefContigID2)]
 smap[, QryContigID:= sub("^", qry_prefix, QryContigID)][]
}

plot_xmap<-function(xmap, i, j=NULL){
 aln<-xmap
 aln[RefContigID == i]->z
 if(!is.null(j)){
  z[QryContigID %in% j]->z
 }
 unique(z$RefLen)->reflen
 z[order(position_ref), off_ref := position_ref[1], by=QryContigID]
 z[, query_orientation := Orientation[1] == "+", by=QryContigID]
 z[, off_qry := ifelse(query_orientation, position_query[1], position_query[.N]), by=QryContigID]
 z[, qry_start := off_ref[1]-off_qry[1], by=QryContigID]
 z[, qry_end := off_ref[1]-off_qry[1]+QryLen[1], by=QryContigID]
 z[, y := -.GRP+1, by=QryContigID]

 par(mar=c(1,11,3,1))
 z[, plot(c(min(0, qry_start), max(reflen, qry_end)), c(min(y)-0.1,1.1), col=0, xlab="", ylab="", yaxt='n', axes=F, xaxt='n')]
 rect(0, 1, reflen, 1.1)
 #text(reflen/2, 1.2, i, xpd=NA)
 axis(2, 1-0.05, i, las=1, tick=F)
 z[!duplicated(QryContigID), axis(2, y-0.05, QryContigID, las=1, tick=F, line=0)]
 z[, line := 1:.N]
 z[, lines(c(off_ref+ position_query - off_qry, position_ref), c(y,1)), by=line]
 z[, rect(qry_start[1], y[1]-0.1, qry_end[1], y[1]) , by=QryContigID]
 #axis(3, c(0, reflen))
 axis(3)
}

read_xmap<-function(file, refmap, qrymap){
 fread(cmd=paste("grep -v '^#'", file), header=F)->xmap
 setnames(xmap, fread(header=F, cmd=paste("grep -m 1 '^#h'", file, "| tr ' \\t' '\\n' | tail -n +2"))$V1)
 xmap[,  tstrsplit(unlist(strsplit(gsub("\\(", "", Alignment), ")")), ","), by="XmapEntryID"]->aln
 setnames(aln, c("V1", "V2"), c("site_id_ref", "site_id_query"))
 aln[, site_id_query := as.integer(site_id_query)]
 aln[, site_id_ref := as.integer(site_id_ref)]
 setkey(aln, "XmapEntryID")
 setkey(xmap[, .(XmapEntryID, QryContigID, RefContigID, QryLen, RefLen, Confidence, Orientation)], XmapEntryID)[ setkey(aln, "XmapEntryID")]->aln
 ref_prefix <- sub("_[0-9]+$", "_", head(n=1, refmap$opt_contig))
 qry_prefix <- sub("_[0-9]+$", "_", head(n=1, qrymap$opt_contig))
 aln[, QryContigID := sub("^", qry_prefix, QryContigID)]
 aln[, RefContigID := sub("^", ref_prefix, RefContigID)]
 setkey(qrymap[, .(QryContigID=opt_contig, site_id_query = site_id, position_query = position)], QryContigID, site_id_query)[setkey(aln, QryContigID, site_id_query)]->aln
 setkey(refmap[, .(RefContigID=opt_contig, site_id_ref = site_id, position_ref = position)], RefContigID, site_id_ref)[setkey(aln, RefContigID, site_id_ref)]
}

read_optical_map<-function(file, prefix){
 fread(cmd=paste("grep -v ^#", file))->m
 as.character(fread(cmd=paste("grep -m 1 '^#h'", file), head=F))->header
 setnames(m, c("opt_contig", "contig_length", "num_sites", "site_id", "label_channel", "position", "sd", "coverage", "occurrence", header[10:length(header)]))
 m[, opt_contig:= sub("^", paste("opt_", prefix ,"_contig_", sep=""), opt_contig)]
 list(map=m, info=unique(m[, 1:3, with=F]))
}

make_optical_joins<-function(opt_aln, bacgraph, info, min_conf){
 mem <- bacgraph$membership
 cl <- bacgraph$cluster_info[, list(cluster, chr, cM)]
 setkey(cl, "cluster")

 x<-opt_aln
 setkey(x, "bac")
 mem[, list(bac, cluster, bin)]->y
 setkey(y, "bac")

 y[x][!is.na(cluster)]->x
 setorder(x, -confidence)
 x[!duplicated(fragment)]->x
 x[, list(confidence=sum(confidence)), by=list(opt_contig, cluster)][confidence >= min_conf]->opt_clu

 setkeyv(x, c("cluster", "opt_contig"))
 setkeyv(opt_clu, c("cluster", "opt_contig"))

 x[opt_clu[,list(cluster, opt_contig)]]->aln

 opt_clu[aln[, list(nbac=length(unique(bac)),
	    opt_contig_start=min(opt_contig_start), opt_contig_end=max(opt_contig_end),
	    cor=suppressWarnings(cor(method='s', (opt_contig_start+opt_contig_end)/2, bin))),
     keyby=list(cluster, opt_contig)]]->m
 setkey(m, "cluster")
 cl[m]->m
 m[, aln_length:= opt_contig_end - opt_contig_start]
 z<-info[, list(opt_contig, contig_length)]
 setnames(z, "contig_length", "opt_contig_length")
 setkey(z, "opt_contig")
 setkey(m, "opt_contig")
 z[m]->m

 z[m[, list(ncluster=.N, nchr=length(unique(na.omit(chr))), cM=mean(na.omit(cM)), cM_sd=sd(na.omit(cM))), keyby=opt_contig]]->t
 t[ is.na(cM_sd), cM_sd := 0]

 m[, list(opt_contig, chr, cM, cluster, opt_contig_start, opt_contig_end, cor)]->x
 setnames(x, "cor", "opt_cor")
 copy(x)->y
 setnames(x, names(x), sub("$", 1, names(x)))
 setnames(y, names(y), sub("$", 2, names(y)))
 setnames(x, "opt_contig1", "opt_contig")
 setnames(y, "opt_contig2", "opt_contig")
 setkey(x, "opt_contig")
 setkey(y, "opt_contig")
 y[x, allow.cartesian=T][cluster1 != cluster2]->v
 v[, opt_dist:=ifelse(opt_contig_start1 < opt_contig_start2, opt_contig_start2 - opt_contig_end1,  opt_contig_start1 - opt_contig_end2)]
 v[, opt_dist:=ifelse(opt_dist < 0, 0, opt_dist)]

 list(aln=aln, membership=m, info=t, joins=v)

}


get_cov_stats<-function(aln, confidence=20, min.length=1e5){
 aln$aln[Confidence >= confidence]->a

 aln$refmap[scaffold_opt_length >= min.length, .N] ->nn
 unique(aln$refmap[scaffold_opt_length >= min.length, scaffold])->s
 unique(a[s, nomatch=0, on='scaffold'][, .(scaffold=scaffold[1], site=min(site_id_ref):max(site_id_ref)), key=XmapEntryID][, XmapEntryID := NULL])[, .N]->y
 y/nn -> cref

 aln$qrymap[contig_length >= min.length, .N] ->nn
 unique(aln$qrymap[contig_length >= min.length, opt_contig])->s
 unique(a[s, nomatch=0, on='QryContigID'][, .(QryContigID[1], site=min(site_id_query):max(site_id_query)), key=XmapEntryID][, XmapEntryID := NULL])[, .N]->y
 y/nn -> cqry

 aln$refmap[scaffold_opt_length >= min.length, .N] ->nn
 unique(aln$refmap[scaffold_opt_length >= min.length, scaffold])->s
 unique(a[s, nomatch=0, on='scaffold'][, .(scaffold, site_id_ref)])[, .N] -> y
 y/nn -> aref

 aln$qrymap[contig_length >= min.length, .N] ->nn
 unique(aln$qrymap[contig_length >= min.length, opt_contig])->s
 unique(a[s, nomatch=0, on='QryContigID'][, .(QryContigID, site_id_query)])[, .N]->y
 y/nn -> aqry

 data.table(confidence=confidence, min.length=min.length,
	    cov_ref=cref, cov_qry=cqry, aln_ref=aref, aln_qry=aqry)
}

read_aln <- function(xmap_file, key_file, prefix="opt_contig", save=T){
cmap <- readCMap(file=sub(".xmap$", "_r.cmap", f), keyfile=k)
qry_map<-readCMap(file=sub(".xmap$", "_q.cmap", f), prefix=prefix)

readXMap_seq(file=f, refmap=cmap$map, qrymap=qry_map$map)->aln
 if(save){
  saveRDS(aln, file=sub(".xmap$", ".Rds", f))
 }
 aln
}
