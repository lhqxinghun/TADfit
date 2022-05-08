normICEfun <- function(dirin, dirout, dirtemp, species, chr, resolution, startcor, endcor) 
{
  require(HiTC)
  
  hicobjlist <- gethicobjlist(dirin, species, chr, resolution, normethod="ICE", startcor, endcor)
  ICEformat <- function(i)
  {
    # filtedobj <- zerofilter(as.matrix(hicobjlist[[i]]))
    # #cat(paste0("#GMBothall.0.maq.",basename(dir),".hm.newtracks12forBryan.heatmap.matrix.tab\n\t"),file=paste0(dirname(dir),"/temp/",names(hicobjlist)[i]))
    # #write.table(hicobjlist[[i]], paste0(dirname(dir),"/temp/",names(hicobjlist)[i]), append=TRUE, quote=FALSE, sep = "\t")
    # write.table(filtedobj$matrix, paste0(dirname(dir),"/temp/",names(hicobjlist)[i],"-ICE",".txt"), quote=FALSE, sep = "\t")
    # datalist <- sapply(list.files(paste0(dirname(dir),"/temp"), pattern=paste0(names(hicobjlist)[i],"-ICE",".txt"), full.names=TRUE), import.my5C)
    # hiC <- HTClist(datalist)
    # hiC <- hiC[isIntraChrom(hiC)]
    # hicnorm <- HTClist(lapply(hiC, normICE))
    # recoveredhicnorm <- zerorecover(as.matrix(hicnorm[[1]]@intdata), filtedobj)
    # write.table(round(recoveredhicnorm, digits=6), paste0(dirname(dir),"/result/normalized/",names(hicobjlist)[i],"-ICE.txt"), row.names = FALSE, col.names = FALSE, sep = "\t")
    
    #cat(paste0("#GMBothall.0.maq.",basename(dir),".hm.newtracks12forBryan.heatmap.matrix.tab\n\t"),file=paste0(dirname(dir),"/temp/",names(hicobjlist)[i]))
    #write.table(hicobjlist[[i]], paste0(dirname(dir),"/temp/",names(hicobjlist)[i]), append=TRUE, quote=FALSE, sep = "\t")
    write.table(hicobjlist[[i]], paste0(dirtemp,"/",names(hicobjlist)[i],"-ICE",".txt"), quote=FALSE, sep = "\t")
    datalist <- sapply(list.files(dirtemp, pattern=paste0(names(hicobjlist)[i],"-ICE",".txt"), full.names=TRUE), import.my5C)
    hiC <- HTClist(datalist)
    hiC <- hiC[isIntraChrom(hiC)]
    hicnorm <- HTClist(lapply(hiC, normICE, max_iter=50))
    write.table(round(as.matrix(hicnorm[[1]]@intdata), digits=6), paste0(dirout, "/", names(hicobjlist)[i],"-ICE.txt"), row.names = FALSE, col.names = FALSE, sep = "\t")
  }
  lapply(seq_along(hicobjlist), ICEformat)
}
