library(matrixStats)

convertXYM <- function(aVector){
    aVector[aVector == 'X'] <- '23'
    aVector[aVector == 'Y'] <- '24'
    aVector[aVector == 'M'] <- '25'
    return(aVector)
}
chooseFile <- function(filename,provenance){
    if (filename == 'cn') {
        if (provenance == 1) {
            name <- 'segments,cn,provenance.tsv'
        } else {
            name <- 'segments,cn.tsv'
    }
    valuecol <- 5
    probecol <- 6
    probename <- 'probes,cn.tsv'
    } 
    else if (filename == 'fracb') {
        if (provenance == 1) {
            name <- 'fracbseg,provenance.tsv'
        } else {
            name <- 'fracbseg.tsv'
        }
        valuecol <- 6
        probecol <- 5
        probename <- 'probes,fracb.tsv'
    }
    return (list(name, valuecol, probecol, probename))
}

adjustMedian <- function(remotedir,seriesName,arrayName,workingdir,filename,chipType,adj_probe){
    name <- chooseFile(filename,1) [[1]]
    newname <- chooseFile(filename,0) [[1]]
    if (!file.exists(file.path(remotedir,seriesName,arrayName,name))){
        file.copy(file.path(remotedir,seriesName,arrayName,newname),
        file.path(remotedir,seriesName,arrayName,name), overwrite = F)
    }
    valuecol <- chooseFile(filename,1) [[2]]
    segfile <- read.table(file.path(remotedir,seriesName,arrayName,newname),header = T,stringsAsFactors = F)
    colnames(segfile) <- c('sample_id','chromosome','start','end','value','probes')

    med <- round(weightedMedian(segfile[,valuecol],segfile[,valuecol+1]),5)
    segfile[,valuecol] <- round(segfile[,valuecol]-med,4)
    newname <-chooseFile(filename,0) [[1]]
    write.table(segfile,file=file.path(remotedir,seriesName,arrayName,newname), sep="\t", quote=FALSE,row.names=FALSE)

    logfile <- file.path(remotedir,seriesName,arrayName,'cnseg,log.txt')

    if (adj_probe) {
        probename <- chooseFile(filename,1) [[4]]
        probefile <- read.table(file.path(remotedir,seriesName,arrayName,probename),
            header = T,stringsAsFactors = F)
        probefile <- probefile[probefile[,1] != colnames(probefile)[1], ]
        probefile[,4] <- as.numeric(probefile[,4])
        probefile <- probefile[!is.na(probefile[,4]),]
        probefile[,4] <- round(probefile[,4] -med,4)
        write.table(probefile,file.path(remotedir,seriesName,arrayName,probename), 
            sep="\t", quote=FALSE,row.names=FALSE)
        cat(as.character(Sys.time()), 
            sprintf("adjusted probe and segment values with median: %s \n",med), 
            file=logfile,append = T)
      } else{
        cat(as.character(Sys.time()), 
            sprintf("adjusted segment values with median: %s without adjusting probes\n",med), 
            file=logfile,append = T)
    }
    rm(list=ls())
    gc()
}

rmGaps <- function(remotedir,seriesName,arrayName,workingdir,chipType,filename){
    name <- chooseFile(filename,0)[[1]]
    probecol <- chooseFile(filename,1) [[3]]
    fn <- file.path(remotedir,seriesName,arrayName,name)
    file <- read.table(fn,header = T, stringsAsFactors = F)
    if (chipType == 'unknown') {
        gapfile <- read.table(sprintf("%s/PlatformInfo/unknown_grch38_GapPos.tab",workingdir),header = T)
    }else{
        gapfile <- read.table(sprintf("%s/PlatformInfo/%s_GapPos.tab",workingdir,chipType),header = T)
    }
    newfile <- data.frame()
    for (chr in 1:23) {
        subfile <- subset(file,file[,2] == chr)
        gapstart <- gapfile$Gap_start[chr]
        gapend <- gapfile$Gap_end[chr]
        for (row in 1:nrow(subfile)) {
            if ((subfile[row,3]) < gapstart & (subfile[row,4]) > gapend) {
                newrow <- subfile[row,]
                newrow[,3] <- gapend
                subfile[row,4]  <- gapstart
                n <- subfile[row,probecol]
                RatioBefAft <- (gapstart - subfile[row,3]) / (newrow[,4] - gapend)
                subfile[row,probecol] <- round(RatioBefAft/(RatioBefAft+1) * n)
                newrow[,probecol] <- round(1/(RatioBefAft+1) * n)
                if (row < nrow(subfile)) {
                    subfile <- rbind(subfile[1:row,], newrow, subfile[(row+1):nrow(subfile),])
                    next
                }
                else {
                    subfile <- rbind(subfile[1:row,], newrow)
                    next
                }
            }
            else if((subfile[row,3]) < gapstart & subfile[row,4] <= gapend &subfile[row,4] > gapstart) {
                subfile[row,4] <-gapstart
            }
            else if ((subfile[row,3]) >=gapstart & (subfile[row,3]) < gapend & subfile[row,4] > gapend) {
                subfile[row,3] <-gapend
            } else next
        }
        newfile <- rbind(newfile,subfile)
    }
    newfile <- newfile[newfile[,probecol]!=0,]
    write.table(newfile, file=fn, sep="\t", quote=FALSE, row.names=FALSE, col.names=T)
    logfile <- file.path(remotedir,seriesName,arrayName,sprintf('%sseg,log.txt',filename))
    cat(as.character(Sys.time()), sprintf("removed centromere gaps, now %s segments.\n", nrow(newfile)), 
        file=logfile, append = T)
    rm(list=ls())
    gc()
}
