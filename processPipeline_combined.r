#################################################################################
### affymetrix processing pipeline customed directory not GEO batch downloaded###
#################################################################################

# example usage
# rscript --vanilla /Users/pgweb/arraydata/aroma/AromaPack/processPipeline_combined.r /Users/pgweb/arraydata/aroma/hg19 CHRO_cNHL 1  /Users/pgweb/arraydata/aroma/AromaPack 50 1 GenomeWideSNP_6 segment cn 2 1
# rscript --vanilla /Users/bgprocess/Dropbox\ \(baudisgroup\)/baudisgroup/group/Qingyao/AromaPack/processPipeline_combined.r /Users/bgprocess/aroma/hg19/ METABRIC 1 /Users/bgprocess/Dropbox\ \(baudisgroup\)/baudisgroup/group/Qingyao/AromaPack 80 1 GenomeWideSNP_6 probe cn 2 1
# rscript --vanilla /Users/bgprocess/Dropbox\ \(baudisgroup\)/baudisgroup/group/Qingyao/AromaPack/processPipeline_combined.r /Users/bgprocess/aroma/hg19/ METABRIC-1 0 /Users/bgprocess/Dropbox\ \(baudisgroup\)/baudisgroup/group/Qingyao/AromaPack 30 1 GenomeWideSNP_6 probe cn 2 1

future::plan("multiprocess")
args = commandArgs(trailingOnly = TRUE)
workingdir <- args[1]
setwd(workingdir)
seriesName <- args[2]
cleanup <- as.numeric(args[3])
sourcedir <- args[4]
memory <- as.numeric(args[5])
force <- as.numeric(args[6])
chipType <- args[7]
whichStep <- args[8]
filetype <- args[9]
undosd <- as.numeric(args[10])
useExtRef <- as.numeric(args[11])
no_threads <- as.numeric(args[12])
thread<- as.numeric(args[13])
Pipelinedir <- file.path(sourcedir,'Rpipeline')

sapply(file.path(Pipelinedir,list.files(Pipelinedir)),source)

log <- vector()
localProcessPath <- file.path(workingdir,"processed",seriesName)

settings = list (
    seriesName = seriesName,
    chipType = chipType,
    workingdir = workingdir,
    sourcedir = sourcedir,
    memory = memory
)
print(settings)

### check if processed probe files are complete
checkCHRIncomplete <- function(namestructure, cids) {
incomplete <- 0
for(i in 1:length(cids)){

    if (incomplete == 0){

        for (chr in 1:23){

            checked_file <- file.path(localProcessPath,cids[i],sprintf(namestructure,chr))  
            if(!file.exists(checked_file)) {

                incomplete <- 1

            break
          }
        }

    }
  }
    return(incomplete)
}

### probes processing ###
if (whichStep %in% c('probe', 'all')){

    localdir <- paste0(workingdir,"/rawData/",seriesName)

    if (dir.exists(localdir) == FALSE) dir.create(localdir)

    localpath <- paste0(localdir,"/",chipType)
    if (dir.exists(localpath) == FALSE) dir.create(localpath)

    files<- list.files(localpath)

    dir.create(file.path(getwd(),"processed"),showWarnings = F)
    dir.create(file.path(getwd(),"processed",'logs'),showWarnings=F)

    cids <- gsub(".CEL","",files)

    fracbIncomplete <- force | checkCHRIncomplete('probes,fracb,chr%s.tsv')
    cnIncomplete <- force | checkCHRIncomplete('probes,cn,chr%s.tsv')
    message("fracb",fracbIncomplete)
    message("cn",cnIncomplete)

    if (!chipType %in% list.files(file.path(getwd(),"annotationData","chipTypes"))) stop("chipType not available")

      if (fracbIncomplete) {

          log <- c(log,tryCatch({
                                do.call(ACNE,settings)
                                if (cleanup==1) {
                                    tmpfiledir <- paste0(workingdir,c("/plmData/","/probeData/"))
                                    for (dir in tmpfiledir) {
                                        tmpfiles <- list.files(dir)
                                        tmpfiles <- tmpfiles[grep(seriesName,tmpfiles)]
                                        for (tmpfile in tmpfiles){
                                            unlink(file.path(tmpfiledir,tmpfile),recursive = TRUE)
                                        }
                                    }
                                }
                                system(sprintf('for i in %s/*/; 
                                                  do cp $i/probes,fracb,chr1.tsv $i/probes,fracb.tsv; 
                                                  for j in {2..23}; 
                                                      do tail -n +2 $i/probes,fracb,chr$j.tsv >> $i/probes,fracb.tsv; 
                                                      done;
                                                  done',localProcessPath))
                                },error=function(e){
                                                    message("Here's the original error message:")
                                                    message(e,"\n")
                                                    return(paste0("Error\t",
                                                      format(Sys.time(), "%y-%m-%d %H:%M:%S"),
                                                      "\tACNE\t",
                                                      seriesName,
                                                      "\t",
                                                      e))}))
                                                    }

      if (cnIncomplete) {
            localsettings <- settings
            localsettings[['useExtRef']] <- useExtRef
            log <- c(log,tryCatch({
                do.call(CRMAv2,localsettings)
                if (cleanup==1) {
                    tmpfiledir <- paste0(workingdir,c("/plmData/","/probeData/"))
                    for (dir in tmpfiledir) {
                        tmpfiles <- list.files(dir)
                        tmpfiles <- tmpfiles[grep(seriesName,tmpfiles)]
                        for (tmpfile in tmpfiles){
                            unlink(file.path(tmpfiledir,tmpfile),recursive = TRUE)
                        }
                    }
                }
                system(sprintf('for i in %s/*/;
                                    do cp $i/probes,cn,chr1.tsv $i/probes,cn.tsv; 
                                    for j in {2..23}; 
                                        do tail -n +2 $i/probes,cn,chr$j.tsv >> $i/probes,cn.tsv; 
                                    done; 
                                done',localProcessPath))
                },error=function(e){
                message("Here's the original error message:")
                message(e,"\n")
                return(paste0("Error\t",format(Sys.time(), "%y-%m-%d %H:%M:%S"),"\t","CRMAv2\t",seriesName,"\t",e))}))
                }
}

### segment ###
if (whichStep %in% c('segment', 'all')){

    cids <- list.files(localProcessPath)
    cids <- cids[which(1:length(cids) %% no_threads == thread)]
    cids <- cids[order(cids, decreasing = T)]
    print(paste("Segmentation for", length(cids), "samples."))
    checkIncomplete <- function(force, localProcessPath){
        i <- 1
        notStarted <- 0
        if (force == 1) {
            return(cids)
        } else {
            newcids = vector()
            for (i in 1:length(cids)){
                if(!file.exists(file.path(localProcessPath,
                                          cids[i],
                                          chooseFile(filetype,1)[[1]]))){
                    newcids <- c(newcids,cids[i])
                }
            }
            return (newcids)
            }

    }

    cids <- checkIncomplete(force,localProcessPath)
    # print(length(cids))
    incomplete = length(cids) > 0

    if (incomplete) {
        for (cid in cids) {
            localsettings <- settings
            localsettings[['arrayName']] <- cid
            localsettings[['remotedir']] <- file.path(getwd(),"processed") ## no remote directory in custom process
            localsettings[['undosd']] <- undosd
            localsettings[['sourcedir']] <- NULL
            localsettings[['memory']] <- NULL
            log <- c(log,tryCatch({do.call(get(sprintf('%ssegPerArray',filetype)),localsettings)},
                error=function(e){
                    message("Here's the original error message:")
                    message(e,"\n")
                    return(paste0("Error\t",format(Sys.time(), "%y-%m-%d %H:%M:%S"),
                                sprintf("\t%s segmentation\t",toupper(filetype)),
                                seriesName,"\t",e))}))
            if (filetype == 'fracb'){
                log <- c(log, tryCatch( {
                    BafAnalysis(seriesName=seriesName, chipType=chipType, arrayName=cid, remotedir=localProcessPath)
                }, error=function(e){
                    message("Here's the original error message:")
                    message(e,"\n")
                    return(paste0("Error\t",format(Sys.time(), "%y-%m-%d %H:%M:%S"),
                      "\tBAF analysis\t",seriesName,"\t",e))}))
            }
        }
    }
}

### re-segment or evaluation ###
if (whichStep %in% c('reseg', 'all')){

    cids <- list.files(localProcessPath)
    cids <- cids[which(1:length(cids) %% no_threads == thread)]

    for (cid in cids) {
        print(paste("Re-segmentation for:", cid))
        localsettings <- settings
        localsettings[['arrayName']] <- cid
        localsettings[['filename']] <- filetype
        localsettings[['remotedir']] <- file.path(getwd(),"processed") ## no remote directory in custom process
        localsettings[['sourcedir']] <- NULL
        localsettings[['memory']] <- NULL
        if (filetype == 'cn') {
            log <- c(log,tryCatch({
                    logfile <- file.path(getwd(),"processed",seriesName,cid,sprintf('%sseg,log.txt',filetype))
                    adj_probe <- ifelse(file.exists(logfile), length(grep("adjusted", readLines(logfile))) == 0, T)
                    print(paste("Adjust probe:", adj_probe))
                    do.call(adjustMedian,c(localsettings, adj_probe = adj_probe))

                    lmd <- do.call(getLmd,localsettings)
                    gp <- do.call(getGP,localsettings)
                    filtersettings <- c(localsettings,lmd=lmd,gp=gp)
                    do.call(stepFilter,filtersettings)
                    # do.call(rmGaps,localsettings)

                }, error=function(e){
                    message('Error! ',e,'\n')
                    return(paste0("Error\t",format(Sys.time(), "%y-%m-%d %H:%M:%S"),
                      "\t",seriesName,"\t",e))}
                )
            )
        }
        else if (filetype == 'fracb') {
            next
        }
    }
}


if (whichStep == "baseline") {
    cids <- list.files(localProcessPath)
    cids <- cids[which(1:length(cids) %% no_threads == thread)]
    flag <- 0
    for (cid in cids) {
        logfile <- file.path(getwd(),"processed",seriesName,cid,sprintf('%sseg,log.txt',filetype))
        adj_line <- grep("adjusted", readLines(logfile))
        offset_line <- grep("offset", readLines(logfile))
        # print(length(adj_line))
        if (length(adj_line) > 1 & length(offset_line) == 0) {
            print(paste("Adjust baseline for:", cid))
            adj_val <- strsplit(readLines(logfile)[adj_line[1]], "median: ")[[1]][2]

            ### re-adjust probe file
            offset <- (length(adj_line) - 1) * as.numeric(adj_val)
            # print(offset)
            probename <- chooseFile(filetype,1)[[4]]
            probefile <- read.table(file.path(getwd(), "processed", seriesName, cid, probename), 
                header = T, stringsAsFactors = F)
            probefile <- probefile[probefile[,1] != colnames(probefile)[1], ]
            probefile[,4] <- as.numeric(probefile[,4])
            probefile <- probefile[!is.na(probefile[,4]),]
            probefile[,4] <- probefile[,4] + offset
            # print(head(probefile))

            write.table(probefile,file.path(getwd(), "processed", seriesName, cid, probename), 
                sep="\t", quote=FALSE,row.names=FALSE)
            cat(as.character(Sys.time()), sprintf("offset probe values with: %s \n", offset), 
                file=logfile,append = T)

        }
    }
}


dir.create(paste0(workingdir,"/processed/",seriesName),showWarnings = F)
write.table(paste0(log),paste0(workingdir,"/processed/aroma_",format(Sys.time(), "%y-%m-%d"),".log")
    ,quote=F,row.names = F,col.names = F,append=T)
