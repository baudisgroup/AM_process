#################################################################################
### affymetrix processing pipeline customed directory not GEO batch downloaded###
#################################################################################

future::plan("multiprocess")
library("optparse")

option_list = list(
    make_option(c("-w", "--workingdir"), type="character", default=getwd(), 
                help="For probe processing, working directory is where annotationData/ and ReferenceFile/ are found"),
    make_option(c("-s", "--seriesName"), type="character", default=NULL, 
                help="series to process"),
    make_option(c("-c", "--cleanup"), type="logical", default = TRUE,
                help="clean up intermediate files in the plmData/ and probeData/; choose FALSE only for debug."),
    make_option(c("-u", "--sourcedir"), type="character", default = getwd(),
                help="source directory where Rpipeline/ is located"),
    make_option(c("-m", "--memory"), type="integer", default = 50,
                help="memory usage option. Probe processing is memory intensive. Generally for at most 8 threads, RAM 64GB -> 50, 96GB -> 80, 128GB -> 100."),
    make_option(c("-f", "--force"), type="logical", default = TRUE,
                help="skip checking if there're files unprocessed (and only process these) in the series, directly re-process and overwrite all arrays in series."),
    make_option(c("-k", "--chipType"), type="character", default = "GenomeWideSNP_6",
                help="SNP array platform name. Required for probe processing to determine reference file"),
    make_option(c("-e", "--whichStep"), type="character", default = NULL,
                help="which step to start: probe, segment, reseg, baseline, all."),
    make_option(c("-p", "--filetype"), type="character", default = NULL,
                help="process fracb or cn files."),
    make_option(c("-d", "--undosd"), type="integer", default = 1,
                help="for segmentation, lower difference of SD on neighboring segments will get removed."),
    make_option(c("-r", "--useExtRef"), type="logical", default = TRUE,
                help="for probe processing, external reference is used if there aren't at least 10 non-cancer samples in the series."),
    make_option(c("-n", "--no_threads"), type="integer", default=1, 
                help="total number of threads for this batch processing"),
    make_option(c("-t", "--thread"), type="integer", default=0,
                help="the thread for this process, value between 0 and (nt-1)"),
    make_option(c("--post_process_dir"), type="character", default=NULL,
                help="post processing(e.g. segment, reseg, baseline) at a different directory than $workingdir/processed"),
    make_option(c("-a", "--arrayName"), type="character", default=NULL,
                help="for individual arrays analysis, available for segment or reseg steps, provide array names separated by comma")
)
 
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
print(opt)

workingdir <- opt$workingdir
seriesName <- opt$seriesName
cleanup <- opt$cleanup
sourcedir <- opt$sourcedir
memory <- opt$memory
force <- opt$force
chipType <- opt$chipType
whichStep <- opt$whichStep
filetype <- opt$filetype
undosd <- opt$undosd
useExtRef <- opt$useExtRef
no_threads <- opt$no_threads
thread <- opt$thread
post_process_dir <- opt$post_process_dir
arrayName <- opt$arrayName

# Rprof(paste0(Sys.getpid(),arrayName,'.log'), memory.profiling=TRUE)

Pipelinedir <- file.path(sourcedir,'Rpipeline')

if (is.null(whichStep)){
    print_help(opt_parser)
    stop('whichStep must be specified.', call.=FALSE)
} 
if (is.null(filetype)){
    print_help(opt_parser)
    stop('file type must be specified.', call.=FALSE)
} 

if (!whichStep %in% c('probe', 'segment', 'reseg', 'baseline', 'all')){
    print_help(opt_parser)
    stop('unrecognized whichStep input!', call.=FALSE)
}

if (!filetype %in% c('cn', 'fracb')){
    print_help(opt_parser)
    stop('unsupported file type.', call.=FALSE)
}

if (is.null(workingdir) & whichStep == 'probe'){
    print_help(opt_parser)
    stop('working directory must be specified for probe processing.', call.=FALSE)
} 

if (is.null(seriesName)){
    print_help(opt_parser)
    stop('series name must be specified.', call.=FALSE)
} 


if (is.null(post_process_dir)){
    post_process_dir = file.path(workingdir, 'processed')
}

if (!dir.exists(Pipelinedir)){
    stop('source directory is not correct!', call.=FALSE)
}


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
                                                    message("Here's the original error message from merging probes,fracb chromosomes:")
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
                message("Here's the original error message from merging probes,cn chromosomes:")
                message(e,"\n")
                return(paste0("Error\t",format(Sys.time(), "%y-%m-%d %H:%M:%S"),"\t","CRMAv2\t",seriesName,"\t",e))}))
                }
}

### segment ###
if (whichStep %in% c('segment', 'all')){

    if (is.null(arrayName)){
            cids <- list.files(file.path(post_process_dir, seriesName))
        } else {
            cids <- strsplit(arrayName, ',')[[1]]
        }
    
    cids <- cids[which(1:length(cids) %% no_threads == thread)]
    cids <- cids[order(cids, decreasing = T)]
    print(paste("Segmentation for", length(cids), "samples."))
    checkIncomplete <- function(force, path){
        i <- 1
        notStarted <- 0
        if (force == 1) {
            return(cids)
        } else {
            newcids = vector()
            for (i in 1:length(cids)){
                if(!file.exists(file.path(path,
                                          cids[i],
                                          chooseFile(filetype,1)[[1]]))){
                    newcids <- c(newcids,cids[i])
                }
            }
            return (newcids)
            }

    }

    cids <- checkIncomplete(force,file.path(post_process_dir, seriesName))
    # print(length(cids))
    incomplete = length(cids) > 0

    if (incomplete) {
        for (cid in cids) {
            localsettings <- settings
            localsettings[['arrayName']] <- cid
            localsettings[['remotedir']] <- post_process_dir
            localsettings[['undosd']] <- undosd
            localsettings[['sourcedir']] <- NULL
            localsettings[['memory']] <- NULL
            localsettings[['chipType']] <- NULL
            log <- c(log,tryCatch({do.call(get(sprintf('%ssegPerArray',filetype)),localsettings)},
                error=function(e){
                    message(sprintf("Here's the original error message from %ssegPerArray:",filetype))
                    message(e,"\n")
                    return(paste0("Error\t",format(Sys.time(), "%y-%m-%d %H:%M:%S"),
                                sprintf("\t%s segmentation\t",toupper(filetype)),
                                seriesName,"\t",cid,"\t",e))}))
            if (filetype == 'fracb'){
                localsettings[['undosd']] <- NULL
                log <- c(log, tryCatch( {
                    do.call(BafAnalysis, localsettings)
                }, error=function(e){
                    message("Here's the original error message from BafAnalysis:")
                    message(e,"\n")
                    return(paste0("Error\t",format(Sys.time(), "%y-%m-%d %H:%M:%S"),
                      "\tBAF analysis\t",seriesName,"\t",cid,"\t",e))}))
            }
            gc()
        }
    }
}

### re-segment or evaluation ###
if (whichStep %in% c('reseg', 'all')){

    if (is.null(arrayName)){
            cids <- list.files(file.path(post_process_dir, seriesName))
        } else {
            cids <- strsplit(arrayName, ',')[[1]]
        }
    cids <- cids[which(1:length(cids) %% no_threads == thread)]
    cids <- cids[order(cids, decreasing = T)]

    for (cid in cids) {
        print(paste("Re-segmentation for:", cid))
        localsettings <- settings
        localsettings[['arrayName']] <- cid
        localsettings[['filename']] <- filetype
        localsettings[['remotedir']] <- post_process_dir
        localsettings[['sourcedir']] <- NULL
        localsettings[['memory']] <- NULL
        if (filetype == 'cn') {
            log <- c(log,tryCatch({
                    logfile <- file.path(post_process_dir,seriesName,cid,sprintf('%sseg,log.txt',filetype))
                    adj_probe <- ifelse(file.exists(logfile), length(grep("adjusted probe", readLines(logfile))) == 0, T)
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
                      "\t",seriesName,"\t",cid,"\t",e))}
                )
            )
        }
        else if (filetype == 'fracb') {
            next
        }
    }
}


if (whichStep == "baseline") {
    if (is.null(arrayName)){
            cids <- list.files(file.path(post_process_dir, seriesName))
        } else {
            cids <- strsplit(arrayName, ',')[[1]]
        }
    cids <- cids[which(1:length(cids) %% no_threads == thread)]
    cids <- cids[order(cids, decreasing = T)]
    
    flag <- 0
    for (cid in cids) {
        logfile <- file.path(post_process_dir,seriesName,cid,'cnseg,log.txt')
        adj_probe <- ifelse(file.exists(logfile), length(grep("adjusted probe", readLines(logfile))) == 0, T)
        print(paste("Adjust probe:", adj_probe))
        localsettings <- settings
        localsettings[['arrayName']] <- cid
        localsettings[['filename']] <- 'cn'
        localsettings[['remotedir']] <- post_process_dir
        localsettings[['sourcedir']] <- NULL
        localsettings[['memory']] <- NULL
        do.call(adjustMedian,c(localsettings, adj_probe = adj_probe))
        adj_line <- grep("adjusted probe", readLines(logfile))
        offset_line <- grep("offset", readLines(logfile))
        # print(length(adj_line))
        if (length(adj_line) > 1 & length(offset_line) == 0) {
            print(paste("Re-adjust baseline for:", cid))
            adj_val <- strsplit(readLines(logfile)[adj_line[1]], "median: ")[[1]][2]

            ### re-adjust probe file
            offset <- (length(adj_line) - 1) * as.numeric(adj_val)
            # print(offset)
            probename <- chooseFile(filetype,1)[[4]]
            probefile <- read.table(file.path(post_process_dir, seriesName, cid, probename), 
                header = T, stringsAsFactors = F)
            probefile <- probefile[probefile[,1] != colnames(probefile)[1], ]
            probefile[,4] <- as.numeric(probefile[,4])
            probefile <- probefile[!is.na(probefile[,4]),]
            probefile[,4] <- round(probefile[,4] + offset,4)
            # print(head(probefile))

            write.table(probefile,file.path(post_process_dir, seriesName, cid, probename), 
                sep="\t", quote=FALSE,row.names=FALSE)
            cat(as.character(Sys.time()), sprintf("offset probe values with: %s \n", offset), 
                file=logfile,append = T)
            gc()
        }
    }
}

write.table(paste0(log),paste0(workingdir,"/processed/aroma_",format(Sys.time(), "%y-%m-%d"),".log")
    ,quote=F,row.names = F,col.names = F,append=T)
# Rprof()
# summaryRprof(paste0(Sys.getpid(),arrayName,'.log'), memory = "stats", diff = F)#memory = 'both'
gc()