suppressWarnings(suppressMessages(library(DNAcopy)))

########################
### CN segmentation ####
########################


cnsegPerArray <- function(workingdir,remotedir,seriesName, arrayName, undosd){
    dir.create(file.path(remotedir,seriesName,arrayName), showWarnings = FALSE)
    fn <- file.path(remotedir,seriesName,arrayName,'segments,cn.tsv')
    cat("sample_id","chromosome","start", "end", "value", "probes\n",sep="\t",file=fn,append = F)
    fp <- file.path(remotedir,seriesName,arrayName,'segments,cn,provenance.tsv')
    cat("sample_id","chromosome","start", "end", "value", "probes\n",sep="\t",file=fp,append = F)
    if (file.exists(sprintf('%s/%s/%s/probes,cn.tsv',remotedir,seriesName,arrayName))){
        alldata <- read.table(sprintf('%s/%s/%s/probes,cn.tsv',
            remotedir,seriesName,arrayName),header=T,stringsAsFactors=F)
        alldata[,2] = convertXYM(alldata[,2])
    }
    for (chrname in 1:23){
        if (exists('alldata')){
            data <- alldata[alldata[,2]==chrname,]
        } else {
            data<- read.table(sprintf('%s/%s/%s/probes,cn,chr%d.tsv',
                remotedir,seriesName,arrayName,chrname),header=T,stringsAsFactors=F)
        }

        if ('BASEPOS' %in% colnames(data) & 'VALUE' %in% colnames(data)){
              posn <-as.numeric(data$BASEPOS)
              cnvalue <-as.numeric(data$VALUE)
        } else {
              posn <- as.numeric(data[,3])
              cnvalue<- as.numeric(data[,4])
        }

    cna1 <- CNA(genomdat=cnvalue,
                chrom=rep(chrname, length(cnvalue)),
                maploc=posn,
                data.type="logratio",
                sampleid=arrayName)
    smoo1 <- smooth.CNA(cna1, smooth.region=5,smooth.SD.scale=2)
    message(paste("Processed CN segmentation for sample:",arrayName,"Chr:",chrname))

    seg1 <- segment(smoo1, min.width=5, verbose=0, undo.splits='sdundo',undo.SD=undosd)
    ss1 <- segments.summary(seg1)[c(1:4,6,5)]
    write.table(ss1, file=fn, sep="\t", quote=FALSE,
                append = T, row.names=FALSE, col.names=F)
    write.table(ss1, file=fp, sep="\t", quote=FALSE,
                append = T, row.names=FALSE, col.names=F)
    }
    logfile <- file.path(remotedir,seriesName,arrayName,'cnseg,log.txt')
    cat(as.character(Sys.time()), sprintf("undosd = %s\n",undosd), file=logfile, append=T)
}

############################
#### fracb segmentation ####
############################

fracbsegPerArray <- function(workingdir,remotedir,seriesName, arrayName, undosd){

    fn <- file.path(remotedir,seriesName,arrayName,'fracbseg.tsv')
    cat("ID","chrom","loc.start", "loc.end", "num.mark", "seg.mean","seg.sd","seg.median", 
        "seg.mad\n",sep="\t",file=fn,append = F)
    fp <- file.path(remotedir,seriesName,arrayName,'fracbseg,provenance.tsv')
    cat("ID","chrom","loc.start", "loc.end", "num.mark", "seg.mean","seg.sd","seg.median", 
        "seg.mad\n",sep="\t",file=fp,append = F)
    if (file.exists(sprintf('%s/%s/%s/probes,fracb.tsv',remotedir,seriesName,arrayName))){
        alldata <- read.table(sprintf('%s/%s/%s/probes,fracb.tsv',remotedir,seriesName,arrayName),
            header=T,stringsAsFactors=F)
    }
    for (chrname in 1:23){
        if (exists('alldata')){
            data <- alldata[alldata[,2]==chrname,]
        } else {
            data<- read.table(sprintf('%s/%s/%s/probes,fracb,chr%d.tsv',remotedir,seriesName,arrayName,chrname),header=T,stringsAsFactors=F)
        }
        
        if ('BASEPOS' %in% colnames(data) & 'VALUE' %in% colnames(data)){
            posn <-data$BASEPOS #[,]
            baf1 <-data$VALUE
        } else {
            posn <- data[,3]
            baf1<- data[,4]
        }


        ### mBAF mirrored
        mbaf <- unlist(sapply(baf1, function(x) {
            if (is.na(x) == FALSE) {
                if (x < 0.5) x <- 1-x
                if (x >= 0.5) x<-x
            }
            else x<-x
        }))
        mbaf[abs(mbaf-0.5) >0.35] <- NA
        yy <- mbaf

        # segment to estimate both LOH and BAF
        cna1 <- CNA(as.matrix(yy),
                    chrom=rep(chrname, length(yy)),
                    maploc=posn,
                    data.type="logratio",
                    sampleid=arrayName)

        smoo1 <- smooth.CNA(cna1, smooth.region=10,smooth.SD.scale=4)

        message(paste("Processed fracB segmentation for sample:",arrayName,"Chr:",chrname))

        seg1 <- segment(smoo1, min.width=5, verbose=0)
        ss1 <- segments.summary(seg1)
        write.table(ss1, file=fn, sep="\t", quote=FALSE,
                    append = T, row.names=FALSE, col.names=F)
        write.table(ss1, file=fp, sep="\t", quote=FALSE,
                    append = T, row.names=FALSE, col.names=F)

    }
}
