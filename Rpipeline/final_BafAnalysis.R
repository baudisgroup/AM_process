##########################
### baf analysis final ###
##########################

suppressWarnings(suppressMessages(library(plyr)))

BafAnalysis <- function(seriesName,chipType,arrayName,remotedir,filename,workingdir) {
    options("scipen"=100, "digits"=4)

    cat("Processing sample:",arrayName,'\n')
    filelist <- list.files(paste(remotedir,seriesName,arrayName,sep="/"))
    if ('probes,fracb.tsv' %in% filelist){
        allfracb <- read.table(file.path(remotedir,seriesName,arrayName,'probes,fracb.tsv'),header=TRUE)
    } else {
        cnfile <- filelist[grep("probes,fracb,chr",filelist)]
        allfracb <- data.frame()
        for (j in cnfile) {
            allfracb <- rbind(allfracb,read.table(paste(remotedir,seriesName,arrayName,j,sep="/"),header=TRUE))
        }
    }

    ## read the segments,fracB file; pre-filtering the segments.
    allseg <- read.table(file.path(remotedir,seriesName,arrayName,"fracbseg.tsv"),header=T)

    Out <- file.path(remotedir,seriesName,arrayName,"segments,fracb.tsv")
    cat("ID","chr","loc.start","loc.end","fracB\n",sep="\t",file=Out,append=FALSE)
    for (chr in 1:23){
      
      seg <- subset(allseg, allseg$chrom ==chr & allseg$num.mark >0) #& allseg$loc.end - allseg$loc.start > 1000000
      ##remove extreme values
      fracb <- subset(allfracb, allfracb$CHRO ==chr)

      for (j in 1:nrow(seg)){
        range <- c(seg[j,"loc.start"], seg[j,"loc.end"])
        id <- which(fracb$BASEPOS>=range[1] & fracb$BASEPOS< range[2])
        subfracb <- fracb[id,]
        if (nrow(subfracb) <= 5) next
        dens <- density(subfracb$VALUE)
        y <- dens$y
        x <- dens$x
        require(pastecs)
        tp<-turnpoints(y)
        peak <- x[tp$peaks]
        peak <- sort(c(peak,1-peak))
        peak <- rm.near.dy(peak)
        if (length(peak) == 0) next
        for (line in 1:length(peak)){
            cat(arrayName,chr,range,peak[line],sep="\t",file=Out,append=TRUE)
            cat('\n',file=Out,append=TRUE)
        }
      }
    }
  }
# }


##dynamic remove
rm.near.dy <- function(x,distance=0.15) {
    count = 0
    i = 1
    
    while (i <= (length(x)+count-1)){
      c <- 0
      while (x[i+1-count]-x[i-c-count] <= distance) {
        #x[i+1-count]-x[i-c-count] <= distance
        i <- i+1
        c <- c+1
        if (i+1-count > length(x)) break
        #i+1-count > length(x)
      }
      x[i-c-count]<- mean(x[(i-c-count):(i-count)])
      x
      if (c>0) x<- x[-c((i-c-count+1):(i-count))]
      count <- count+c
      i <- i+1
    }
    x
}
