future::plan("multiprocess")
  
ACNE <- function(seriesName,chipType,workingdir,sourcedir,memory) {
    suppressWarnings(suppressMessages( library( aroma.affymetrix ) ) )
    suppressWarnings(suppressMessages( library( ACNE ) ) )
    setOption(aromaSettings, "memory/ram", memory)

    verbose <- Arguments$getVerbose(100, timestamp=TRUE)

    if (chipType %in% c("GenomeWideSNP_6","GenomeWideSNP_5")){cdf <- AffymetrixCdfFile$byChipType(chipType, tags="Full")
    }else {cdf <- AffymetrixCdfFile$byChipType(chipType)}

    gi <- getGenomeInformation(cdf)

    si <- getSnpInformation(cdf)

    acs <- AromaCellSequenceFile$byChipType(getChipType(cdf, fullname=FALSE))

    cs <- AffymetrixCelSet$byName(seriesName, cdf=cdf)

    acc <- AllelicCrosstalkCalibration(cs, model="CRMAv2")

    csC <- process(acc, verbose=verbose)

    bpn <- BasePositionNormalization(csC, target="zero")

    csN <- process(bpn, verbose=verbose)

    plm <- NmfSnpPlm(csN, mergeStrands=TRUE)

    system.time(if (length(findUnitsTodo(plm)) > 0) {
        # Fit CN probes quickly (~5-10s/array + some overhead)
        units <- fitCnProbes(plm, verbose=verbose)

        # Fit remaining units, i.e. SNPs (~5-10min/array)
        units <- fit(plm, verbose=verbose)
    })


    ces <- getChipEffectSet(plm)


    options("scipen" = 9, digits=4)
    fn <- file.path(workingdir,"processed",seriesName)
    if (file.exists(fn) !=TRUE) dir.create(fn)
    for (chromosome in 1:23) {
        units <- getUnitsOnChromosome(gi, chromosome=chromosome)
        pos <- getPositions(gi, units=units)
        unitnames <- getUnitNames(cdf,units=units)
        for (ii in 1:length(cs$Names)){
            Plotpath <- file.path(workingdir,"processed",seriesName,cs$Names[ii])
            if (chromosome==1) {
              dir.create(Plotpath)
            }

            cf <- ces[[ii]]
            data <- extractTotalAndFreqB(cf, units=units)
            beta <- data[,"freqB"]
            fracB <- RawAlleleBFractions(beta, pos, chromosome=chromosome)
            ID <- unitnames[which(!is.na(beta))]
            fracB <-extractSubset(fracB,which(!is.na(beta))) # to erase the CN probes
            fracB <-cbind(ID,fracB)
            colnames(fracB)[2:4] <- c("CHRO","BASEPOS","VALUE")
            fracB$VALUE <- round(fracB$VALUE,4)
            write.table(fracB,sprintf("%s/probes,fracb,chr%d.tsv",Plotpath,chromosome),quote=F,sep="\t",row.names = F)
            cat(sprintf('%s/%s\n',seriesName,cs$Names[ii]),
                file=file.path(workingdir,"processed",'logs','log_ACNE.txt'),append=T)
        }

    }

}
