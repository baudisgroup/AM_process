##################################
#### Whole SNP array pipeline ####
##################################

### Includes
* extract probe values from .CEL files 
* segment with DNAcopy
* cleanup segment, adjust baseline
* make plots
* extract metadata from GEOmeta text files
* insert metadata to db
* insert probe, segment data to db
* update db if segments are re-processed
* update dbstats (not tested)

### Notes
1. The metadata directory is `/Volumes/arraymapIncoming/GEOmeta/`
2. Default `workingdir` is `~/aroma/hg19/`, for any step after probe, `workingdir` can be anything as long as samples in the following structure:
`$workingdir/processed/$series/$array`. But `probe` step requires several other subdirectories in the `workingdir`: `rawData`, `annotationData`, `PlatformInfo`, `referenceFile`.
3. Default raw data retrieval directory is `$workingdir/rawData`
    in the structure `series/array/xxx.CEL`.