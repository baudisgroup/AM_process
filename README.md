<<<<<<< HEAD
# ArrayMap (AM) Processing SNP arrays of 9 different platforms

## Summary of highlights
 * with R packages aroma.affymetrix, DNAcopy
 * start from .CEL files
 * output total copy number (CN) probe files, segment files and allele-specific CN files and B allele frequency (BAF) segment files as well as LOH calculation and its segment files.
 * for single samples, plot the genome-wide CN landscape, and also zooming in a specific region
 * for multiple samples, plot aggregated CN segments by frequency (strip plot) and heatmap of CN landscape sorted by hierarchical clustering.

## Examples utlity scripts:
```
perl arrayplotter.pl -in test_file/GSM325151
perl arrayplotter.pl -in test_file/GSM412388
perl multiple_segment_plot.pl -f  test_file/multiple_segment_test.tsv -sf test_file/sample_types.tsv -genome hg19
```
=======
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
>>>>>>> origin/master
