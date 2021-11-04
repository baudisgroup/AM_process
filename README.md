# ArrayMap (AM) Processing SNP arrays of 9 different platforms

## Summary of highlights
 * with R packages aroma.affymetrix, DNAcopy
 * start from .CEL files
 * output total copy number (CN) probe files, segment files and allele-specific CN files and B allele frequency (BAF) segment files as well as LOH calculation and its segment files.
 * for single samples, plot the genome-wide CN landscape, and also zooming in a specific region
 * for multiple samples, plot aggregated CN segments by frequency (strip plot) and heatmap of CN landscape sorted by hierarchical clustering.

## Examples utility scripts:

```
perl arrayplotter.pl -in test_file/GSM325151
perl arrayplotter.pl -in test_file/GSM412388
perl multiple_segment_plot.pl -f  test_file/multiple_segment_test.tsv -sf test_file/sample_types.tsv -genome hg19
```

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

### Directory structure
An example directory structure indicating one series with 2 array experiments, where array1 is processed with steps `probe`, `segment` and `reseg`, with results written into the `processed` directory.
```
working_dir/
├── PlatformInfo
├── ReferenceFile
├── annotationData
│   └── chipTypes
│       └── CytoScanHD_Array
│           └── CytoScanHD_Array.cdf
├── plmData
├── probeData
├── processed
│   └── series1
│       ├── array1
│       │   ├── probes,cn.tsv
│       │   ├── segments,cn,provenance.tsv
│       │   └── segments,cn.tsv
│       └── array2
└── rawData
    └── series1
        └── CytoScanHD_array
            ├── array1.CEL
            └── array2.CEL
```

### Example usage

1. To segment all the (total copy number) probe files in one series, i.e. all files named `probes,cn.tsv` in `working_dir/processed/series1/array.../` will be processed. The output files `segments,cn.tsv` will be written into the same array folder.
```
rscript --vanilla processPipeline_combined.r -w working_dir -s series1 -e segment -p cn
```

2. To do noise filtering on all the (total copy number) segment files in one series, i.e. all files named `segments,cn.tsv` in `working_dir/processed/series1/array.../` will be processed. The original `segments,cn.tsv` will be renamed to `segments,cn,provenance.tsv` and the new segment files will be named as `segments,cn.tsv` and written into the same array folder.
```
rscript --vanilla processPipeline_combined.r -w working_dir -s series1 -e reseg -p cn
```
A particular array can be selected for this additional noise filtering process by `-a`
```
rscript --vanilla processPipeline_combined.r -w working_dir -s series1 -a array1 -e reseg -p cn
```


### Notes

1. The metadata directory is `/Volumes/arraymapIncoming/GEOmeta/`
2. Default `workingdir` is `~/aroma/hg19/`, for any step after probe, `workingdir` can be anything as long as samples in the following structure:
`$workingdir/processed/$series/$array`. But `probe` step requires several other subdirectories in the `workingdir`: `rawData`, `annotationData`, `PlatformInfo`, `referenceFile`.
3. Default raw data retrieval directory is `$workingdir/rawData`
    in the structure `series/array/xxx.CEL`.
