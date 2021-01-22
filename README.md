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
