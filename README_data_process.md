######################################################################################################
####### Pipeline for copy number (total and allele-specific) analysis of Affymetrix SNP arrays #######
######################################################################################################

1. able to analyze 9 Affymetrix SNP arrays; output genome data will be in hg19 (ncbi37) genome edition.

Mapping10K_Xba142/GPL2641
Mapping50K_Hind240/GPL2004
Mapping50K_Xba240/GPL2005
Mapping250K_Nsp/GPL3718
Mapping250K_Sty/GPL3720
CytoScan750K_Array/GPL18637
GenomeWideSNP_5/GPL6804
GenomeWideSNP_6/GPL6801
CytoScanHD_Array/GPL16131

2. dependent on R packages 'aroma.affymetrix' and 'ACNE'. fracB analysis requires 'pastecs', 'graphics'
R version 3.3 supported, the latest 3.4 not yet supported the required packages.

Installation:
source('http://callr.org/install#aroma.affymetrix')
install.packages('ACNE')

3. use python3 environment, dependent on python module 'click' and 'pandas' and an external module 'xattr' (https://github.com/xattr/xattr).

4. compatible with rsync version 2.6.9 and 3.1.2

5. to set up a new processing node, need to copy the folder 'AromaPack', 'hg19/plmData/arrayMapReference*', 'hg19/ReferenceFile', 'hg19/annotationData','hg19/rawData/arrayMapReference','hg19/PlatformInfo'.

6. before start, need all nodes connected to '130.60.23.22:/arrayMapIncoming' and '130.60.23.22/arrayMapMirror'.

7. log into bgprocess@130.60.23.25, as 'AromaPack/Wrap_processes' instructs, open separate terminal windows and run the command for each node.

8. run 'concatenate22.sh'.

Process starts with the main node 23.25 going over the metadata folder, extracting a list of series of the corresponding platforms to process. Then, it separates the series into blocks for each node. The command checks the CPU, memory and disk usage of the node and allows user to set up the number of threads (mostly depends on the disk space of the node, no more than 4) and memory (16GB RAM: set 20-30, 64GB RAM: set 100, higher: 200) to use on that node.

The process will generate files named 'probes,cn,chr(1-23).tsv', 'fracB,chr(1-23).tsv', 'fracBseg.tab' and 'segments,fracb.tsv' in 'hg19/processed' folder and sync them to '130.60.23.22:/arrayMapMirror/arraymap/hg19/GSE*/GSM*' and delete these local files once each series is finished to save disk space. Finally, 'concatenate22.sh' command will merge chromosome files to 'probes,cn.tsv' and 'probes,fracb.tsv' files directly on 23.22.

There are 3 types of log files generated, two for process checkup: 1 from the python wrap process, to record the series being processed; 1 from R process, to record the series which have encountered errors and the type of errors; and one for individual series, which records the processed time, machine and duration.

Note: all folder directories are prefixed with '/Users/bgprocess/aroma'
