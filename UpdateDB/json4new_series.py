import subprocess as sbp
import os, sys

new_serarr = 'new_serarr_20191120.txt'
with open(new_serarr) as f:
    for l in f:
        i,j = l.strip().split()
        sbp.run('perl pgx_segfileconverter.pl  \
            -f /Volumes/arrayMaster/arraymap/grch38/{0}/{1}/segments,cn.tsv \
            -outdir /Volumes/arrayMaster/arraymap/grch38/{0}/{1}'.format(i,j), shell=True)
