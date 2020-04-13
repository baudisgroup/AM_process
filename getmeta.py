#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 13 15:12:57 2017

@author: pgweb
"""

import os
import sys
normalwords = ("reference", "normal", "healthy", "unaffected","germline","control","peripheral","remission","non-malignant","Healthy", "Normal", "Reference", "Unaffected","Control","Peripheral","Remission")
nameline = ("title","source","characteristics")

#path = "/Users/pgweb/aroma/GEOmeta/"+series
series = sys.argv[1]
platformNeed = sys.argv[2]

def platform2GPL(platform):
    map = {'Mapping10k_Xba142':'GPL2641',
    'Mapping50K_Hind240':'GPL2004',
    'Mapping50K_Xba240':'GPL2005',
    'Mapping250K_Nsp':'GPL3718',
    'Mapping250K_Sty':'GPL3720',
    'CytoScan750k':'GPL18637',
    'GenomeWideSNP_5':'GPL6804',
    'GenomeWideSNP_6':'GPL6801',
    'CytoScanHD':'GPL16131'}
    return map.get(platform)
platformNeed = platform2GPL(platformNeed)

path = "/Volumes/arraymapIncoming/GEOmeta/"+series
output = []
for _,subdir,_ in os.walk(path):
    subdir.sort()
    platform = []
    for filepath in subdir:
        for file in os.listdir(os.path.join(path,filepath)):
            filepath = os.path.join(path,filepath,file)
            if filepath.endswith(".soft") == 1:
                c = 0
                tmp = []
                with open(filepath,encoding='utf-8') as f:
                    content = f.readlines()
                content = [x.strip() for x in content]
                for line in content:
                    info =  line.split()
                    info2 = info[0].split("_")
                    if len(info2) >1: info2 = info2[1]
                    if line.startswith("!Sample_platform_id = GPL"): platform.append(line.split(" = ")[1])
                    if info2 in nameline:
                        tmp.append(line)
                        if any([i in line.split() for i in normalwords]):
                            output.append("Normal")
                            c += 1
                            break
                if c == 0:
                    output.append("Tumor")
idx=[i for i,p in enumerate(platform) if p == platformNeed]
output = output[idx]
print(output)
