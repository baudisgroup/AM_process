### LOH analysis pipeline
### from segments,fracb.tsv extract segments -> LOH file
### merge LOH file similar to the CN provenance file
library(dplyr)
library(pastecs)
load('model_param.RData')

combine_seg_peak_finder <- function(directory){
  allfracb <- read.table(paste(directory,'probes,fracb.tsv',sep ='/'), header = T)
  fracbseg <- read.table(paste(directory,'segments,fracb.tsv',sep ='/'),  header = T)
  sample_id <- fracbseg[1,1]
  fracbseg <- unique(fracbseg[,c(2,3,4)])
  cnseg <- read.table(paste(directory,'segments,cn.tsv',sep ='/'), header = T)
  cnseg <- cnseg[c(2,3,4)]
  combined_seg <- data.frame('chr'=numeric(),'start'=numeric(), 'end'=numeric())
  for (chr in 1:23) {
    starts = NULL
    ends = NULL
    fracbseg_chr = fracbseg[fracbseg$chr == chr,]
    len_fracbseg_chr = nrow(fracbseg_chr)
    cnseg_chr = cnseg[cnseg$chromosome == chr,]
    len_cnseg_chr = nrow(cnseg_chr)
    starts = sort(unique(c(fracbseg_chr$loc.start, cnseg_chr$start, fracbseg_chr$loc.end[len_fracbseg_chr], cnseg_chr$end[len_cnseg_chr] )))
    len_starts = length(starts)
    combined_seg <- rbind(combined_seg, data.frame('chr'=chr, 'start'= starts[1:len_starts-1],'end'=starts[2:len_starts]))
  }
  
  out = data.frame()
  for (chr in 1:23){
    
    seg <- combined_seg[combined_seg$chr ==chr,]
    fracb <- subset(allfracb, allfracb$CHRO ==chr)
    
    for (j in 1:nrow(seg)){
      range <- c(seg[j,'start'], seg[j,"end"])
      id <- which(fracb$BASEPOS>=range[1] & fracb$BASEPOS< range[2])
      subfracb <- fracb[id,]
      if (nrow(subfracb) <= 5) next
      dens <- density(subfracb$VALUE)
      
      y <- dens$y
      x <- dens$x
      tp<-turnpoints(y)
      peak <- x[tp$peaks]
      peak <- sort(c(peak,1-peak))
      peak <- rm.near.dy(peak)
      if (length(peak) == 0) next
      for (line in peak){
        out = rbind(out, data.frame('id'=sample_id,'chr'=chr,'start'=range[1],'end'=range[2],'fracB'=line))
      }
    }
 
  }
  return (out)
}

output_loh <- function(directory,new_seg_fracb){
  
  seg_ori <- new_seg_fracb
  type_map = c(0,0,0,0,0,1,2,2,3,4)
  names(type_map) = c(6,7,8,9,10,2,1,3,4,5)
  seg_ori <- unique(seg_ori)
  colnames(seg_ori) = c('id','chr','start','end','fracb')
  seg <- seg_ori %>% group_by(chr, start, end) %>% summarise(n = n()) 
  seg$type <- type_map[as.character(seg$n)]
  seg <- merge_same_type(seg)
  loh <- seg[seg$type == 1,] %>% select(chr,start,end)
  loh$length <- loh$end - loh$start
  write.table(loh, paste(directory,'segments,loh.tsv', sep = '/'), sep = '\t', row.names = F, quote = F )
}

output_loh_with_model_validation <- function(directory,new_seg_fracb,logfile){

  seg_ori <- new_seg_fracb
  type_map = c(0,0,0,0,0,1,2,2,3,4)
  names(type_map) = c(6,7,8,9,10,2,1,3,4,5)
  seg_ori <- unique(seg_ori)
  colnames(seg_ori) = c('id','chr','start','end','fracb')
  seg <- seg_ori %>% group_by(chr, start, end) %>% summarise(n = n(), min_fracB = min(fracb)) 
  no_joint_seg <- nrow(seg)
  seg$type <- type_map[as.character(seg$n)]
  pred_map = c('2','4','5')
  for (i in 1:nrow(seg)){
    if (seg$type[i] == 1 & seg$min_fracB[i] > 0.2){
      vals = new_seg_fracb[new_seg_fracb[[2]] == seg[[1]][i] & new_seg_fracb[[3]] <= seg[[2]][i] & new_seg_fracb[[4]] >= seg[[3]][i], 5]
      pred = model_fit(list(d2,d4,d5),vals)
      if (pred != 1) {
        seg$type[i] = type_map[pred_map[pred]]
      }
    }
  }
  seg <- merge_same_type(seg)
  loh <- seg[seg$type == 1,] %>% select(chr,start,end)
  write.table(loh, paste(directory,'segments,loh.tsv', sep = '/'), sep = '\t', row.names = F, quote = F )
  write(sprintf('joined segmentation resulted in %d segments', no_joint_seg), logfile, append = T)
  loh$length <- loh$end - loh$start
  for (chro in unique(loh$chr)) {
    loh_len = sum(loh$length[loh$chr == chro])%/%1000
    write(paste('Chromsome', chro, 'LOH length:', ifelse(loh_len > 1000,  paste(loh_len %/% 1000, 'Mb'), paste(loh_len, 'kb'))), logfile, append = T)
  }
  
}

model_fit <- function(model_list,data){
  
  sum <- NULL
  for (m in model_list){
    sum <- c(sum,sum(sapply(data, function(a) approx(m$x, m$y, a)$y), na.rm = T))
  }
  
  return(which.max(sum))
  
}

merge_same_type <- function(seg){
  merge_idx <- NULL
  last_type <- -1
  last_chr <- 0
  for (i in 1:nrow(seg)){
    if (seg$type[i] == last_type & seg$chr[i] == last_chr) merge_idx = c(merge_idx,i)
    last_type = seg$type[i]
    last_chr = seg$chr[i]
  }
  
  new_seg <- data.frame()
  all_idx <- sort(union(merge_idx, merge_idx-1))
  for (i in 1:nrow(seg)){
    if (i %in% merge_idx){
      new_seg[nrow(new_seg), 'end'] = seg$end[i]
      next
    }
    if (i %in% (merge_idx-1)){
      tmp <- seg[i,]
      tmp$end <- seg$start[i+1]
      new_seg <- rbind(new_seg, as.data.frame(tmp))
      next
    }
    if (!i %in% all_idx) {
      new_seg <- rbind(new_seg, as.data.frame(seg[i,]))
    }
  }
  return(new_seg)
}

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
