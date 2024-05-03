library(ggplot2)
library(tidyr)
library(dplyr)
library(tidyverse)
library(data.table)
library(stringr)
library(ggpubr)
library(relaimpo)
library(caret)
library(randomForest)
library(varImp)
library(vegan)
library(pheatmap)
library(patchwork)
library(data.table)
library(caret)
library(dplyr)
library(coin)
library(tidyr)


# load data abs per ko per sample
data = read.table("kpath_number_per_sample.tab",sep = '\t',header = F,quote = "")
names(data) = c("ko","sample","count","desc")
data$sample = tstrsplit(data$sample, "_")[[2]]
data$sample = str_replace_all(data$sample,"b",'')

# load marker copy number per samplpe
mabs = read.table("../../maker_gene_count/maker_count_per_sample.tab", header = F, sep = "\t")
names(mabs) = c("sample","mean_n_mar","median_n_mar","total_n_mar")
mabs$sample = tstrsplit(mabs$sample, "_")[[2]]
mabs$sample = str_replace_all(mabs$sample,"b",'')
mabs = mabs %>% 
  group_by(sample) %>% 
  dplyr::select(sample,mean_n_mar) %>% 
  unique()

data = data %>% 
  left_join(mabs,by = "sample")


# load metadata
metadata = read.table("metadata.tab",header = T,sep = '\t')
metadata$sample = paste(metadata$Tube.number,sep  = '')
data = merge(data,metadata,by = "sample")

# copy number per sample
data$copy_n_per_sample = data$count / data$mean_n_mar

# select columns
data = data %>% 
  dplyr::select(ko,copy_n_per_sample,remark,sample)


lvs = levels(as.factor(data$remark))
lvs = lvs[ !lvs == 'control']
i = "antibiotics"
for (i in lvs) {
  
  print(i)

  dspread = data %>%
    filter(remark %in% c("control",i)) %>% 
    spread(key = ko,value = copy_n_per_sample) 
  
  
  #  Wilcoxon test, fdr multitesting correction 
  mean_comparison = function(x) p.adjust((wilcox.test(x ~ remark, data=dspread)$p.value),method = 'fdr',n = ncol(dspread)-2)
  
  p_values = lapply(dspread[3:(ncol(dspread))],mean_comparison)
  p_values = as.matrix(p_values)
  p_values = as.data.frame(p_values)
  p_values$ko = rownames(p_values)
  names(p_values) = c("pvalue","ko")
  p_values$treat = i
  
  name = paste("DA_results_copyN/",i,"_vs_control.tab",sep = '')
  write.table(apply(p_values,2,as.character),name,row.names = F,sep = '\t')
  
  
}


# fold change 
foldch_res = data.frame(cluster = character(),foldc = numeric(),remark = character(),stringsAsFactors =FALSE)
for (i in lvs) {
  
  print(i)

  dspread = data %>%
    filter(remark %in% c("control",i)) %>% 
    spread(key = ko,value = copy_n_per_sample) 
  head(dspread[1:4])
  
  for (j in c(3:(ncol(dspread)))){
    abs_j = filter(dspread,remark == i)[,j]
    abs_ctr = filter(dspread,remark == "control")[,j]
    qj <- quantile(log10(abs_j), probs=0.5)
    qctr <- quantile(log10(abs_ctr), probs=0.5)
    change = sum(qj - qctr)/length(qctr)
    new_row = c(as.character(names(dspread)[j]),as.numeric(change),as.character(i))
    foldch_res[(nrow(foldch_res) + 1), ] = new_row
  }
}


write.table(foldch_res,"fold_change_kpaths_copyN.tab",row.names = F,sep = '\t')

