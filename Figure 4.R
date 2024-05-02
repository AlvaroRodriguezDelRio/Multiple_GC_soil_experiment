# funct combined plot 
library(ggplot2)
library(patchwork)
library(dplyr)
library(data.table)
library(tidyr)
library(tidyverse)
library(pheatmap)
library(ggrepel)
library(ggalt)
library(ggpubr)


####
# composition plot
####


# load metadata
metadata = read.table("metadata.tab",header = T,sep = '\t')
metadata$Tube.number
metadata$sample = paste(metadata$Tube.number,sep  = '')

# load ko number per sample
data = read.table("KOs/ko_number_per_sample.desc.tab",header = F, sep = '\t',quote = "")
names(data) = c("desc","ko","sample","n_genes")
data$sample = tstrsplit(data$sample, "_")[[2]]
data$sample = str_replace_all(data$sample,"b",'')


data = data %>% 
  left_join(metadata,by = "sample")

# load marker copy number per samplpe
mabs = read.table("../maker_gene_count/maker_count_per_sample.tab", header = F, sep = "\t")
names(mabs) = c("sample","mean_n_mar","median_n_mar","total_n_mar")
mabs$sample = tstrsplit(mabs$sample, "_")[[2]]
mabs$sample = str_replace_all(mabs$sample,"b",'')
mabs = mabs %>% 
  group_by(sample) %>% 
  dplyr::select(sample,mean_n_mar) %>% 
  unique()

data = data %>% 
  left_join(mabs,by = "sample")


data$copy_n_per_cell = data$n_genes / data$mean_n_mar

# select columns
data = data %>% 
  dplyr::select(ko,copy_n_per_cell,remark,sample) %>% 
  filter(sample != 85) %>% 
  unique()

dspread = data %>%
  spread(key = ko,value = copy_n_per_cell) 

rownames(dspread) = dspread$sample
dspread$sample = NULL
dspread$remark = NULL
dspread = (as.matrix(dspread))

d_m = vegdist(dspread, method="bray")
pcoa = ape::pcoa(d_m)
pcoa_plot = as.data.frame(pcoa$vectors[,1:2])

pcoa_plot$sample = row.names(pcoa_plot)
pcoa_plot = pcoa_plot %>% 
  left_join(metadata,by = "sample")

custom_palette <- c("#999999", "#377EB8", "#F781BF","#4DAF4A", "#FF7F00",
                    "#E41A1C", "#A65628", "#984EA3", "#FFFF33", "#66C2A5", "black", "#8EBA42")



p_pcoa = ggplot(pcoa_plot)+
  geom_point(aes(x = Axis.1, y = Axis.2,color = Treatment,shape = as.factor(Lv)),size =3)+
  theme_classic()+
  xlab(paste("PCoA 1 (",round(pcoa$values$Relative_eig[1]*100,digits = 2),"%)",sep = ""))+
  ylab(paste("PCoA 2 (",round(pcoa$values$Relative_eig[2]*100,digits = 2),"%)",sep = ""))+
  geom_text_repel(data = pcoa_plot[pcoa_plot$sample %in% c(124,127,128),],
                  aes(label = sample,x = Axis.1,y = Axis.2),
                  box.padding   = 1.5, 
                  point.padding = 0.5,
                  segment.color = 'grey50')+
  labs(shape = "Number of factors",color = "Treatment")+
  theme(legend.position = "none")+
  scale_color_manual(values = custom_palette)

###
# DA plots
###

data = read.table("kpaths/kpath_DA_data_plot.tab")
data = data %>%
  arrange(desc(fold_change))
y_order = unique(data$desc_)
data$desc_ =  factor(data$desc_, levels = y_order)

p_DA_cols = ggplot(data %>% 
                arrange(desc(fold_change))) +
  geom_dumbbell(aes(x=0, xend=fold_change, y = desc_, color = general_),
                size=1.2) +
  theme_minimal()+
  theme(axis.text.x=element_text(angle=90, hjust=1),
        #legend.position = "bottom",
        strip.text.y.left = element_text(angle = 0))+
  xlab("Fold Change")+
  ylab("KEGG pathway")+
  scale_color_brewer(palette = "Paired")+
  labs(color="General pathway classification")

#####
# correlations
####

# load data ko number genes per sample
data_c = read.table("KOs/ko_number_per_sample.desc.tab",header = F, sep = '\t',quote = "")
names(data_c) = c("desc","ko","sample","n_genes")
data_c$sample = tstrsplit(data_c$sample, "_")[[2]]
data_c$sample = str_replace_all(data_c$sample,"b",'')

# load marker copy number per samplpe
mabs = read.table("../../maker_gene_count/maker_count_per_sample.tab", header = F, sep = "\t")
names(mabs) = c("sample","mean_n_mar","median_n_mar","total_n_mar")
mabs$sample = tstrsplit(mabs$sample, "_")[[2]]
mabs$sample = str_replace_all(mabs$sample,"b",'')
mabs = mabs %>% 
  group_by(sample) %>% 
  dplyr::select(sample,mean_n_mar) %>% 
  unique()

data = merge(data_c,metadata,by = "sample")
data = merge(data,mabs,by = "sample")

# flagelum
head(data)
d = data %>% 
  filter(grepl("flg",desc))

d$desc =substr(d$desc, 1, 4)
d = d %>% 
  group_by(sample) %>% 
  mutate(n = sum(n_genes / mean_n_mar)) %>% 
  dplyr::select(WSA,n,Lv,remark) %>% 
  unique()

corr_flg = ggplot(d,aes(x = WSA,y = n))+
  theme_classic()+
  geom_smooth(method = "lm")+
  geom_point(aes(shape = as.factor(Lv),color = remark),size = 2)+
  theme(axis.text.x=element_text(angle=90, hjust=1))+
  ylab("Copy number per cell")+
  labs(color = "Treatment",shape = "Number of factors")+
  ggtitle("Flg operon")+
  scale_color_manual(values = custom_palette)

corr_flg

# respiration
names(data)
d = data %>% 
  filter(sample != 85) %>% 
  filter(ko %in% c("K02274","K02275","K02276","K02277")) %>% 
  group_by(sample) %>% 
  mutate(n = mean(n_genes / mean_n_mar)) %>% 
  dplyr::select(CO2_ppm_ThirdWeek,n,Lv,remark) %>% 
  unique()


corr_resp = ggplot(d,aes(x = CO2_ppm_ThirdWeek,y = n))+
  theme_classic()+
  geom_smooth(method = "lm")+
  geom_point(aes(shape = as.factor(Lv),color = remark),size = 2)+
  theme(axis.text.x=element_text(angle=90, hjust=1))+
  ylab("Copy number per cell")+
  labs(color = "Treatment",shape = "Number of factors")+
  ggtitle("cox operon")+
  xlab("CO2 ppm")+
  scale_color_manual(values = custom_palette)

####
# combine 
###


layout <- "AACCC
           DDDDD"

corr_com = corr_flg  / corr_resp 

pdf("Figure 4.pdf", width=10, height=8)
free(p_pcoa) + corr_com  + p_DA_cols + plot_layout(design = layout,guides = "collect")+
  plot_annotation(tag_levels = 'A',title = 'Figure 4')
graphics.off()
