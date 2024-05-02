library(ggplot2)
library(tidyr)
library(dplyr)
library(tidyverse)
library(data.table)
library(stringr)
library(stringr)
library(ggpubr)
library(patchwork)

setwd("~/analysis/Berlin/soil 2019 metagenomes/gene_prediction/CARD/")

metadata = read.table("../../metadata.tab",header = T,sep = '\t')
metadata$Tube.number
metadata$sample = paste(metadata$Tube.number,sep  = '')


####
# load data
####

# load 
data = read.table("hits.1e7_id80_cov75.n_genes_ass.gene_name.CARD_code.tab",header = F, sep = '\t',quote = "")
names(data) = c("sample","gene","cat","name","n","ncard","ncard_plasmids","nmob")
data$sample = tstrsplit(data$sample, "_")[[2]]
data$sample = str_replace_all(data$sample,"b",'')

data = data %>% 
  left_join(metadata,by = "sample")

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


# load composition per sample
#comp = read.table("../../tax_profile/COMBINED/composition.tab",header = T)
comp = read.table("../../MAGs/semibin/dRep/coverM/composition.no85.tab",header = T)
comp$sample = tstrsplit(comp$sample, " ")[[1]]
comp$sample = str_replace_all(comp$sample,"b",'')
comp = comp %>% 
  dplyr::select(sample,V1,V2)
data = data %>% 
  left_join(comp,by = "sample")

data$remark = factor(data$remark,levels = c("control","antibiotics", "copper", "drought",
                                            "microplastic", "Ndep", "salinity", "temp",
                                            "fungicide", "glyphosate","insecticide","Level 8"
))


#####
# copies per cell, per category general 
####

names(data)
d = data %>% 
  filter(sample != 85) %>% 
  group_by(sample,cat) %>% 
  mutate(n  = sum( ncard / mean_n_mar)) %>% 
  dplyr::select(n,sample,cat,Lv,remark,Salinity..9.) %>% 
  unique()

dcontrol = d %>% 
  filter(remark == 'control') %>%
  group_by(cat) %>% 
  mutate(median_control = median(n)) %>% 
  dplyr::select(cat,median_control) %>% 
  ungroup() %>% 
  unique()


d = d %>% 
  left_join(dcontrol,by = "cat")

p_general = ggplot(d %>% 
                     filter(!cat %in% c("antibiotic target alteration","antibiotic efflux","reduced permeability to antibiotic","antibiotic target alteration;antibiotic target replacement")),
                   aes(x = remark, y = n))+
  geom_boxplot(aes(color = as.factor(Lv),fill = as.factor(Lv)),outlier.shape = NA,alpha = 0.3)+
  coord_flip()+
  facet_wrap(~cat,scales = "free_x",ncol = 2)+
  stat_compare_means(label = "p.signif", method = "wilcox.test",
                     ref.group = "control",hide.ns = TRUE, tip.length = 0)+
  theme_classic()+
  geom_jitter(aes(color = as.factor(Lv),fill = as.factor(Lv)),shape=16, position=position_jitter(0.2),alpha = 0.8)+
  ylab("Copies per cell")+
  theme_minimal()+
  geom_hline(aes(yintercept= median_control), linetype="dashed", 
             color = "grey", size=1)+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2")+
  labs(fill = "Number of factors",color = "Number of factors")

####
# particular genes boxplot
####

dcontrol = data %>% 
  mutate(n  =  ncard / mean_n_mar) %>% 
  filter(remark == 'control') %>%
  group_by(name) %>% 
  mutate(median_control = median(n)) %>% 
  dplyr::select(name,median_control) %>% 
  ungroup() %>% 
  unique()

d = data %>% 
  left_join(dcontrol,by = "name")


p_genes = ggplot(d %>% 
                   filter(name %in% c("BJP-1","RbpA")) %>% 
                   filter(sample != 85),aes(x = remark, y = ncard / mean_n_mar))+
  geom_boxplot(aes(color = as.factor(Lv),fill = as.factor(Lv)),outlier.shape = NA,alpha = 0.3)+
  coord_flip()+
  facet_wrap(~name,scales = "free_x",ncol = 2)+
  stat_compare_means(label = "p.signif", method = "wilcox.test",
                     ref.group = "control",hide.ns = TRUE, tip.length = 0)+
  theme_classic()+
  geom_jitter(aes(color = as.factor(Lv),fill = as.factor(Lv)),shape=16, position=position_jitter(0.2),alpha = 0.8)+
  ylab("Copies per cell")+
  theme_minimal()+
  geom_hline(aes(yintercept= median_control), linetype="dashed", 
             color = "grey", size=1)+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2")+
  labs(fill = "Number of factors",color = "Number of factors")

#####
# corr with composition 
#####


# load collapsed per sample
data = read.table("~/analysis/Berlin/soil 2019 metagenomes/gene_prediction/stats_per_sample.1e7_id80_cov75.tab",header = F, sep = '\t')
names(data) = c("sample","n_genes","ncard","ncard_plasmids","div_card","nres","nres_plasmids","div_res","nmob","div_mob")
data$sample = tstrsplit(data$sample, "_")[[2]]
data$sample = str_replace_all(data$sample,"b",'')

data = data %>% 
  left_join(metadata,by = "sample") %>% 
  left_join(mabs,by = "sample") %>% 
  left_join(comp,by = "sample") %>%
  filter(sample != 85)

custom_palette <- c("#999999", "#377EB8", "#F781BF","#4DAF4A", "#FF7F00",
                    "#E41A1C", "#A65628", "#984EA3", "#FFFF33", "#66C2A5", "black", "#8EBA42")


ggplot(data,aes(y = ncard / mean_n_mar, x = remark))+
  geom_boxplot(aes(color = remark))+
  stat_compare_means(aes(group=remark,y = ncard / mean_n_mar, x = remark),label = "p.signif", method = "t.test",
                     ref.group = "Level 8",hide.ns = TRUE, tip.length = 0)
  
pcorr = ggplot(data ,aes(x = ncard / mean_n_mar, y = V1))+
  geom_point(alpha = 0.5,size = 3,aes(color = remark))+
  geom_smooth(method = 'lm')+
  theme_minimal()+
  ylab("Bacterial composition")+
  xlab("ARG copy number per cell")+
  scale_color_manual(values = custom_palette)


pcorr2 = ggplot(data,aes(x = ncard / mean_n_mar, y = nmob / mean_n_mar))+
  geom_point(alpha = 0.5,size = 3,aes(color = remark))+
  geom_smooth(method = 'lm')+
  theme_minimal()+
  ylab("MGE copy number per cell")+
  xlab("ARG copy number per cell")+
  scale_color_manual(values = custom_palette)


#####
# combine 
#####

layout <- "AACC
           BBDD"

p_general + p_genes  + pcorr + pcorr2 + plot_layout(guides = "collect",design = layout)+# design = layout
  plot_annotation(tag_levels = 'A',title = 'Figure 5')

pdf("~/analysis/Berlin/soil 2019 metagenomes/Figures/Figure 5.pdf", width=15, height=8)
p_general + p_genes  + pcorr + pcorr2 + plot_layout(guides = "collect",design = layout)+# design = layout
  plot_annotation(tag_levels = 'A',title = 'Figure 5')

graphics.off()
