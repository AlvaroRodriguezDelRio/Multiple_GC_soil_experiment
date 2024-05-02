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

setwd("~/analysis/Berlin/soil 2019 metagenomes/gene_prediction/")

###
# heatmap
###


data_d = read.table("kpaths/heatmap_kpath_general_data_plot.tab")
data_d = data_d[ order(as.numeric(tstrsplit(rownames(data_d), " ")[[1]])), ]
names(data_d) = str_replace_all(names(data_d),"\\."," ")

head(data_d)
color <- colorRampPalette((c("skyblue", "white", "orange")))(50)
p_heat = pheatmap(data_d, cluster_rows = T, show_rownames = T,scale = "column",show_colnames = T,
         labels_row=rownames(data_d),color = color,border_color = FALSE,
         angle_col = 45,cutree_rows = 2,main = "Figure S5")# ,annotation_col=BarMeta, ,fontsize_row = 5


p_heat
pdf("~/analysis/Berlin/soil 2019 metagenomes/Figures/Figure SX heatmap.pdf")
p_heat
graphics.off()



#####
# heatmap not collapsed  
####

data = read.table("kpaths/data_heatmap_kpath_all.tab")

data_d = data %>% 
  filter(sample != 85) %>% 
  filter(!grepl("Cancer",general)) %>% 
  filter(!grepl("Chronology",general)) %>% 
  filter(!grepl("Target-based classification",general)) %>% 
  filter(!grepl("Immune",general)) %>% 
  filter(!grepl("Chemical structure transformation maps",general)) %>% 
  filter(!grepl("Digestive",general)) %>% 
  filter(!grepl("Development",general)) %>% 
  filter(!grepl("Global",general)) %>% 
  filter(!grepl("Endocrine",general)) %>% 
  filter(!grepl("Circulatory",general)) %>% 
  filter(!grepl("Sensory",general)) %>% 
  filter(!grepl("Nervous",general)) %>% 
  filter(!grepl("Cardiovascular",general)) %>% 
  filter(!grepl("dependence",general)) %>% 
  filter(!grepl("NA",general)) %>% 
  filter(!grepl("Aging",general)) %>%  
  mutate(cpc = n_genes / mean_n_mar) %>% 
  group_by(sample) %>% 
  dplyr::select(kpath,cpc,sample) %>% 
  spread(kpath,cpc) 

head(data_d)

dheat = as.data.frame(data_d)
rownames(dheat) = dheat$sample

d_heat_n = merge(dheat,metadata,by = 'sample')

row_annot = data.frame(as.character(d_heat_n$Lv),as.character(d_heat_n$remark))
rownames(row_annot) = rownames(dheat)
names(row_annot) = c("Lv","Treatment")
row_annot$Treatment  = factor(row_annot$Treatment,levels = c("control","antibiotics", "copper", "drought",
                                                             "microplastic", "Ndep", "salinity", "temp",
                                                             "fungicide", "glyphosate","insecticide","Level 8"
))

col_annot = data %>% 
  dplyr::select(kpath,general) %>% 
  unique()

rownames(col_annot) = col_annot$kpath
col_annot$kpath = NULL
names(col_annot) = c("Category")
dheat$sample = NULL

p_heat2 = pheatmap(dheat , annotation_row = row_annot,cutree_rows = 3,show_rownames = T,
         scale = "column",show_colnames = F)#annotation_col = col_annot



####
# pcoa 
####

pcoa_plot = read.table("kos/kos_PCOA_copy_per_cell.plot.tab")
head(pcoa_plot)
pcoa_plot$Treatment = as.factor(pcoa_plot$remark)

pcoa_plot$Treatment = factor(pcoa_plot$Treatment,levels = c("control","antibiotics", "copper", "drought",
                                            "microplastic", "Ndep", "salinity", "temp",
                                            "fungicide", "glyphosate","insecticide","Level 8"
))

custom_palette <- c("#999999", "#377EB8", "#F781BF","#4DAF4A", "#FF7F00",
                    "#E41A1C", "#A65628", "#984EA3", "#FFFF33", "#66C2A5", "black", "#8EBA42")



p_pcoa = ggplot(pcoa_plot)+
  geom_point(aes(x = Axis.1, y = Axis.2,color = Treatment,shape = as.factor(Lv)),size =3)+
  theme_classic()+
  xlab(paste("PCoA 1 (35.92 %)",sep = ""))+
  ylab(paste("PCoA 2 (9.02 %)",sep = ""))+
  geom_text_repel(data = pcoa_plot[pcoa_plot$sample %in% c(124,127,128),],
                  aes(label = sample,x = Axis.1,y = Axis.2),
                  box.padding   = 1.5, 
                  point.padding = 0.5,
                  segment.color = 'grey50')+
  labs(shape = "Number of factors",color = "Treatment")+
  theme(legend.position = "none")+
  scale_color_manual(values = custom_palette)


p_pcoa

###
# DA plots
###

# from kpath_plot_DA_kpaths_as_CRC.r
# q-value < 0.01

data = read.table("kpaths/kpath_DA_data_plot.tab")
nrow(data)
data = data %>%
  arrange(desc(fold_change))
y_order = unique(data$desc_)
data$desc_ =  factor(data$desc_, levels = y_order)

#cpalette <- brewer.pal(12, "Set3")

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

p_DA_cols

#####
# correlations
####

# load data ko number genes per sample
data_c = read.table("KOs/ko_number_per_sample.desc.tab",header = F, sep = '\t',quote = "")
head(data_c)
names(data_c) = c("desc","ko","sample","n_genes")
data_c$sample = tstrsplit(data_c$sample, "_")[[2]]
data_c$sample = str_replace_all(data_c$sample,"b",'')
head(data_c)

data = merge(data_c,metadata,by = "sample")
data = merge(data,mabs,by = "sample")
data$total_n_reads = 31310395
head(data)
data$remark = factor(data$remark,levels = c("control","antibiotics", "copper", "drought",
                                      "microplastic", "Ndep", "salinity", "temp",
                                      "fungicide", "glyphosate","insecticide","Level 8"
))



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

correlation_coef <- d %>%
  ungroup() %>% 
  summarize(correlation = cor(WSA, n,method = "spearman"),
            p_value = cor.test(WSA, n)$p.value)


correlation_coef
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


cor.test(d$CO2_ppm_ThirdWeek,d$n)

correlation_coef <- d %>%
  ungroup() %>% 
  summarize(correlation = cor(CO2_ppm_ThirdWeek, n),
            p_value = cor.test(CO2_ppm_ThirdWeek, n)$p.value)

correlation_coef
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

corr_resp


# Soil aggregation genes / correlation with WSA 
head(data)
d = data %>% 
  filter(ko %in% c("K01991", # wza, gfcE; polysaccharide biosynthesis/export protein, https://www.biorxiv.org/content/10.1101/860932v2.full
                   "K09688", # capsular polysaccharide export system (Kps) https://www.biorxiv.org/content/10.1101/860932v2.full
                   "K09689", # capsular polysaccharide export system (Kps) https://www.biorxiv.org/content/10.1101/860932v2.full
                   "K10107",# capsular polysaccharide export system (Kps) https://www.biorxiv.org/content/10.1101/860932v2.full
                   "K07091",# lipopolysaccharide export proteins (Lpt)
                   "K09774", # lipopolysaccharide export proteins (Lpt)
                   "K11719", # lipopolysaccharide export proteins (Lpt)
                   "K11720", # lipopolysaccharide export proteins (Lpt)
                   "K04077"#HSP60/GroEL
  )) 

d$remark = factor(d$remark,levels = c("control","antibiotics", "copper", "drought",
                                      "microplastic", "Ndep", "salinity", "temp",
                                      "fungicide", "glyphosate","insecticide","Level 8"
))

d$desc =   substr(d$desc, 1, 4)

p_corr_WSA_genes = ggplot(d,aes(x = WSA,y = n_genes / mean_n_mar))+
  theme_classic()+
  geom_smooth(method = "lm")+
  geom_point(aes(color = as.factor(Lv)))+
  theme(axis.text.x=element_text(angle=90, hjust=1))+
  ylab("Copy number per cell")+
  facet_wrap(~desc,scales = "free")+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2")+
  labs(fill = "Number of factors",color = "Number of factors")+
  ggtitle("Figure S6")


p_corr_WSA_genes
pdf("~/analysis/Berlin/soil 2019 metagenomes/Figures/Figure SX WSA genes correlation.pdf")
p_corr_WSA_genes
graphics.off()

####
# combine 
###


layout <- "AACCC
           DDDDD"

corr_com = corr_flg  / corr_resp 

pdf("~/analysis/Berlin/soil 2019 metagenomes/Figures/Figure 4 functional genes.pdf", width=10, height=8)
free(p_pcoa) + corr_com  + p_DA_cols + plot_layout(design = layout,guides = "collect")+
  plot_annotation(tag_levels = 'A',title = 'Figure 4')
graphics.off()

#######

layout <- "AAAACCC"
corr_com = corr_flg  / corr_resp + plot_layout(guides = "collect")

pdf("~/analysis/Berlin/soil 2019 metagenomes/Figures/Figure 4 functional genes.pdf", width=10, height=5)
free(p_pcoa) + corr_com  + plot_layout(design = layout,guides = "collect")+
  plot_annotation(tag_levels = 'A')
graphics.off()


#####

corr_com = (corr_flg + corr_resp) + plot_layout(guides = "collect")
layout <- "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAABBBBBBBBB
           AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAABBBBBBBBB
           AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAABBBBBBBBB
           CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCDDDDDDDDDDDDDDD
           CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCDDDDDDDDDDDDDDD"


as_ggplot(p_heat$gtable) + p_DA  + corr_com + free(p_pcoa) + plot_layout(design = layout)+
  plot_annotation(tag_levels = 'A')


pdf("~/analysis/Berlin/soil 2019 metagenomes/Figures/Figure 4 functional genes.pdf", width=10, height=8)
as_ggplot(p_heat$gtable) + p_DA  + corr_com + free(p_pcoa) + plot_layout(design = layout)+
  plot_annotation(tag_levels = 'A')
graphics.off()



layout <- "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAABBBBBBBBBBBBBBBBBBBBB
           AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAABBBBBBBBBBBBBBBBBBBBB
           AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAABBBBBBBBBBBBBBBBBBBBB
           CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCDDDDDDDDDDDDDDDDDDD
           CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCDDDDDDDDDDDDDDDDDDD"


as_ggplot(p_heat$gtable) + free(p_pcoa) + p_DA + corr_com + plot_layout(design = layout)+
  plot_annotation(tag_levels = 'A')




