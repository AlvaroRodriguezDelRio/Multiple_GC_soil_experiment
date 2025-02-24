library(ggplot2)
library(tidyr)
library(dplyr)
library(tidyverse)
library(data.table)
library(stringr)
library(stringr)
library(ggpubr)
library(patchwork)
library(scales)
library(vegan)
library(ggrepel)

setwd("~/analysis/Berlin/soil 2019 metagenomes/Figures/Source_data/")

d = read.csv('6A.csv')

head(d)


p_general = ggplot(d %>% 
                     filter(!cat %in% c("antibiotic target alteration","antibiotic efflux","reduced permeability to antibiotic","antibiotic target alteration;antibiotic target replacement")),
                   aes(x = remark, y = n))+
  geom_boxplot(aes(color = as.factor(Lv),fill = as.factor(Lv)),outlier.shape = NA,alpha = 0.3)+
  coord_flip()+
  facet_wrap(~cat,scales = "free_x",ncol = 2)+
  stat_compare_means(label = "p.signif", method = "wilcox.test",
                     ref.group = "Control",hide.ns = TRUE, tip.length = 0)+
  geom_jitter(aes(color = as.factor(Lv),fill = as.factor(Lv)),shape=16, position=position_jitter(0.2),alpha = 0.8)+
  ylab("Copies per cell")+
  theme_minimal()+
  geom_hline(aes(yintercept= median_control), linetype="dashed", 
             color = "grey", size=1)+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2")+
  labs(fill = "Number of factors",color = "Number of factors")+
  xlab("")+
  theme(text = element_text(size = 25))+
  scale_y_continuous(breaks = pretty_breaks(n = 2))


p_general



# 6B

d = read.csv('6B.csv')

p_genes = ggplot(d %>% 
                   filter(name %in% c("BJP-1","RbpA")) %>% 
                   filter(sample != 85),aes(x = remark, y = ncard / mean_n_mar))+
  geom_boxplot(aes(color = as.factor(Lv),fill = as.factor(Lv)),outlier.shape = NA,alpha = 0.3)+
  coord_flip()+
  facet_wrap(~name,scales = "free_x",ncol = 2)+
  stat_compare_means(label = "p.signif", method = "wilcox.test",
                     ref.group = "Control",hide.ns = TRUE, tip.length = 0)+
  theme_classic()+
  geom_jitter(aes(color = as.factor(Lv),fill = as.factor(Lv)),shape=16, position=position_jitter(0.2),alpha = 0.8)+
  ylab("Copies per cell")+
  theme_minimal()+
  geom_hline(aes(yintercept= median_control), linetype="dashed", 
             color = "grey", size=1)+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2")+
  labs(fill = "Number of factors",color = "Number of factors")+
  xlab("")+
  theme(text = element_text(size = 25))+
  scale_y_continuous(breaks = pretty_breaks(n = 2))



pbox = p_general / p_genes 
pbox


# 6C

custom_palette <- c("#999999", "#377EB8", "#F781BF","#4DAF4A", "#FF7F00",
                    "#E41A1C", "#A65628", "#984EA3", "#FFFF33", "#66C2A5","#8EBA42","black")



data = read.csv('6C.csv')

pcorr = ggplot(data ,aes(x = ncard / mean_n_mar, y = V1))+
  geom_point(alpha = 0.5,size = 4,aes(color = remark))+
  geom_smooth(method = 'lm')+
  theme_minimal()+
  ylab("Bacterial composition")+
  xlab("ARG copy number per cell")+
  scale_color_manual(values = custom_palette,name = "Treatment")+
  theme(text = element_text(size = 25))

pcorr

#6D

cont =read.csv('6D.csv')

corr_plasmids = ggplot(cont,aes(x = ncard / mean_n_mar, y =  nplasmids_strict / n_contigs))+
  geom_point(alpha = 0.5,size = 4,aes(color = remark))+
  geom_smooth(method = 'lm')+
  theme_minimal()+
  ylab("Plasmid frequency")+
  xlab("ARG copy number per cell")+
  scale_color_manual(values = custom_palette,name = "Treatment")+
  theme(text = element_text(size = 25))

corr_plasmids

# 6E

corr_phages = ggplot(cont,aes(x = ncard / mean_n_mar, y =  nphage_strict / n_contigs))+
  geom_point(alpha = 0.5,size = 4,aes(color = remark))+
  geom_smooth(method = 'lm')+
  theme_minimal()+
  ylab("Phage frequency")+
  xlab("ARG copy number per cell")+
  scale_color_manual(values = custom_palette,name = "Treatment")+
  theme(text = element_text(size = 25))



layout <- "AACC
           AACC
           AADD
           BBDD
           BBEE
           BBEE"

p_general + p_genes  + pcorr + corr_plasmids + corr_phages + plot_layout(guides = "collect",design = layout)+# design = layout
  plot_annotation(tag_levels = 'A')# ,title = 'Figure 5'
