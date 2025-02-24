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

setwd("~/analysis/Berlin/soil 2019 metagenomes/Figures/Source_data/")


data = read.csv('5C.csv')

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
  #  scale_color_brewer(palette = "Paired")+
  scale_color_manual(values = c(
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b",
    "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#98df8a", "#ff9896"
  )) +
  labs(color="General pathway classification")+
  guides(color = guide_legend(override.aes = list(size = 4)))+   # Adjust the size of the points in the legend
  theme(text = element_text(size = 25))


p_DA_cols

# 5A

d = read.csv('5A.csv')
head(d)

corr_flg = ggplot(d,aes(x = WSA,y = n))+
  theme_classic()+
  geom_smooth(method = "lm")+
  geom_point(size = 5,alpha = 0.5)+ # aes(shape = as.factor(Lv),color = remark)
  theme(axis.text.x=element_text(angle=90, hjust=1))+
  ylab("Copy # per cell")+
  labs(color = "Treatment",shape = "Number of factors")+
  ggtitle("Flg operon")+
  guides(color = guide_legend(override.aes = list(size = 4)),
         size = guide_legend(override.aes = list(size = 4)))+
  theme(text = element_text(size = 25))


corr_flg

# 5B

d = read.csv('5B.csv')
corr_resp = ggplot(d,aes(x = CO2_ppm_ThirdWeek,y = n))+
  theme_classic()+
  geom_smooth(method = "lm")+
  geom_point(size = 5,alpha = 0.5)+ #aes(shape = as.factor(Lv),color = remark)
  theme(axis.text.x=element_text(angle=90, hjust=1))+
  ylab("Copy # per cell")+
  labs(color = "Treatment",shape = "Number of factors")+
  ggtitle("cox operon")+
  xlab("CO2 ppm")+
  guides(color = guide_legend(override.aes = list(size = 4)),
         size = guide_legend(override.aes = list(size = 4)))+
  theme(text = element_text(size = 25))
# Adjust the size of the points in the legend


corr_resp


layout <- "ACC
           BCC"


corr_com = free(corr_flg)  + free(corr_resp)


corr_com +  p_DA_cols + plot_layout(design = layout,guides = "collect")+
  plot_annotation(tag_levels = 'A')



