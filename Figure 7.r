library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)
library(data.table)
library(stringr)
library(ggpubr)
library(ggrepel)


setwd('~/analysis/Berlin/soil 2019 metagenomes/Figures/Source_data/')

# 7B

custom_palette <- c("#999999", "#377EB8", "#F781BF","#4DAF4A", "#FF7F00",
                    "#E41A1C", "#A65628", "#984EA3", "#FFFF33", "#66C2A5","#8EBA42","black")



data_pcoa = read.csv('7B.csv')
pcoa_plot = ggplot(data_pcoa)+
  geom_point(aes(x = V1, y = V2,color = remark),size = 5)+
  theme_classic()+
  xlab('PCoA 1(21.06%)')+
  ylab('PCoA 2 (13.31%)')+
  #xlab(paste("PCoA 1 (",round(pcoa$values$Relative_eig[1]*100,digits = 2),"%)",sep = ""))+
  #ylab(paste("PCoA 2 (",round(pcoa$values$Relative_eig[2]*100,digits = 2),"%)",sep = ""))+
  # theme(text = element_text(size = 20))+
  geom_text_repel(data = data_pcoa[data_pcoa$s %in% c("124 Level 8","127 Level 8","128 Level 8"),],
                  aes(label = sample,x = V1,y = V2),
                  box.padding   = 3, 
                  point.padding = 0.5,
                  segment.color = 'grey50')+
  scale_color_manual(values = custom_palette,name = "Treatment")+
  theme(text = element_text(size = 25))



pcoa_plot

# 7A
d = read.csv('7A.csv')

p1 = ggplot(d)+
  geom_histogram(aes(x = log10(size), fill = novel),bins = 10)+
  theme_classic()+
  scale_fill_manual(values = c("red","grey"))+
  xlab("Log10 (# genes per family)")+
  theme(legend.position = "none")+
  theme(text = element_text(size = 25))


d1 = d %>% 
  group_by(novel) %>% 
  mutate(tot = sum(size)) %>% 
  dplyr::select(tot,novel) %>% 
  unique()

p2 = ggplot(d1)+
  geom_histogram(aes(x = tot,y = novel,fill = novel),position = "dodge",stat = "identity")+
  ylab("")+
  xlab("# genes")+
  theme_bw()+
  scale_fill_manual(values = c("red","grey"))+
  theme(legend.position = "none")+
  scale_x_continuous(breaks = c(0,10000000,20000000))+
  theme(text = element_text(size = 20))

p2
options(scipen = 999)
p_num_novel = p1 + inset_element(p2,  0.4, 0.6, 1, 0.9)
p_num_novel


top_plot      <- wrap_elements((p_num_novel) + 
                                 plot_annotation(title = "A"))
bottom_plot   <- wrap_elements(pcoa_plot + 
                                 plot_annotation(title = "B"))
top_plot + bottom_plot 

