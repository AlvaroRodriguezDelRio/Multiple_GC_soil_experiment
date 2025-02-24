library(ggplot2)
library(tidyr)
#install.packages("compositions")
library(dplyr)
library(patchwork)
library(vegan)
#install.packages("ecodist")
library(ecodist)
#install.packages("MASS")
library(MASS)
library(data.table)
library(stringr)
#install.packages("ggsignif")
library(ggsignif)
#library(devtools) # Load the devtools package
#install_github("microbiome/microbiome") # Install the package
library(microbiome)
library(ggpubr)
library(relaimpo)
library(randomForest)
library(ggrepel)
library(forcats)



setwd("~/analysis/Berlin/soil 2019 metagenomes/Figures/Source_data/")

# 4A

d = read.csv('4A.csv')
d$remark = factor(d$remark,levels = c("Control","Warming","Drought","Nitrogen deposition",
                                      "Salinity","Heavy metals","Microplastics",
                                      "Antibiotics","Fungicide", "Herbicide",
                                      "Insecticide","8 factors"
))

# plot profiled flipped 
p_prof = ggplot(d)+
  theme_void()+
  geom_histogram(aes(y = sample,x = rel_sp *100 / tot,fill = class),stat = "identity")+
  theme(text = element_text(size = 25))+
  theme(axis.text.x=element_text(angle=90, hjust=1,size = 15),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "bottom",
        strip.text.y.left = element_text(angle = 0),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20))+
  facet_wrap(remark~.,scales = 'free_y',ncol = 1,strip.position = "right")+
  scale_fill_brewer(palette="Spectral",name = "")+
  ylab("Relative abundance (%)")+
  xlab("")+
  scale_y_discrete(position = "left")+
  guides(fill=guide_legend(nrow=2, byrow=TRUE))


p_prof


# 4B

div = read.csv('4B.csv')

control = div %>% 
  filter(Lv == 0)
control_median_rich = median(control$shannon)

df.summary <- div %>%
  group_by(remark,Lv) %>%
  summarise(
    sd = sd(shannon, na.rm = TRUE),
    len = median(shannon)
  )
df.summary


df.summary$Lv = as.factor(df.summary$Lv)

p2 = ggplot(div)+
  geom_jitter(aes(y = shannon, x = fct_rev(remark)),position = position_jitter(0.2), color = "darkgray",alpha = 0.2) + 
  geom_errorbar(aes(ymin = len-sd, ymax = len+sd,x = fct_rev(remark),color = Lv),data = df.summary,
                position = position_dodge(0.3), width = 0.2)+
  geom_point(aes(x = fct_rev(remark), y = len,color = Lv), data = df.summary,
             position = position_dodge(0.3)) +
  theme_minimal()+
  xlab("")+
  stat_compare_means(aes(group=fct_rev(remark), x = fct_rev(remark), y = shannon),label = "p.signif", method = "wilcox.test",
                     ref.group = "Control",hide.ns = TRUE, tip.length = 0)+
  coord_flip()+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2")+
  geom_hline(yintercept=control_median_rich, linetype="dashed", 
             color = "grey", size=1)+
  theme(axis.text.y = element_blank(),   # Remove x-axis ticks
        strip.background = element_blank(),  # Remove facet labels
        strip.text = element_blank(),
        legend.position = "none")+  # Remove facet labels
  xlab("")+
  ylab("Shannon \ndiversity")+
  theme(text = element_text(size = 25))+
  scale_y_continuous(breaks = c(4.5,5.5))



p2 

# 4C

pcoa_plot = read.csv('4C.csv')

custom_palette <- c("#999999", "#377EB8", "#F781BF","#4DAF4A", "#FF7F00",
                    "#E41A1C", "#A65628", "#984EA3", "#FFFF33", "#66C2A5", "#8EBA42","black")

pcoa_plot$remark = factor(pcoa_plot$remark,levels = c("Control","Warming","Drought","Nitrogen deposition",
                                                      "Salinity","Heavy metals","Microplastics",
                                                      "Antibiotics","Fungicide", "Herbicide",
                                                      "Insecticide","8 factors"
))

p_cor = ggplot(pcoa_plot)+
  geom_point(aes(x = V1, y = V2,color = remark),size = 4)+
  theme_classic()+
  xlab('PCoA 1 (47.34%)')+
  ylab('PCoA 2 (14.68%)')+
  #xlab(paste("PCoA 1 (",round(pcoa$values$Relative_eig[1]*100,digits = 2),"%)",sep = ""))+
  #ylab(paste("PCoA 2 (",round(pcoa$values$Relative_eig[2]*100,digits = 2),"%)",sep = ""))+
  #  theme(text = element_text(size = 20))+
  geom_text_repel(data = pcoa_plot[pcoa_plot$s %in% c("124 Level 8","127 Level 8","128 Level 8"),],
                  aes(label = sample,x = V1,y = V2),
                  box.padding   = 3, 
                  point.padding = 0.5,
                  segment.color = 'grey50')+
  theme(legend.position = "bottom", 
        legend.text = element_text(size=20),
        legend.title = element_blank())+
  guides(shape = FALSE,fill=guide_legend(nrow=2, byrow=TRUE),
         color = guide_legend(nrow = 3, byrow = TRUE))+ 
  scale_color_manual(values = custom_palette)+
  theme(text = element_text(size = 25))


p_cor


layout <- "
AAABCCC
"
p_prof + p2 + p_cor +  plot_layout(design = layout)+plot_annotation(tag_levels = 'A')# , title = 'Figure 3'

