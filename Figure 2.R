library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)
library(ggrepel)
library(ggalt)
library(patchwork)
library(ggpubr)

setwd("~/analysis/Berlin/soil 2019 metagenomes/Figures/Source_data/")

d = read.csv("2A.csv")
d$remark = factor(d$remark,levels = c("Control","Warming","Drought","Nitrogen deposition",
                                      "Salinity","Heavy metals","Microplastics",
                                      "Antibiotics","Fungicide", "Herbicide",
                                      "Insecticide","8 factors"
))

# plot profile flipped 
p = ggplot(d)+
  theme_void()+
  geom_histogram(aes(y = sample,x = rel *100 / tot,fill = class),stat = "identity")+
  theme(axis.text.x=element_text(angle=90, hjust=1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "bottom",
        strip.text.y.left = element_text(angle = 0),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        axis.text=element_text(size=14))+
  facet_wrap(remark~.,scales = 'free_y',ncol = 1,strip.position = "right")+
  scale_fill_brewer(palette="Spectral",name = "")+
  ylab("Relative abundance (%)")+
  xlab("")+
  scale_y_discrete(position = "left")+
  guides(fill=guide_legend(nrow=3, byrow=TRUE))+
  theme(text = element_text(size = 25))

p

########
# 2B
#########

div = read.csv('2B.csv')
div$remark = factor(div$remark,levels = c("Control","Warming","Drought","Nitrogen deposition",
                                          "Salinity","Heavy metals","Microplastics",
                                          "Antibiotics","Fungicide", "Herbicide",
                                          "Insecticide","8 factors"
))


control = div %>% 
  filter(Lv == 0)
control_median_rich = median(control$shannon)

div$Lv = as.factor(div$Lv)
div$remark <- factor(div$remark, levels = rev(levels(div$remark)))

df.summary <- div %>%
  group_by(remark,Lv) %>%
  summarise(
    sd = sd(shannon, na.rm = TRUE),
    len = median(shannon)
  )


p2 = ggplot(div)+
  geom_jitter(aes(y = shannon, x = remark),position = position_jitter(0.2), color = "darkgray",alpha = 0.2) + 
  geom_errorbar(aes(ymin = len-sd, ymax = len+sd,x = remark,color = Lv),data = df.summary,
                position = position_dodge(0.3), width = 0.2)+
  geom_point(aes(x = remark, y = len,color = Lv), data = df.summary,
             position = position_dodge(0.3)) +
  theme_minimal()+
  xlab("")+
  stat_compare_means(aes(group=remark, x = remark, y = shannon),label = "p.signif", method = "wilcox.test",
                     ref.group = "Control",hide.ns = TRUE, tip.length = 0,paired = F)+
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
  scale_y_continuous(breaks = c(3,3.5))



p2 

# 2C

data = read.csv('2C.csv')

control = data %>% 
  filter(Lv == 0)
control_median_rich = median(control$n)
data$control_n  = control_median_rich

df.summary <- data %>%
  group_by(remark,Lv) %>%
  summarise(
    sd = sd(n, na.rm = TRUE),
    len = median(n)
  )

data$Lv = as.factor(data$Lv)
df.summary$Lv = as.factor(df.summary$Lv)
p3 = ggplot(data)+
  geom_jitter(aes(y = n, x = remark),position = position_jitter(0.2), color = "darkgray",alpha = 0.2) + 
  geom_errorbar(aes(ymin = len-sd, ymax = len+sd,x = remark,color = Lv),data = df.summary,
                position = position_dodge(0.3), width = 0.2,)+
  geom_point(aes(x = remark, y = len,color = Lv), data = df.summary,
             position = position_dodge(0.3)) +
  theme_minimal()+
  xlab("")+
  stat_compare_means(aes(group=remark, x = remark, y = n),label = "p.signif", method = "wilcox.test",paired = F,
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
  ylab("Genome \nsize (Mbps)")+
  theme(text = element_text(size = 25))+
  scale_y_continuous(breaks = c(5,10))



p3

# 2D

pcoa_plot = read.csv('2D.csv')
pcoa_plot$remark = factor(pcoa_plot$remark,levels = c("Control","Warming","Drought","Nitrogen deposition",
                                                      "Salinity","Heavy metals","Microplastics",
                                                      "Antibiotics","Fungicide", "Herbicide",
                                                      "Insecticide","8 factors"
))


custom_palette <- c("#999999", "#377EB8", "#F781BF","#4DAF4A", "#FF7F00",
                    "#E41A1C", "#A65628", "#984EA3", "#FFFF33", "#66C2A5", "#8EBA42","black")

head(pcoa_plot)
p4 = ggplot(pcoa_plot)+
  geom_point(aes(x = V1, y = V2,color = remark),size = 4)+
  theme_classic()+
  xlab('PCoA 1 (42.43%)')+
  ylab('PCoA 2 (17.66%)')+
  #xlab(paste("PCoA 1 (",round(pcoa$values$Relative_eig[1]*100,digits = 2),"%)",sep = ""))+
  #ylab(paste("PCoA 2 (",round(pcoa$values$Relative_eig[2]*100,digits = 2),"%)",sep = ""))+
  geom_text_repel(data = pcoa_plot[pcoa_plot$sample %in% c(124,127,128),],
                  aes(label = sample,x = V1,y = V2),
                  box.padding   = 1, 
                  point.padding = 0.5,
                  segment.color = 'grey50')+
  theme(legend.position = "bottom", 
        legend.text = element_text(size=20),
        legend.title=element_blank())+
  guides(shape = FALSE,fill=guide_legend(nrow=2, byrow=TRUE),
         color = guide_legend(nrow = 3, byrow = TRUE))+ 
  scale_color_discrete(name = "none")+
  scale_color_manual(values = custom_palette)+
  theme(text = element_text(size = 25))


p4

##
# plot
##

layout <- "
AAAABCDDDD
"

p_all = p + p2 + p3 + p4 +  plot_layout(design = layout)+
  plot_annotation(tag_levels = 'A')# , title = 'Figure 1'
p_all

