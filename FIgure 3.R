#install.packages("tidytree")
library(dplyr)
library(data.table)
library(stringr)
library(dplyr)
library(tidyverse)
library(ggridges)
library(ggpubr)
library(patchwork)
library(scales)


setwd("~/analysis/Berlin/soil 2019 metagenomes/Figures/Source_data/")

#3A

data = read.csv('3A.csv')

mm = data %>% 
  filter(remark == 'control') %>%
  group_by(remark) %>% 
  mutate(m = median(rel)) %>% 
  dplyr::select(m) %>% 
  unique()

data$Lv = as.character(data$Lv)
p_myco = ggplot(data,
                aes(x = Lv, y = rel, fill = as.factor(Lv)))+
  #  ggridges::geom_density_ridges()
  geom_boxplot(outlier.shape = NA,alpha = 0.8)+
  stat_compare_means(aes(group = as.factor(Lv)),label = "p.signif", method = "wilcox.test",
                     ref.group = "0",hide.ns = TRUE, tip.length = 0)+
  #coord_flip()+
  theme_minimal()+
  geom_hline(yintercept= mm$m, linetype="dashed", 
             color = "grey", size=1)+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2")+
  geom_jitter(aes(color = as.factor(Lv)),shape=16, position=position_jitter(0.2),alpha = 0.8)+
  ylab("Relative abudance (%)")+
  xlab("Number of factors")+
  labs(fill = "Number of factors",color = "Number of factors")+
  theme(legend.position = "none")+
  coord_flip()+
  theme(text = element_text(size = 19))


p_myco

# 3B

d = read.csv('3B.csv')

mm = d %>% 
  filter(remark == 'control') %>%
  group_by(remark,species) %>% 
  mutate(m = median(rel)) %>%
  ungroup() %>% 
  dplyr::select(m,species) %>% 
  unique()




d$remark = factor(d$remark,levels = c("Control","Warming","Drought","Nitrogen deposition",
                                      "Salinity","Heavy metals","Microplastics",
                                      "Antibiotics","Fungicide", "Herbicide",
                                      "Insecticide","8 factors"
))

p_motus = ggplot(d,aes(x = remark, y = rel))+
  geom_boxplot(outlier.shape = NA,alpha = 0.8,aes(fill = as.factor(Lv)))+
  stat_compare_means(aes(group = remark),label = "p.signif", method = "wilcox.test",
                     ref.group = "Control",hide.ns = TRUE, tip.length = 0)+
  facet_wrap(~species,nrow = 1,scales = "free_x")+
  coord_flip()+
  theme_minimal() +
  geom_hline(aes(yintercept = m), linetype="dashed", 
             color = "grey", size=1)+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2")+
  geom_jitter(aes(color = as.factor(Lv)),shape=16, position=position_jitter(0.2),alpha = 0.8)+
  ylab("Relative abundance (%)")+
  xlab("")+
  theme(text = element_text(size = 20))+
  labs(fill = "Number of factors",color = "Number of factors")+
  theme(legend.position = "none",
        strip.text.x = element_text(
          size = 10
        ))+  # Remove legend title
  theme(
    strip.text.x = element_text(size = 16), # For column facet titles
    strip.text.y = element_text(size = 15))+
  theme(axis.text.y = element_text(size=20))+
  scale_y_continuous(breaks = pretty_breaks(n = 2))



p_motus

#####
# combine
#####


layout <- "
AAAAA
BBBBB
BBBBB
"

p_all = p_myco + p_motus + 
  plot_layout(design = layout)+
  plot_annotation(tag_levels = 'A')# , title = 'Figure 2'


p_all




