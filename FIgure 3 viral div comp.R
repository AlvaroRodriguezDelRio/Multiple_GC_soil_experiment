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


setwd("~/analysis/Berlin/soil 2019 metagenomes/viral_MAGs/drep/")

metadata = read.table("../../metadata.tab",header = T,sep = '\t')
metadata$Tube.number
head(metadata)
metadata$sample = paste(metadata$Tube.number,sep  = '')
head(metadata)

####
# load data
####

# load collapsed per tax data for plotting profiles
data = read.table("coverM/coverM.rel_abs.tab",header = F, sep = '\t')
names(data) = c("bin","species","sample","rel","rel_tot")
data$sample = tstrsplit(data$sample, "_")[[2]]
data$sample = str_replace_all(data$sample,"b",'')

data = data %>%
  left_join(metadata,by = "sample") %>% 
  filter(!is.na(species))

data[data$species=="NULL",]$species <- "Unclassified"
head(data)


###
# plot total vir per treat
###

d = data %>%
  group_by(sample) %>% 
  mutate(tot = sum(rel)) %>% 
  dplyr::select(sample,remark,tot) %>% 
  unique()

d$remark = factor(d$remark,levels = c("control","antibiotics", "copper", "drought",
                                      "microplastic", "Ndep", "salinity", "temp",
                                      "fungicide", "glyphosate","insecticide","Level 8"
))


ggplot(d)+
  geom_boxplot(aes(x = remark, y = tot))


#####
# plot viral abs per species
####

names(data)
d = data %>%
  group_by(sample) %>% 
  mutate(tot = sum(rel)) %>%
  group_by(species,sample) %>% 
  mutate(rel = sum(rel)) %>% 
  dplyr::select(sample,remark,tot,species,rel) %>% 
  unique()

d$remark = factor(d$remark,levels = c("control","antibiotics", "copper", "drought",
                                      "microplastic", "Ndep", "salinity", "temp",
                                      "fungicide", "glyphosate","insecticide","Level 8"
))

####
# plot profile for all species,
####



data$basal = paste(tstrsplit(data$species, ";")[[1]],tstrsplit(data$species, ";")[[2]],tstrsplit(data$species, ";")[[3]],sep = ";")
data[data$basal=="Unclassified;NA;NA",]$basal <- "Unclassified"


d = data %>%
  filter(sample != 85) %>% 
  group_by(sample) %>%
  mutate(tot = sum(rel)) %>%
  group_by(sample,basal) %>%
  mutate(rel_sp = sum(rel)) %>%
  ungroup() %>% 
  dplyr::select(tot,sample,remark,rel_sp,basal,Lv) %>% 
  unique()



d$remark = factor(d$remark,levels = c("control","antibiotics", "copper", "drought",
                                      "microplastic", "Ndep", "salinity", "temp",
                                      "fungicide", "glyphosate","insecticide","Level 8"
))


# plot profiled flipped 
p_prof = ggplot(d)+
  theme_void()+
  geom_histogram(aes(y = sample,x = rel_sp *100 / tot,fill = basal),stat = "identity")+
  theme(axis.text.x=element_text(angle=90, hjust=1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "bottom",
        strip.text.y.left = element_text(angle = 0),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8))+
  facet_wrap(remark~.,scales = 'free_y',ncol = 1,strip.position = "right")+
  scale_fill_brewer(palette="Spectral",name = "")+
  ylab("Relative abundance (%)")+
  xlab("")+
  scale_y_discrete(position = "left")+
  guides(fill=guide_legend(nrow=3, byrow=TRUE))


p_prof

# abs of unknown viruses 
dunk = d %>% 
  filter(sample != "85 salinity" & basal == 'Unclassified')
ggplot(dunk,
       aes(x = remark, y = rel_sp *100 / tot, fill = as.factor(Lv),color = as.factor(Lv)))+
  geom_boxplot(outlier.shape = NA,alpha = 0.8)+
  stat_compare_means(label = "p.signif", method = "wilcox.test",
                     ref.group = "control",hide.ns = TRUE, tip.length = 0)+
  coord_flip()+
  theme_minimal()+
  #  geom_hline(aes(yintercept= control_n), linetype="dashed", 
  #             color = "grey", size=1)+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2")+
  geom_jitter(shape=16, position=position_jitter(0.2),alpha = 0.8)+
  ylab("unknown virus proportion (%)")+
  xlab("Remark")+
  labs(fill = "Number of factors",color = "Number of factors")

# test unknown species 
wilcox.test(dunk[dunk$Lv == 8,]$rel_sp / dunk[dunk$Lv == 8,]$tot,dunk[dunk$Lv == 0,]$rel_sp)


#####
# Diversity 
####

head(data)


data_d = data %>% 
  filter(sample != 85) %>% 
  dplyr::select(bin ,rel,sample,remark) %>% 
  mutate(sample = paste(sample,remark)) %>% 
  group_by(bin,sample) %>% 
  mutate(rel_sp = sum(rel)) %>% 
  group_by(bin ) %>% 
  mutate(nsamples = sum(rel_sp > 0)) %>%
  #  filter(nsamples > 65) %>% 
  dplyr::select(-remark) %>%
  dplyr::select(-rel) %>%
  dplyr::select(-nsamples) %>% 
  unique() %>% 
  spread(sample,rel_sp)

data_d[is.na(data_d)] = 0
head(data_d)
r = (data_d$bin )
data_d$bin  = NULL
s = names(data_d)
t = data_d
head(t[,1:4])

div = microbiome::diversity(t,index = "shannon") # no jaccard 

div = data.frame(s,div)
div$treat = gsub("[0-9]* ","",div$s)
div = div %>% 
  mutate(lv = case_when(treat == "control" ~0,
                        treat == "Level8"~8,
                        .default = 1))



control = div %>% 
  filter(lv == 0)
control_median_rich = median(control$shannon)
div$control_n  = control_median_rich

div$treat = factor(div$treat,levels = c("control","antibiotics", "copper", "drought",
                                        "microplastic", "Ndep", "salinity", "temp",
                                        "fungicide", "glyphosate","insecticide","Level8"
))



levels(as.factor(d$species))
ggplot(div %>% 
         filter(s != "85 salinity"),
       aes(x = treat, y = shannon, fill = as.factor(lv),color = as.factor(lv)))+
  geom_boxplot(outlier.shape = NA,alpha = 0.8)+
  stat_compare_means(label = "p.signif", method = "wilcox.test",
                     ref.group = "control",hide.ns = TRUE, tip.length = 0)+
  coord_flip()+
  theme_minimal()+
  geom_hline(aes(yintercept= control_n), linetype="dashed", 
             color = "grey", size=1)+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2")+
  geom_jitter(shape=16, position=position_jitter(0.2),alpha = 0.8)+
  ylab("Viral diversity")+
  xlab("Remark")+
  labs(fill = "Number of factors",color = "Number of factors")



div$sample =  tstrsplit(div$s, " ")[[1]]
div = div %>% 
  left_join(metadata,by = "sample")

control = div %>% 
  filter(lv == 0)
control_median_rich = median(control$shannon)

df.summary <- div %>%
  group_by(remark,Lv) %>%
  summarise(
    sd = sd(shannon, na.rm = TRUE),
    len = median(shannon)
  )
df.summary


div$remark = factor(div$remark,levels = c("control","antibiotics", "copper", "drought",
                                          "microplastic", "Ndep", "salinity", "temp",
                                          "fungicide", "glyphosate","insecticide","Level 8"
))
div$remark <- factor(div$remark, levels = rev(levels(div$remark)))


df.summary$Lv = as.factor(df.summary$Lv)
p2 = ggplot(div)+
  geom_jitter(aes(y = shannon, x = (remark)),position = position_jitter(0.2), color = "darkgray",alpha = 0.2) + 
  geom_errorbar(aes(ymin = len-sd, ymax = len+sd,x = remark,color = Lv),data = df.summary,
                position = position_dodge(0.3), width = 0.2)+
  geom_point(aes(x = remark, y = len,color = Lv), data = df.summary,
             position = position_dodge(0.3)) +
  theme_minimal()+
  xlab("")+
  stat_compare_means(aes(group=remark, x = remark, y = shannon),label = "p.signif", method = "wilcox.test",
                     ref.group = "control",hide.ns = TRUE, tip.length = 0)+
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
  ylab("Sannon diversity")

p2 


####
# relevance of different factors on diversity 
###


facts = c("Drought..1.", "N.dep..2.", "Temp..3.", "Microplastic..4.", "Glyphosate..5.", "Antibiotics.oxytetra..6.", "Fungicide..7.","Cu..8.", "Salinity..9.", "Insect.imidacloprid..10.", "Lv")

c_div = div %>% 
  dplyr::select(shannon,facts) %>% 
  mutate(Lv = case_when(Lv == 8~1,
                        .default = 0))

regressor <- randomForest(shannon ~ ., data = c_div, importance=TRUE) # fit the random forest with default parameter
relImp = as.data.frame(importance(regressor))# %IncMSE is the most robust and informative measure. It is the increase in mse of predictions
#(estimated with out-of-bag-CV) as a result of variable j being permuted(values randomly shuffled).
# the higher number, the more important
names(relImp)
relImp$treat = rownames(relImp)

ggplot(relImp)+
  geom_bar(aes(x = `%IncMSE`, y = treat),stat = "identity")

relImp

#####
# PCOA per bin
####
head(data)
data_d = data %>% 
  filter(sample != 85) %>% 
  dplyr::select(bin,rel,sample,remark) %>% 
  mutate(sample = paste(sample,remark)) %>% 
  group_by(bin,sample) %>% 
  mutate(rel_sp = sum(rel)) %>% 
  group_by(bin) %>% 
  mutate(nsamples = sum(rel_sp > 0)) %>%
  #  filter(nsamples > 65) %>% 
  dplyr::select(-remark) %>%
  dplyr::select(-rel) %>%
  dplyr::select(-nsamples) %>% 
  unique() %>% 
  spread(sample,rel_sp)

dpcoa = data_d
rownames_ = as.character(dpcoa$bin)
dpcoa$bin = NULL
s_names = names(dpcoa)
dpcoa = as.matrix(dpcoa)
dpcoa = t(dpcoa)

rownames(dpcoa) = s_names
sample_order = gsub("[0-9]* ","",rownames(dpcoa)) 
dpcoa[is.na(dpcoa)] = 0

d_m = vegdist(dpcoa, method="bray")
pcoa = ape::pcoa(d_m)
pcoa_plot = as.data.frame(pcoa$vectors[,1:2])

pcoa_plot$remark = sample_order
pcoa_plot = pcoa_plot %>% mutate(lv = case_when(remark == "control" ~0,
                                                remark == "Level8"~8,
                                                .default = 1))

names(pcoa_plot) = c("V1","V2","remark","lv")
pcoa_plot$lv = as.factor(pcoa_plot$lv)
pcoa_plot$s = rownames(pcoa_plot)
pcoa_plot$sample = tstrsplit(pcoa_plot$s, " ")[[1]]
pcoa_plot$`Number of factors` = pcoa_plot$lv

head(pcoa_plot)
pcoa_plot$remark = factor(pcoa_plot$remark,levels = c("control","antibiotics", "copper", "drought",
                                            "microplastic", "Ndep", "salinity", "temp",
                                            "fungicide", "glyphosate","insecticide","Level8"
))

custom_palette <- c("#999999", "#377EB8", "#F781BF","#4DAF4A", "#FF7F00",
                    "#E41A1C", "#A65628", "#984EA3", "#FFFF33", "#66C2A5", "black", "#8EBA42")



p_cor = ggplot(pcoa_plot)+
  geom_point(aes(x = V1, y = V2,color = remark,shape = `Number of factors`),size = 3)+
  theme_classic()+
  xlab(paste("PCoA 1 (",round(pcoa$values$Relative_eig[1]*100,digits = 2),"%)",sep = ""))+
  ylab(paste("PCoA 2 (",round(pcoa$values$Relative_eig[2]*100,digits = 2),"%)",sep = ""))+
  #  theme(text = element_text(size = 20))+
  geom_text_repel(data = pcoa_plot[pcoa_plot$s %in% c("124 Level 8","127 Level 8","128 Level 8"),],
                  aes(label = sample,x = V1,y = V2),
                  box.padding   = 3, 
                  point.padding = 0.5,
                  segment.color = 'grey50')+
  theme(legend.position = "bottom", 
        legend.text = element_text(size=9),
        legend.title = element_text(size=9))+
  guides(shape = FALSE)+
  scale_color_manual(values = custom_palette)

p_cor

#####
# plot combined
####

layout <- "
AAABCCC
"
p_prof + p2 + p_cor +  plot_layout(design = layout)+plot_annotation(tag_levels = 'A', title = 'Figure 3')

pdf("~/analysis/Berlin/soil 2019 metagenomes/Figures/Figure 3 viruses.pdf", width=20, height=8)
p_prof + p2 + p_cor +  plot_layout(design = layout)+plot_annotation(tag_levels = 'A',, title = 'Figure 3')
graphics.off()
