library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)
library(vegan)
library(ecodist)
library(MASS)
library(data.table)
library(stringr)
library(ggsignif)
library(microbiome)
library(ggpubr)
library(relaimpo)
library(randomForest)
library(ggrepel)

metadata = read.table("../../metadata.tab",header = T,sep = '\t')
metadata$sample = paste(metadata$Tube.number,sep  = '')

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
data$basal = paste(tstrsplit(data$species, ";")[[1]],tstrsplit(data$species, ";")[[2]],tstrsplit(data$species, ";")[[3]],sep = ";")
data[data$basal=="Unclassified;NA;NA",]$basal <- "Unclassified"

#####
# plot profile
####

d = data %>%
  filter(sample != 85) %>% 
  group_by(sample) %>%
  mutate(tot = sum(rel)) %>%
  group_by(sample,basal) %>%
  mutate(rel_sp = sum(rel)) %>%
  ungroup() %>% 
  dplyr::select(tot,sample,remark,rel_sp,basal,Lv) %>% 
  unique()

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



#####
# add diversity plot 
####

data_d = data %>% 
  filter(sample != 85) %>% 
  dplyr::select(bin ,rel,sample,remark) %>% 
  mutate(sample = paste(sample,remark)) %>% 
  group_by(bin,sample) %>% 
  mutate(rel_sp = sum(rel)) %>% 
  group_by(bin ) %>% 
  mutate(nsamples = sum(rel_sp > 0)) %>%
  dplyr::select(-remark) %>%
  dplyr::select(-rel) %>%
  dplyr::select(-nsamples) %>% 
  unique() %>% 
  spread(sample,rel_sp)

data_d[is.na(data_d)] = 0
r = (data_d$bin )
data_d$bin  = NULL
s = names(data_d)
t = data_d

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

#####
# Add composition plot bin
####

data_d = data %>% 
  filter(sample != 85) %>% 
  dplyr::select(bin,rel,sample,remark) %>% 
  mutate(sample = paste(sample,remark)) %>% 
  group_by(bin,sample) %>% 
  mutate(rel_sp = sum(rel)) %>% 
  group_by(bin) %>% 
  mutate(nsamples = sum(rel_sp > 0)) %>%
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


#####
# plot combined
####

layout <- "
AAABCCC
"
p_prof + p2 + p_cor +  plot_layout(design = layout)+plot_annotation(tag_levels = 'A', title = 'Figure 3')

pdf("Figure 3.pdf", width=20, height=8)
p_prof + p2 + p_cor +  plot_layout(design = layout)+plot_annotation(tag_levels = 'A',, title = 'Figure 3')
graphics.off()
