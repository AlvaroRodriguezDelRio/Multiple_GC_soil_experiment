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
library(ggalt)

metadata = read.table("metadata.tab",header = T,sep = '\t')
metadata$Tube.number
head(metadata)
metadata$sample = paste(metadata$Tube.number,sep  = '')


####
# load data
####

# load collapsed per tax data for plotting profiles
data = read.table("coverM/coverM.rel_abs.per_tax.tab",header = F, sep = '\t')
head(data)
names(data) = c("species","sample","rel")
data$sample = tstrsplit(data$sample, "_")[[2]]
data$sample = str_replace_all(data$sample,"b",'')

data = data %>%
  left_join(metadata,by = "sample")

####
# plot profile for all species, no leggends, to get an idea of composition, combine with diversity and genome size 
####

d = data %>%
  filter(sample != 85) %>% 
  filter(grepl("c__",species)) %>%
  filter(!grepl("o__",species)) %>% 
  group_by(sample) %>% 
  mutate(tot = sum(rel))

d$remark = factor(d$remark,levels = c("control","antibiotics", "copper", "drought",
                                      "microplastic", "Ndep", "salinity", "temp",
                                      "fungicide", "glyphosate","insecticide","Level 8"
))


# plot profile flipped 
p = ggplot(d)+
  theme_void()+
  geom_histogram(aes(y = sample,x = rel *100 / tot,fill = species),stat = "identity")+
  theme(axis.text.x=element_text(angle=90, hjust=1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "bottom",
        strip.text.y.left = element_text(angle = 0),
        legend.text = element_text(size=6),
        legend.title = element_text(size=6))+
  facet_wrap(remark~.,scales = 'free_y',ncol = 1,strip.position = "right")+
  scale_fill_brewer(palette="Spectral",name = "")+
  ylab("Relative abundance (%)")+
  xlab("")+
  scale_y_discrete(position = "left")+
  guides(fill=guide_legend(nrow=5, byrow=TRUE))

##
# add diversity plot
##

div = read.table("coverM/diversity.tab",sep = '\t')
div$sample = tstrsplit(div$s, " ")[[1]]
div = div %>% 
  filter(sample!= 85) %>% 
  dplyr::select(sample,shannon) %>% 
  left_join(metadata,by = 'sample')


control = div %>% 
  filter(Lv == 0)
control_median_rich = median(control$shannon)

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
                     ref.group = "control",hide.ns = TRUE, tip.length = 0,paired = F)+
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


##
# add genome size plot 
##

# load collapsed per tax data for plotting profiles
data = read.table("../genome2nbps.MQ.tab",header = F, sep = '\t')
names(data) = c("bin","n")
data$sample = tstrsplit(data$bin, "_")[[2]]
data$sample = str_replace_all(data$sample,"b",'')
data = data %>%
  left_join(metadata,by = "sample")

data %>% 
  filter(remark == 'control') %>%
  group_by(remark) %>% 
  mutate(m = median(n)) %>% 
  dplyr::select(m)

data$remark <- factor(data$remark, levels = rev(levels(data$remark)))
data$n = data$n / 1000000

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


p3 = ggplot(data)+
  geom_jitter(aes(y = n, x = remark),position = position_jitter(0.2), color = "darkgray",alpha = 0.2) + 
  geom_errorbar(aes(ymin = len-sd, ymax = len+sd,x = remark,color = Lv),data = df.summary,
                position = position_dodge(0.3), width = 0.2,)+
  geom_point(aes(x = remark, y = len,color = Lv), data = df.summary,
             position = position_dodge(0.3)) +
  theme_minimal()+
  xlab("")+
  stat_compare_means(aes(group=remark, x = remark, y = n),label = "p.signif", method = "wilcox.test",paired = F,
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
  ylab("Genome size (Mbps)")


p3


##
# add composition plot
##

data = read.table("coverM/coverM.rel_abs.tab",header = F, sep = '\t')
names(data) = c("genome","species","sample","rel","rel_tot")
data$sample = tstrsplit(data$sample, "_")[[2]]
data$sample = str_replace_all(data$sample,"b",'')
data = data %>%
  left_join(metadata,by = "sample")

data_d = data %>% 
  filter(sample != 85) %>% 
  filter(!grepl("UNCLASSIFIED",species)) %>% 
  dplyr::select(genome,rel,sample,remark) %>% 
  mutate(sample = paste(sample,remark)) %>% 
  dplyr::select(-remark) %>% 
  spread(sample,rel)

dpcoa = data_d
rownames_ = as.character(dpcoa$genome)
dpcoa$genome = NULL
s_names = names(dpcoa)
dpcoa = as.matrix(dpcoa)
dpcoa = t(dpcoa)

rownames(dpcoa) = s_names
sample_order = gsub("[0-9]* ","",rownames(dpcoa)) 
head(dpcoa[,1:4])
dpcoa[is.na(dpcoa)] = 0

d_m = vegdist(dpcoa, method="bray")
pcoa = ape::pcoa(d_m)
pcoa_plot = as.data.frame(pcoa$vectors[,1:2])


pcoa_plot$remark = sample_order
pcoa_plot = pcoa_plot %>% mutate(lv = case_when(remark == "control" ~0,
                                                remark == "Level8"~8,
                                                .default = 1))

names(pcoa_plot) = c("V1","V2","remark","lv")
pcoa_plot$sample = tstrsplit(rownames(pcoa_plot), " ")[[1]]
pcoa_plot$`Number of factors` = as.factor(pcoa_plot$lv)
pcoa_plot$remark = as.factor(pcoa_plot$remark)


pcoa_plot$remark = factor(pcoa_plot$remark,levels = c("control","antibiotics", "copper", "drought",
                                      "microplastic", "Ndep", "salinity", "temp",
                                      "fungicide", "glyphosate","insecticide","Level8"
))

pcoa_plot$`Number of factors` = as.factor(pcoa_plot$lv)

custom_palette <- c("#999999", "#377EB8", "#F781BF","#4DAF4A", "#FF7F00",
                    "#E41A1C", "#A65628", "#984EA3", "#FFFF33", "#66C2A5", "black", "#8EBA42")

p4 = ggplot(pcoa_plot)+
  geom_point(aes(x = V1, y = V2,color = remark,shape = `Number of factors`),size = 4)+
  theme_classic()+
  xlab(paste("PCoA 1 (",round(pcoa$values$Relative_eig[1]*100,digits = 2),"%)",sep = ""))+
  ylab(paste("PCoA 2 (",round(pcoa$values$Relative_eig[2]*100,digits = 2),"%)",sep = ""))+
  #theme(text = element_text(size = 20))+
  geom_text_repel(data = pcoa_plot[pcoa_plot$sample %in% c(124,127,128),],
                  aes(label = sample,x = V1,y = V2),
                  box.padding   = 1, 
                  point.padding = 0.5,
                  segment.color = 'grey50')+
  theme(legend.position = "bottom", 
        legend.text = element_text(size=9))+
  guides(shape = FALSE)+
  scale_color_discrete(name = "")+
  scale_color_manual(values = custom_palette)



p4



##
# plot
##

layout <- "
AAAABCDDDD
"

p_all = p + p2 + p3 + p4 +  plot_layout(design = layout)+plot_annotation(tag_levels = 'A', title = 'Figure 1')
p_all

pdf("~/analysis/Berlin/soil 2019 metagenomes/Figures/Figure 1.pdf", width=20, height=8)
print(p_all)
graphics.off()

