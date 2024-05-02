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
metadata$Tube.number
metadata$sample = paste(metadata$Tube.number,sep  = '')

####
# Composition plot 
####

# load collapsed per tax data for plotting profiles
data = read.table("n_genes_per_sample.fams_no_emapper.fams_more_50_samples.tab",header = F, sep = '\t')
names(data) = c("species","sample","abs")
data$sample = as.character(data$sample)

# quit samples with no data (not in the 70 samples sequenced here)
samples = data %>% 
  filter(abs > 0) %>% 
  dplyr::select(sample) %>% 
  unique()

data =  data %>% 
  filter(sample %in% samples$sample)

# joim with metadata
data = data %>%
  left_join(metadata,by = "sample")

# load marker copy number per samplpe
mabs = read.table("../../maker_gene_count/maker_count_per_sample.tab", header = F, sep = "\t")
names(mabs) = c("sample","mean_n_mar","median_n_mar","total_n_mar")
mabs$sample = tstrsplit(mabs$sample, "_")[[2]]
mabs$sample = str_replace_all(mabs$sample,"b",'')
mabs = mabs %>% 
  group_by(sample) %>% 
  dplyr::select(sample,mean_n_mar) %>% 
  unique()
mabs$sample = as.character(mabs$sample)
data = data %>% 
  left_join(mabs,by = "sample")
data$cnpc = data$abs / data$mean_n_mar

data_d = data %>% 
  filter(sample != "85") %>% 
  dplyr::select(species,cnpc,sample,remark) %>% 
  mutate(sample = paste(sample,remark)) %>% 
  dplyr::select(-remark) %>% 
  spread(sample,cnpc)

dpcoa = data_d
rownames_ = as.character(dpcoa$species)
dpcoa$species = NULL
s_names = names(dpcoa)
dpcoa = as.matrix(dpcoa)
dpcoa = t(dpcoa)

rownames(dpcoa) = s_names
sample_order = gsub("[0-9]* ","",rownames(dpcoa)) 

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


data_pcoa = pcoa_plot
pcoa_plot = ggplot(pcoa_plot)+
  geom_point(aes(x = V1, y = V2,color = remark,shape = `Number of factors`),size = 5)+
  theme_classic()+
  xlab(paste("PCoA 1 (",round(pcoa$values$Relative_eig[1]*100,digits = 2),"%)",sep = ""))+
  ylab(paste("PCoA 2 (",round(pcoa$values$Relative_eig[2]*100,digits = 2),"%)",sep = ""))+
 # theme(text = element_text(size = 20))+
  geom_text_repel(data = pcoa_plot[pcoa_plot$s %in% c("124 Level 8","127 Level 8","128 Level 8"),],
                  aes(label = sample,x = V1,y = V2),
                  box.padding   = 3, 
                  point.padding = 0.5,
                  segment.color = 'grey50')+
  scale_color_manual(values = custom_palette)


#####
# n novel fams 
####

d = read.table("clusters.sizes.tab")
names(d) = c("fam","size")
n = read.table("clusters.folded.no_emapper.reps.txt")
names(n) = c("fam")

d = d %>% 
  mutate(novel = case_when(fam %in% n$fam ~ "Novel",
                           .default =  "Not novel"))

p1 = ggplot(d)+
  geom_histogram(aes(x = log10(size), fill = novel),bins = 10)+
  theme_classic()+
  scale_fill_manual(values = c("red","grey"))+
  xlab("Log10 (# genes per family)")+
  theme(legend.position = "none")


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
  scale_x_continuous(breaks = c(0,10000000,20000000))

p2
options(scipen = 999)
p_num_novel = p1 + inset_element(p2,  0.4, 0.6, 1, 0.9)
p_num_novel


########
# combine 
#######

top_plot      <- wrap_elements((p_num_novel) + 
                                 plot_annotation(title = "A"))
bottom_plot   <- wrap_elements(pcoa_plot + 
                                 plot_annotation(title = "B"))
top_plot + bottom_plot 

pdf("~/analysis/Berlin/soil 2019 metagenomes/Figures/Figure 6 novel families.pdf",height = 7,width = 15)
top_plot + bottom_plot + plot_annotation(title = 'Figure 6')
graphics.off()
