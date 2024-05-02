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


setwd("~/analysis/Berlin/soil 2019 metagenomes/gene_prediction/novel_families/")

metadata = read.table("../../metadata.tab",header = T,sep = '\t')
metadata$Tube.number
head(metadata)
metadata$sample = paste(metadata$Tube.number,sep  = '')



####
# pcoa
####

# load collapsed per tax data for plotting profiles
data = read.table("n_genes_per_sample.fams_no_emapper.fams_more_50_samples.tab",header = F, sep = '\t')
head(data)
names(data) = c("species","sample","abs")
head(data)
data$sample = as.character(data$sample)

# quit samples with no data (not in the 70 samples used)
samples = data %>% 
  filter(abs > 0) %>% 
  dplyr::select(sample) %>% 
  unique()


nrow(data)
data =  data %>% 
  filter(sample %in% samples$sample)
nrow(data)

# joim with metadata
head(metadata)
data = data %>%
  left_join(metadata,by = "sample")

head(data)

# load marker copy number per samplpe
mabs = read.table("../../maker_gene_count/maker_count_per_sample.tab", header = F, sep = "\t")
head(mabs)
names(mabs) = c("sample","mean_n_mar","median_n_mar","total_n_mar")
mabs$sample = tstrsplit(mabs$sample, "_")[[2]]
mabs$sample = str_replace_all(mabs$sample,"b",'')
head(mabs)
mabs = mabs %>% 
  group_by(sample) %>% 
  dplyr::select(sample,mean_n_mar) %>% 
  unique()
mabs$sample = as.character(mabs$sample)
data = data %>% 
  left_join(mabs,by = "sample")
nrow(data)
head(data)
data$cnpc = data$abs / data$mean_n_mar


#####
# PCOA
####

head(data)

data_d = data %>% 
  filter(sample != "85") %>% 
  dplyr::select(species,cnpc,sample,remark) %>% 
  mutate(sample = paste(sample,remark)) %>% 
  dplyr::select(-remark) %>% 
  spread(sample,cnpc)

dpcoa = data_d
head(dpcoa)

rownames_ = as.character(dpcoa$species)
dpcoa$species = NULL
s_names = names(dpcoa)
dpcoa = as.matrix(dpcoa)
dpcoa = t(dpcoa)

dim(dpcoa)
length(rownames_)
head(dpcoa)

rownames(dpcoa) = s_names
sample_order = gsub("[0-9]* ","",rownames(dpcoa)) 
head(dpcoa[,1:4])
dim(dpcoa)
dpcoa[is.na(dpcoa)] = 0

#d_m = dist(dpcoa)
d_m = vegdist(dpcoa, method="bray")
pcoa = ape::pcoa(d_m)
pcoa_plot = as.data.frame(pcoa$vectors[,1:2])
nrow(pcoa_plot)

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


pcoa_plot

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
head(d)

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

######
# PcoA axis 1 in different treatments 
#####

pcoa_mags = read.table("~/analysis/Berlin/soil 2019 metagenomes/MAGs/semibin/dRep/coverM/composition.no85.tab",
                       header = T)
pcoa_mags = pcoa_mags %>% 
  dplyr::select(sample,V1) %>% 
  mutate(origin = "MAGs")

pcoa_kos = read.table("~/analysis/Berlin/soil 2019 metagenomes/gene_prediction/KOs/kos_PCOA_copy_per_cell.plot.tab",
                      header = T)
pcoa_kos = pcoa_kos %>% 
  dplyr::select(sample,Axis.1) %>% 
  mutate(origin = "KOs")
head(pcoa_kos)
names(pcoa_kos) = c("sample","V1","origin")

data_pcoa_novel = data_pcoa %>% 
  dplyr::select(sample,V1) %>% 
  mutate(origin = "Novel families")

d_all = rbind(pcoa_mags,pcoa_kos,data_pcoa_novel)
d_all = d_all %>% 
  left_join(metadata,by = "sample")


d_all$Lv = as.factor(d_all$Lv)
ggplot(d_all,aes(x = as.factor(Lv),y = V1))+
  geom_boxplot()+
  coord_flip()+
  facet_wrap(~origin)+
  stat_compare_means(aes(group=Lv, x = Lv, y = V1),label = "p.signif", method = "t.test",
                     ref.group = "0",hide.ns = TRUE, tip.length = 0)


######
# corr lv unknonw / n novel fams
#####

f = read.table("num_nfams_per_bin.unk_lev.tab",sep = '\t')
head(f)
names(f) = c("b","t","unk","n","tot_genes")

ff = f %>% 
  mutate(prop = n / tot_genes) %>% 
  group_by(unk) %>% 
  summarise(p = mean(prop))

ggplot(ff)+
  geom_bar(aes(x = p, y = unk),stat = "identity")


f %>% 
  mutate(prop = n / tot_genes) %>% 
  filter(prop > 0.2)

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
