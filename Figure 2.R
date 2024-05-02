
library(tidytree)
library(ggtree)
library(dplyr)
library(data.table)
library(stringr)
library(dplyr)
library(tidyverse)
library(ggridges)
library(ggpubr)
library(patchwork)


tree <- read.tree("alignment.tr.oneline.rep.faa.treefile.txt",tree.names = T)
a = as_tibble(tree) 

# format tree tip names 
tree$tip.label = tstrsplit(tree$tip.label, "\\.")[[1]]
tree$tip.label = substr(tree$tip.label, 1, nchar(tree$tip.label) - 1)
tips = tree$tip.label

# add genome annotations to tree
tax = read.table("../../genome2taxonomy.tab",sep = '\t')
names(tax) =  c("genome","tax")
s = tibble(label = tax$genome,tax = tax$tax)
tree = full_join(tree, s, by = 'label')


# add myco - non-myco  to tree
tax = tax %>% 
  mutate(path = case_when(grepl("Mycobacterium",tax)~"Mycobacterium",
                   .default = ""))

s = tibble(label = tax$genome,path = tax$path)
tree = full_join(tree, s, by = 'label')

# get VF count per genome
abs = read.table("VFDB_per_bin.curated.tab",sep ='\t',header = F)
names(abs) = c("tip","tax","gene","gene_v","desc")
head(abs)
abs$tax = NULL
abs$gene = NULL
abs$desc = NULL
abs$gene_v = tstrsplit(abs$gene_v, "\\(")[[1]]

abs_t = abs %>%
  group_by(tip,gene_v) %>% 
  mutate(n = n()) %>% 
  unique() %>% 
  spread(key = gene_v,value = n)

abs_t = as.data.frame(abs_t)
rownames(abs_t) = abs_t$tip
abs_t$tip = NULL
head(abs_t)
abs_t[is.na(abs_t)] = 0
dim(abs_t)
tips
abs_t

# add missing genomes 
for (row_name in tips) {
  if (!(row_name %in% rownames(abs_t))) {
    d = data.frame(matrix(0, nrow = 1, ncol = ncol(abs_t)))
    print(row_name)
    rownames(d) = row_name
    names(d) = names(abs_t)
    abs_t <- rbind(abs_t,d)
  }
}


# plot tree
custom_palette <- c("white", "lightblue", "orange")
p = ggtree(tree)+# layout = "equal_angle",branch.length = 'none'
  geom_tippoint(aes(color=path), size=3, alpha=1)+ # our MAGs / no our MAGs
  scale_color_manual(name = "",values=custom_palette) 
  

# add heatmap
p1 = gheatmap(p, abs_t,
         colnames_angle=90, 
         hjust=2, font.size=0,high = 'blue',low = 'white',offset=1,legend_title = "Number",width = 10)


#####
# add mycobacterium abundance plot 
#####

metadata = read.table("metadata.tab",header = T,sep = '\t')
metadata$Tube.number
metadata$sample = paste(metadata$Tube.number,sep  = '')

data = read.table("../coverM/coverM.rel_abs.per_tax.tab",header = F, sep = '\t')
names(data) = c("species","sample","rel")
data$sample = tstrsplit(data$sample, "_")[[2]]
data$sample = str_replace_all(data$sample,"b",'')

data = data %>%
  left_join(metadata,by = "sample")


# plot myco relative abundance per treatmen
data_ = data %>% 
  filter(grepl("Mycobacterium",species)) %>% 
  filter(!grepl("s__",species))

mm = data_ %>% 
  filter(remark == 'control') %>%
  group_by(remark) %>% 
  mutate(m = median(rel)) %>% 
  dplyr::select(m) %>% 
  unique()

p_myco = ggplot(data_,
       aes(x = Lv, y = rel, fill = as.factor(Lv)))+
  geom_boxplot(outlier.shape = NA,alpha = 0.8)+
  stat_compare_means(aes(group = as.factor(Lv)),label = "p.signif", method = "wilcox.test",
                     ref.group = "0",hide.ns = TRUE, tip.length = 0)+
  theme_minimal()+
  geom_hline(yintercept= mm$m, linetype="dashed", 
             color = "grey", size=1)+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2")+
  geom_jitter(aes(color = as.factor(Lv)),shape=16, position=position_jitter(0.2),alpha = 0.8)+
  ylab("Relative abudance (%)")+
  xlab("Number of factors")+
  labs(fill = "Number of factors",color = "Number of factors")+
  theme(legend.position = "none")  # Remove legend title

p_myco

#####
# add plot motus mycobacterium species abundance
#####

setwd("~/analysis/Berlin/soil 2019 metagenomes/tax_profile/prok/")
data = read.table("motus_merged.rarefied.tsv",header = T, sep = '\t') 

data_g = data %>% 
  gather(sample,"n_reads",-mOTUs2_clade) %>% 
  left_join(metadata,by = "sample")

data_g$species = tstrsplit(data_g$mOTUs2_clade, "\\|")[[7]]
data_g$species = str_replace_all(data_g$species,"s__","")

d_bacteria = data_g %>% 
  filter(mOTUs2_clade == 'k__Bacteria') %>% 
  dplyr::select(sample,n_reads)

names(d_bacteria) = c("sample","tot_bac_reads")

d = data_g %>% 
  filter(grepl("Mycobacterium",mOTUs2_clade) & grepl("s__",mOTUs2_clade))

# group by species
d = d %>% 
  left_join(d_bacteria,by = "sample") %>% 
  group_by(sample,species) %>% 
  mutate(rel = sum(n_reads)) %>% 
  dplyr::select(species,sample,rel,remark,Lv) %>% 
  unique()


# filter for significant
ds <- d %>%
  filter(remark %in% c("control","Level 8")) 
d$rel = d$rel *100 / 31310395

results = ds %>% 
  group_by(species) %>%
  summarise(p_value = wilcox.test(rel ~ remark)$p.value)

significant_factors <- results %>%
  filter(p_value < 0.05) %>%
  filter(!species %in% c("Mycobacterium sp. 852002-50816_SCH5313054-b [ref_mOTU_v31_04182]"))

d = d %>% 
  filter(species %in% significant_factors$species)

# plot significant
mm = d %>% 
  filter(remark == 'control') %>%
  group_by(remark,species) %>% 
  mutate(m = median(rel)) %>%
  ungroup() %>% 
  dplyr::select(m,species) %>% 
  unique()

d = d %>% 
  left_join(mm,by = "species") %>% 
  filter(sample != "X85")

p_motus = ggplot(d,aes(x = remark, y = rel))+
  geom_boxplot(outlier.shape = NA,alpha = 0.8,aes(fill = as.factor(Lv)))+
  stat_compare_means(aes(group = remark),label = "p.signif", method = "wilcox.test",
                     ref.group = "control",hide.ns = TRUE, tip.length = 0)+
  facet_wrap(~species,nrow = 1)+
  coord_flip()+
  theme_minimal() +
  geom_hline(aes(yintercept = m), linetype="dashed", 
             color = "grey", size=1)+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2")+
  geom_jitter(aes(color = as.factor(Lv)),shape=16, position=position_jitter(0.2),alpha = 0.8)+
  ylab("Percentage of reads (%)")+
  xlab("")+
  labs(fill = "Number of factors",color = "Number of factors")+
  theme(legend.position = "none",
        strip.text.x = element_text(
          size = 7
        ))  # Remove legend title

p_motus

#####
# combine
#####


layout <- "
ABBBBB
ABBBBB
CCCCCC
"

p_all = p_myco + p1 + p_motus + 
  plot_layout(design = layout)+
  plot_annotation(tag_levels = 'A', title = 'Figure 2')

p_all
  
