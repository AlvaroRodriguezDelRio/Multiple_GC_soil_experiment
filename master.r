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
library(coin)

setwd("~/analysis/Berlin/soil 2019 metagenomes/MAGs/semibin/dRep/coverM/")

# load metadata
metadata = read.table("../../../../metadata.tab",header = T,sep = '\t')
metadata$Tube.number
head(metadata)
metadata$sample = paste(metadata$Tube.number,sep  = '')
metadata %>% 
  filter(remark == 'control') %>% 
  dplyr::select(sample)

######
# Bacteria
######

# diff rel abs of particular lineages
data = read.table("coverM.rel_abs.per_tax.tab",header = F, sep = '\t')
names(data) = c("species","sample","rel")
data$sample = tstrsplit(data$sample, "_")[[2]]
data$sample = str_replace_all(data$sample,"b",'')
data = data %>%
  left_join(metadata,by = "sample") %>% 
  filter(sample != 85)

### c__Actinomycetia, Level 8
d = data %>%
  filter(grepl("c__",species)) %>%
  filter(!grepl("o__",species)) %>% 
  group_by(sample) %>% 
  mutate(tot = sum(rel))

wilcox.test(d[d$species == 'd__Bacteria;p__Actinobacteriota;c__Actinomycetia' & d$remark == 'control',]$rel,
            d[d$species == 'd__Bacteria;p__Actinobacteriota;c__Actinomycetia' & d$remark == 'Level 8',]$rel)

wilcox.test(d[d$species == 'd__Bacteria;p__Actinobacteriota;c__Actinomycetia' & d$remark == 'control',]$rel,
            d[d$species == 'd__Bacteria;p__Actinobacteriota;c__Actinomycetia' & d$remark == 'salinity',]$rel)

d = d %>% 
  filter(species == 'd__Bacteria;p__Actinobacteriota;c__Actinomycetia')

ggplot(d)+
  geom_boxplot(aes(x = remark, y = rel))


#### p__Actinobacteriota, level 8
d = data %>%
  filter(grepl("p__",species)) %>%
  filter(!grepl("c__",species)) %>% 
  group_by(sample) %>% 
  mutate(tot = sum(rel))

wilcox.test(d[d$species == 'd__Bacteria;p__Actinobacteriota' & d$remark == 'control',]$rel,
            d[d$species == 'd__Bacteria;p__Actinobacteriota' & d$remark == 'Level 8',]$rel)

d = d %>% 
  filter(species == 'd__Bacteria;p__Actinobacteriota')

ggplot(d)+
  geom_boxplot(aes(x = remark, y = rel))


# difference in 8 factor samples non salinity treatment vs salinity treatment
wilcox.test(d[d$species == 'd__Bacteria;p__Actinobacteriota' & d$remark == 'control',]$rel,
            d[d$species == 'd__Bacteria;p__Actinobacteriota' & d$remark == 'Level 8' & d$sample %in% c(124, 127, 128),]$rel)

wilcox.test(d[d$species == 'd__Bacteria;p__Actinobacteriota' & d$remark == 'control',]$rel,
            d[d$species == 'd__Bacteria;p__Actinobacteriota' & d$remark == 'Level 8' & !d$sample %in% c(124, 127, 128),]$rel)


ggplot(d %>% 
         filter(sample %in% c(124,127,128) | remark != "Level 8"))+
  geom_boxplot(aes(x = remark, y = rel))

ggplot(d %>% 
         filter(!sample %in% c(124,127,128) | remark != "Level 8"))+
  geom_boxplot(aes(x = remark, y = rel))

##### p__Proteobacteria, level 8 & salinity
d = data %>%
  filter(grepl("p__",species)) %>%
  filter(!grepl("c__",species)) %>% 
  group_by(sample) %>% 
  mutate(tot = sum(rel))

wilcox.test(d[d$species == 'd__Bacteria;p__Proteobacteria' & d$remark == 'control',]$rel,
            d[d$species == 'd__Bacteria;p__Proteobacteria' & d$remark == 'Level 8',]$rel)

wilcox.test(d[d$species == 'd__Bacteria;p__Proteobacteria' & d$remark == 'control',]$rel,
            d[d$species == 'd__Bacteria;p__Proteobacteria' & d$remark == 'salinity',]$rel)


wilcox.test(d[d$species == 'd__Bacteria;p__Proteobacteria' & d$remark == 'control',]$rel,
            d[d$species == 'd__Bacteria;p__Proteobacteria' & d$remark == 'Level 8' & d$sample %in% c(124, 127, 128),]$rel)

wilcox.test(d[d$species == 'd__Bacteria;p__Proteobacteria' & d$remark == 'control',]$rel,
            d[d$species == 'd__Bacteria;p__Proteobacteria' & d$remark == 'Level 8' & !d$sample %in% c(124, 127, 128),]$rel)


d = d %>% 
  filter(species == 'd__Bacteria;p__Proteobacteria')

ggplot(d)+
  geom_boxplot(aes(x = remark, y = rel))

ggplot(d %>% 
         filter(sample %in% c(124,127,128) | remark != "Level 8"))+
  geom_boxplot(aes(x = remark, y = rel))

ggplot(d %>% 
         filter(!sample %in% c(124,127,128) | remark != "Level 8"))+
  geom_boxplot(aes(x = remark, y = rel))

###### Firmicutes, salinity 
d = data %>%
  filter(grepl("p__",species)) %>%
  filter(!grepl("c__",species)) %>% 
  group_by(sample) %>% 
  mutate(tot = sum(rel))


wilcox.test(d[d$species == 'd__Bacteria;p__Firmicutes' & d$remark == 'control',]$rel,
            d[d$species == 'd__Bacteria;p__Firmicutes' & d$remark == 'salinity',]$rel)


##### Mycobacterium, 8 Factor 
d = data %>%
  filter(grepl("g__",species)) %>%
  filter(!grepl("s__",species)) %>% 
  group_by(sample) %>% 
  mutate(tot = sum(rel))


# ctr - level 8
wilcox.test(d[d$species == 'd__Bacteria;p__Actinobacteriota;c__Actinomycetia;o__Mycobacteriales;f__Mycobacteriaceae;g__Mycobacterium' & d$remark == 'control',]$rel,
            d[d$species == 'd__Bacteria;p__Actinobacteriota;c__Actinomycetia;o__Mycobacteriales;f__Mycobacteriaceae;g__Mycobacterium' & d$remark == 'Level 8',]$rel)

# salinity - level 8
wilcox.test(d[d$species == 'd__Bacteria;p__Actinobacteriota;c__Actinomycetia;o__Mycobacteriales;f__Mycobacteriaceae;g__Mycobacterium' & d$remark == 'salinity',]$rel,
            d[d$species == 'd__Bacteria;p__Actinobacteriota;c__Actinomycetia;o__Mycobacteriales;f__Mycobacteriaceae;g__Mycobacterium' & d$remark == 'Level 8',]$rel)

# salinity - control
wilcox.test(d[d$species == 'd__Bacteria;p__Actinobacteriota;c__Actinomycetia;o__Mycobacteriales;f__Mycobacteriaceae;g__Mycobacterium' & d$remark == 'salinity',]$rel,
            d[d$species == 'd__Bacteria;p__Actinobacteriota;c__Actinomycetia;o__Mycobacteriales;f__Mycobacteriaceae;g__Mycobacterium' & d$remark == 'control',]$rel)



# p__Actinobacteriota, Level 8, difference in salinity, non salinity
# (because of pcoa at the p level which indicates that p__Actinobacteriota 
#drives the difference among level 8 samples)
d = data %>%
  filter(grepl("p__",species)) %>%
  filter(!grepl("c__",species)) %>% 
  group_by(sample) %>% 
  mutate(tot = sum(rel))


head(d)
ggplot(d)+
  geom_point(aes(x = remark, y = tot, color = Salinity..9.))


####### load reference MAGs relative abundance data
data = read.table("coverM.rel_abs.tab",header = F, sep = '\t')
names(data) = c("genome","species","sample","rel","rel_tot")
data$sample = tstrsplit(data$sample, "_")[[2]]
data$sample = str_replace_all(data$sample,"b",'')
data = data %>%
  left_join(metadata,by = "sample")



# beta diversity, before blocking by salinity
data_d = data %>% 
  filter(sample != 85) %>% 
  filter(!grepl("UNCLASSIFIED",species)) %>% 
  mutate(sample = paste(sample,remark)) %>%
  dplyr::select(genome,rel,sample,remark) %>% 
  dplyr::select(-remark) %>% 
  spread(sample,rel)


data_d[is.na(data_d)] = 0
r = (data_d$genome)
data_d$genome = NULL
s = names(data_d)
t = data_d

# beta diversity differences across groups
d = vegdist(t(t),method = "bray")
g = as.factor(tstrsplit(s, " ")[[2]])
b = betadisper(d, g)
anova(b)

# permutation test on beta diversity
ptest = permutest(b, pairwise = TRUE, permutations = 999)
ptest 

# permutation test on beta diversity, after blocking by salinity
sal = data %>% 
  ungroup() %>% 
  filter(sample != 85) %>%
  filter(!grepl("UNCLASSIFIED",species)) %>% 
  mutate(sample = paste(sample,remark)) %>%
  dplyr::select(-genome) %>% 
  dplyr::select(Salinity..9.,sample) %>% 
  unique()

sal = sal %>% 
  arrange(sample,s)

?permutest
ptest = permutest(b, pairwise = TRUE, permutations = 999,strata = sal$Salinity..9.)
ptest


# alpha diversity
data_d = data %>% 
  filter(sample != 85) %>% 
  filter(!grepl("UNCLASSIFIED",species)) %>% 
  dplyr::select(genome,rel,sample,remark) %>% 
  mutate(sample = paste(sample,remark)) %>% 
  dplyr::select(-remark) %>% 
  spread(sample,rel)

data_d[is.na(data_d)] = 0
head(data_d)
r = (data_d$genome)
data_d$genome = NULL
s = names(data_d)
t = data_d

div = diversity(t,index = "shannon") # no jaccard 
div = data.frame(s,div)
div$treat = gsub("[0-9]* ","",div$s)
div = div %>% 
  mutate(lv = case_when(treat == "control" ~0,
                        treat == "Level8"~8,
                        .default = 1))



# relevance of each factor on diversity 
div$sample  = tstrsplit(div$s, " ")[[1]]
div = div %>% 
  left_join(metadata,by = "sample")

facts = c("Drought..1.", "N.dep..2.", "Temp..3.", "Microplastic..4.", "Glyphosate..5.", "Antibiotics.oxytetra..6.", "Fungicide..7.","Cu..8.", "Salinity..9.", "Insect.imidacloprid..10.", "Lv")

c_div = div %>% 
  dplyr::select(shannon,facts) %>% 
  mutate(Lv = case_when(Lv == 8~1,
                        .default = 0))

set.seed(123)
regressor <- randomForest(shannon ~ ., data = c_div, importance=TRUE) 
relImp = as.data.frame(importance(regressor))
relImp$treat = rownames(relImp)

ggplot(relImp)+
  geom_bar(aes(x = `%IncMSE`, y = treat),stat = "identity")

relImp

#####
# Viral  
##### 

setwd("~/analysis/Berlin/soil 2019 metagenomes/viral_MAGs/drep/coverM/")

# load relative abundance per bin
data = read.table("coverM.rel_abs.tab",header = F, sep = '\t')
names(data) = c("bin","species","sample","rel","rel_tot")
data$sample = tstrsplit(data$sample, "_")[[2]]
data$sample = str_replace_all(data$sample,"b",'')

data = data %>%
  left_join(metadata,by = "sample") %>%  
   filter(!is.na(species))

data[data$species=="NULL",]$species <- "Unclassified"

# unknown viral genomes abundance 
head(data)
d = data %>%
  group_by(sample) %>%
  mutate(tot = sum(rel)) %>%
  group_by(sample,species) %>%
  mutate(rel_sp = sum(rel)) %>%
  ungroup() %>% 
  dplyr::select(tot,sample,remark,rel_sp,species,Lv) %>% 
  unique()

dp = d %>% 
  filter(sample != "85 salinity" & species == 'Unclassified')
wilcox.test(dp[dp$Lv == 8,]$rel_sp ,dp[dp$Lv == 0,]$rel_sp)

# alpha diversity
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


div$sample =  tstrsplit(div$s, " ")[[1]]
div = div %>% 
  left_join(metadata,by = "sample")

# relevance of different factors on viral alpha diversity 
facts = c("Drought..1.", "N.dep..2.", "Temp..3.", "Microplastic..4.", "Glyphosate..5.", "Antibiotics.oxytetra..6.", "Fungicide..7.","Cu..8.", "Salinity..9.", "Insect.imidacloprid..10.", "Lv")

c_div = div %>% 
  dplyr::select(shannon,facts) %>% 
  mutate(Lv = case_when(Lv == 8~1,
                        .default = 0))

set.seed(123)
regressor <- randomForest(shannon ~ ., data = c_div, importance=TRUE) 
relImp = as.data.frame(importance(regressor))
relImp$treat = rownames(relImp)

ggplot(relImp)+
  geom_bar(aes(x = `%IncMSE`, y = treat),stat = "identity")

relImp

# plot variables rel imp
head(c_div)
ggplot(c_div)+
  geom_boxplot(aes(x = as.factor(Antibiotics.oxytetra..6.), y = shannon))

# beta diversity
data_d = data %>% 
  filter(sample != 85) %>% 
  dplyr::select(bin,rel,sample,remark) %>% 
  mutate(sample = paste(sample,remark)) %>% 
  dplyr::select(-remark) %>% 
  spread(sample,rel)

data_d[is.na(data_d)] = 0
head(data_d)
r = (data_d$bin)
data_d$bin = NULL
s = names(data_d)
t = data_d

d = vegdist(t(t),method = "bray")
g = as.factor(tstrsplit(s, " ")[[2]])
b = betadisper(d, g)
anova(b)

ptest = permutest(b, pairwise = TRUE, permutations = 99)
ptest


#####
# genes
#####

setwd("~/analysis/Berlin/soil 2019 metagenomes/gene_prediction/")

# load marker copy number per samplpe
mabs = read.table("../maker_gene_count/maker_count_per_sample.tab", header = F, sep = "\t")
names(mabs) = c("sample","mean_n_mar","median_n_mar","total_n_mar")
mabs$sample = tstrsplit(mabs$sample, "_")[[2]]
mabs$sample = str_replace_all(mabs$sample,"b",'')
mabs = mabs %>% 
  group_by(sample) %>% 
  dplyr::select(sample,mean_n_mar) %>% 
  unique()


# load data abs per ko per sample
data_rel = read.table("KOs/ko_abs_per_sample.desc.tab",sep = '\t',header = F,quote = "")
names(data_rel) = c("desc","ko","sample","count")
head(data_rel)
data_rel$sample = tstrsplit(data_rel$sample, "_")[[2]]
data_rel$sample = str_replace_all(data_rel$sample,"b",'')

# load data ko number genes per sample
data_c = read.table("KOs/ko_number_per_sample.desc.tab",header = F, sep = '\t',quote = "")
head(data_c)
names(data_c) = c("desc","ko","sample","n_genes")
data_c$sample = tstrsplit(data_c$sample, "_")[[2]]
data_c$sample = str_replace_all(data_c$sample,"b",'')
head(data_c)
data_c$desc = NULL

data = merge(data_c,data_rel,by = c("ko","sample"))
data = merge(data,metadata,by = "sample")
data = merge(data,mabs,by = "sample")
data$total_n_reads = 31310395

# respiration 
d = data %>% 
  filter(sample!=85) %>% 
  filter(desc %in% c("coxA, ctaD; cytochrome c oxidase subunit I [EC:7.1.1.9]",
                     "coxB, ctaC; cytochrome c oxidase subunit II [EC:7.1.1.9]",
                     "coxC, ctaE; cytochrome c oxidase subunit III [EC:7.1.1.9]",
                     "coxD, ctaF; cytochrome c oxidase subunit IV [EC:7.1.1.9]")) 

d$remark = factor(d$remark,levels = c("control","antibiotics", "copper", "drought",
                                      "microplastic", "Ndep", "salinity", "temp",
                                      "fungicide", "glyphosate","insecticide","Level 8"
))

d$desc =   substr(d$desc, 1, 4)

p1 = ggplot(d,aes(x = remark,y = count / total_n_reads))+
  theme_classic()+
  geom_boxplot()+
  theme(axis.text.x=element_text(angle=90, hjust=1))+
  facet_wrap(~desc,scales = 'free_x',ncol = 4)+
  stat_compare_means(label = "p.signif", method = "wilcox.test",
                     ref.group = "control",hide.ns = TRUE, tip.length = 0)+
  coord_flip()+
  ylab("Relative abundance")


p2 = ggplot(d,aes(x = remark,y = n_genes / mean_n_mar))+
  theme_classic()+
  geom_boxplot()+
  theme(axis.text.x=element_text(angle=90, hjust=1))+
  facet_wrap(~desc,scales = 'free_x',ncol = 4)+
  stat_compare_means(label = "p.signif", method = "wilcox.test",
                     ref.group = "control",hide.ns = TRUE, tip.length = 0)+
  coord_flip()+
  ylab("Copy number per cell")



p1 + p2 

dc = d %>% 
  group_by(sample) %>% 
  mutate(tot_copyn = mean(n_genes / mean_n_mar)) %>% 
  dplyr::select(sample,remark,tot_copyn) %>% 
  unique()

cnpc = ggplot(dc,aes(x = remark,y = tot_copyn))+
  theme_classic()+
  geom_boxplot()+
  theme(axis.text.x=element_text(angle=90, hjust=1))+
  #  facet_wrap(~desc,scales = 'free',ncol = 1)+
  stat_compare_means(label = "p.signif", method = "wilcox.test",
                     ref.group = "control",hide.ns = TRUE, tip.length = 0)+
  coord_flip()+
  ylab("Copy number per cell")

cnpc
wilcox.test(dc[dc$remark == "control",]$tot_copyn,dc[dc$remark == "Level 8",]$tot_copyn)

# corr
names(d)
p1_corr = ggplot(d,aes(x = CO2_ppm_ThirdWeek,y = count / total_n_reads))+
  theme_classic()+
  geom_smooth(method = "lm")+
  geom_point(aes(color = as.factor(Lv)))+
  theme(axis.text.x=element_text(angle=90, hjust=1))+
  ylab("Relative abundance")+
  facet_wrap(~desc,scales = "free")+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2")+
  labs(fill = "Number of factors",color = "Number of factors")

head(d)
p2_corr = ggplot(d,aes(x = CO2_ppm_ThirdWeek,y = n_genes / mean_n_mar))+
  theme_classic()+
  geom_smooth(method = "lm")+
  geom_point(aes(color = as.factor(Lv)))+
  theme(axis.text.x=element_text(angle=90, hjust=1))+
  ylab("Copy number per cell")+
  facet_wrap(~desc,scales = "free")+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2")+
  labs(fill = "Number of factors",color = "Number of factors")

p1_corr + p2_corr +  plot_layout(guides = "collect")

dc = dc %>% 
  left_join(metadata,by = "sample")
cor.test(dc$tot_copyn ,dc$CO2_ppm_ThirdWeek,method = "spearman")


# flg 
d = data %>% 
  filter(grepl("flg",desc))

d$remark = factor(d$remark,levels = c("control","antibiotics", "copper", "drought",
                                      "microplastic", "Ndep", "salinity", "temp",
                                      "fungicide", "glyphosate","insecticide","Level 8"
))

d$desc =   substr(d$desc, 1, 4)

p1 = ggplot(d,aes(x = remark,y = count / total_n_reads))+
  theme_classic()+
  geom_boxplot()+
  theme(axis.text.x=element_text(angle=90, hjust=1))+
  facet_wrap(~desc,scales = 'free_x',ncol = 4)+
  stat_compare_means(label = "p.signif", method = "wilcox.test",
                     ref.group = "control",hide.ns = TRUE, tip.length = 0)+
  coord_flip()+
  ylab("Relative abundance")


p2 = ggplot(d,aes(x = remark,y = n_genes / mean_n_mar))+
  theme_classic()+
  geom_boxplot()+
  theme(axis.text.x=element_text(angle=90, hjust=1))+
  facet_wrap(~desc,scales = 'free_x',ncol = 4)+
  stat_compare_means(label = "p.signif", method = "wilcox.test",
                     ref.group = "control",hide.ns = TRUE, tip.length = 0)+
  coord_flip()+
  ylab("Copy number per cell")



p1 + p2 

dc = d %>% 
  group_by(sample) %>% 
  mutate(tot_copyn = mean(n_genes / mean_n_mar)) %>% 
  dplyr::select(sample,remark,tot_copyn) %>% 
  unique()

cnpc = ggplot(dc,aes(x = remark,y = tot_copyn))+
  theme_classic()+
  geom_boxplot()+
  theme(axis.text.x=element_text(angle=90, hjust=1))+
  #  facet_wrap(~desc,scales = 'free',ncol = 1)+
  stat_compare_means(label = "p.signif", method = "wilcox.test",
                     ref.group = "control",hide.ns = TRUE, tip.length = 0)+
  coord_flip()+
  ylab("Copy number per cell")

head(dc)
wilcox.test(dc[dc$remark == "control",]$tot_copyn,dc[dc$remark == "Level 8",]$tot_copyn)


# corr
names(d)
p1_corr = ggplot(d,aes(x = WSA,y = count / total_n_reads))+
  theme_classic()+
  geom_smooth(method = "lm")+
  geom_point(aes(color = as.factor(Lv)))+
  theme(axis.text.x=element_text(angle=90, hjust=1))+
  ylab("Relative abundance")+
  facet_wrap(~desc,scales = "free")+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2")+
  labs(fill = "Number of factors",color = "Number of factors")

head(d)
p2_corr = ggplot(d,aes(x = WSA,y = n_genes / mean_n_mar))+
  theme_classic()+
  geom_smooth(method = "lm")+
  geom_point(aes(color = as.factor(Lv)))+
  theme(axis.text.x=element_text(angle=90, hjust=1))+
  ylab("Copy number per cell")+
  facet_wrap(~desc,scales = "free")+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2")+
  labs(fill = "Number of factors",color = "Number of factors")

p1_corr + p2_corr +  plot_layout(guides = "collect")

dc = dc %>% 
  left_join(metadata,by = "sample")
cor.test(dc$tot_copyn ,dc$WSA,method = "spearman")


# nitrification
d = data %>% 
  filter(ko %in% c("K10944", # pmoA-amoA; methane/ammonia monooxygenase subunit A
                   "K10945", #pmoB-amoB
                   "K10946"))# pmoC-amoC

d$remark = factor(d$remark,levels = c("control","antibiotics", "copper", "drought",
                                      "microplastic", "Ndep", "salinity", "temp",
                                      "fungicide", "glyphosate","insecticide","Level 8"
))


d$desc =   substr(d$desc, 1, 4)
p1 = ggplot(d,aes(x = remark,y = count / total_n_reads))+
  theme_classic()+
  geom_boxplot()+
  theme(axis.text.x=element_text(angle=90, hjust=1))+
  facet_wrap(~desc,scales = 'free',ncol = 1)+
  stat_compare_means(label = "p.signif", method = "wilcox.test",
                     ref.group = "control",hide.ns = TRUE, tip.length = 0)+
  coord_flip()+
  ylab("Relative abundance")


p2 = ggplot(d,aes(x = remark,y = n_genes / mean_n_mar))+
  theme_classic()+
  geom_boxplot()+
  theme(axis.text.x=element_text(angle=90, hjust=1))+
  facet_wrap(~desc,scales = 'free',ncol = 1)+
  stat_compare_means(label = "p.signif", method = "wilcox.test",
                     ref.group = "control",hide.ns = TRUE, tip.length = 0)+
  coord_flip()+
  ylab("Copy number per cell")


p1 + p2

dc = d %>% 
  group_by(sample) %>% 
  mutate(tot_copyn = mean(n_genes / mean_n_mar)) %>% 
  dplyr::select(sample,remark,tot_copyn) %>% 
  unique()

ggplot(dc,aes(x = remark,y = tot_copyn))+
  theme_classic()+
  geom_boxplot()+
  theme(axis.text.x=element_text(angle=90, hjust=1))+
  #  facet_wrap(~desc,scales = 'free',ncol = 1)+
  stat_compare_means(label = "p.signif", method = "wilcox.test",
                     ref.group = "control",hide.ns = TRUE, tip.length = 0)+
  coord_flip()+
  ylab("Copy number per cell")

wilcox.test(dc[dc$remark == 'control',]$tot_copyn,dc[dc$remark == 'Level 8',]$tot_copyn)

# Organic nitrogen metabolism
#####
# Organic N metabolism
#####


# https://www.biorxiv.org/content/10.1101/2020.06.23.167700v1.full

d = data %>% 
  filter(ko %in% c("K19823",# NAO
                   "K01915"))# glna 

d$remark = factor(d$remark,levels = c("control","antibiotics", "copper", "drought",
                                      "microplastic", "Ndep", "salinity", "temp",
                                      "fungicide", "glyphosate","insecticide","Level 8"
))


d$desc = substr(d$desc, 1, 4)

p1 = ggplot(d,aes(x = remark,y = count / total_n_reads))+
  theme_classic()+
  geom_boxplot()+
  theme(axis.text.x=element_text(angle=90, hjust=1))+
  facet_wrap(~desc,scales = 'free_x',ncol = 2)+
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "control",hide.ns = TRUE, tip.length = 0)+
  coord_flip()+
  ylab("Relative abundance")


p2 = ggplot(d,aes(x = remark,y = n_genes / mean_n_mar))+
  theme_classic()+
  geom_boxplot()+
  theme(axis.text.x=element_text(angle=90, hjust=1))+
  facet_wrap(~desc,scales = 'free_x',ncol = 2)+
  stat_compare_means(label = "p.signif", method = "wilcox.test",
                     ref.group = "control",hide.ns = TRUE, tip.length = 0)+
  coord_flip()+
  ylab("Copy number per cell")


p1 + p2


# combined into one 
## copy n per cell
dc = d %>% 
  group_by(sample) %>% 
  mutate(tot_copyn = mean(n_genes / mean_n_mar)) %>% 
  dplyr::select(sample,remark,tot_copyn) %>% 
  unique()

cnpc = ggplot(dc,aes(x = remark,y = tot_copyn))+
  theme_classic()+
  geom_boxplot()+
  theme(axis.text.x=element_text(angle=90, hjust=1))+
  #  facet_wrap(~desc,scales = 'free',ncol = 1)+
  stat_compare_means(label = "p.signif", method = "wilcox.test",
                     ref.group = "control",hide.ns = TRUE, tip.length = 0)+
  coord_flip()+
  ylab("Copy number per cell")

## rel abs 
head(d)
dc = d %>% 
  group_by(sample) %>% 
  mutate(tot_copyn = mean(count / total_n_reads)) %>% 
  dplyr::select(sample,remark,tot_copyn) %>% 
  unique()

rel_abs = ggplot(dc,aes(x = remark,y = tot_copyn))+
  theme_classic()+
  geom_boxplot()+
  theme(axis.text.x=element_text(angle=90, hjust=1))+
  #  facet_wrap(~desc,scales = 'free',ncol = 1)+
  stat_compare_means(label = "p.signif", method = "wilcox.test",
                     ref.group = "control",hide.ns = TRUE, tip.length = 0)+
  coord_flip()+
  ylab("Relative abundance")

cnpc + rel_abs


# cazymes
data = read.table("CAZY/cazy_number_per_sample.tab",header = F, sep = '\t',quote = "")
head(data)
names(data) = c("ko","sample","n_genes","desc")
data$sample = tstrsplit(data$sample, "_")[[2]]
data$sample = str_replace_all(data$sample,"b",'')
head(data)

data = data %>% 
  left_join(metadata,by = "sample") %>% 
  left_join(mabs,by = "sample")
nrow(data)


data$remark = factor(data$remark,levels = c("control","antibiotics", "copper", "drought",
                                            "microplastic", "Ndep", "salinity", "temp",
                                            "fungicide", "glyphosate","insecticide","Level 8"
))



d = data %>% 
  group_by(sample) %>% 
  mutate(tot = sum(n_genes)) %>% 
  dplyr::select(tot,remark,sample,mean_n_mar) %>% 
  unique()

ggplot(d,aes(x = remark,y = tot / mean_n_mar))+
  theme_classic()+
  geom_boxplot()+
  theme(axis.text.x=element_text(angle=90, hjust=1))+
  #  facet_wrap(~desc,scales = 'free')+
  stat_compare_means(label = "p.signif", method = "wilcox.test",
                     ref.group = "control",hide.ns = TRUE, tip.length = 0)+
  ylab("CAZy gene copy per cell")

wilcox.test(d[d$remark=='Level 8',]$tot / d[d$remark=='Level 8',]$mean_n_mar,d[d$remark=='control',]$tot)

# correlations

# respiration
names(data)
d = data %>% 
  filter(sample != 85) %>% 
  filter(ko %in% c("K02274","K02275","K02276","K02277")) %>% 
  group_by(sample) %>% 
  mutate(n = mean(n_genes / mean_n_mar)) %>% 
  mutate(n2 = mean( count / total_n_reads)) %>% 
  dplyr::select(CO2_ppm_ThirdWeek,n,Lv,n2,remark) %>% 
  unique()


cor.test(d$CO2_ppm_ThirdWeek,d$n)

correlation_coef <- d %>%
  ungroup() %>% 
  summarize(correlation = cor(CO2_ppm_ThirdWeek, n,method = "spearman"),
            p_value = cor.test(CO2_ppm_ThirdWeek, n, method = "spearman")$p.value)

correlation_coef
corr_resp = ggplot(d,aes(x = CO2_ppm_ThirdWeek,y = n2))+
  theme_classic()+
  geom_smooth(method = "lm")+
  geom_point(aes(shape = as.factor(Lv),color = remark),size = 2)+
  theme(axis.text.x=element_text(angle=90, hjust=1))+
  ylab("Copy number per cell")+
  labs(color = "Treatment",shape = "Number of factors")+
  ggtitle("cox operon")+
  xlab("CO2 ppm")+
  scale_color_manual(values = custom_palette)

corr_resp

# flg
d = data %>% 
  filter(grepl("flg",desc))

d$desc =substr(d$desc, 1, 4)


d = d %>% 
  group_by(sample) %>% 
  mutate(n = sum(n_genes / mean_n_mar)) %>% 
  mutate(n2 = mean( count / total_n_reads)) %>% 
  dplyr::select(WSA,n,Lv,n2,remark) %>% 
  unique()

correlation_coef <- d %>%
  ungroup() %>% 
  summarize(correlation = cor(WSA, n,method = "spearman"),
            p_value = cor.test(WSA, n,method = "spearman")$p.value)


correlation_coef
corr_flg = ggplot(d,aes(x = WSA,y = n2))+
  theme_classic()+
  geom_smooth(method = "lm")+
  geom_point(aes(shape = as.factor(Lv),color = remark),size = 2)+
  theme(axis.text.x=element_text(angle=90, hjust=1))+
  ylab("Copy number per cell")+
  labs(color = "Treatment",shape = "Number of factors")+
  ggtitle("Flg operon")+
  scale_color_manual(values = custom_palette)


corr_flg

#####
# AMR
#####
setwd("~/analysis/Berlin/soil 2019 metagenomes/gene_prediction/")

# load data
data = read.table("CARD/hits.1e7_id80_cov75.n_genes_ass.gene_name.CARD_code.tab",header = F, sep = '\t',quote = "")
names(data) = c("sample","gene","cat","name","n","ncard","ncard_plasmids","nmob")
data$sample = tstrsplit(data$sample, "_")[[2]]
data$sample = str_replace_all(data$sample,"b",'')

data = data %>% 
  left_join(metadata,by = "sample") %>% 
  left_join(mabs,by = "sample")

# load composition per sample
comp = read.table("../MAGs/semibin/dRep/coverM/composition.no85.tab",header = T)
comp$sample = tstrsplit(comp$sample, " ")[[1]]
comp$sample = str_replace_all(comp$sample,"b",'')
comp = comp %>% 
  #  filter(dataset == "SILVA" & tax_group == "Bacteria") %>% 
  dplyr::select(sample,V1)
data = data %>% 
  left_join(comp,by = "sample")


data$remark = factor(data$remark,levels = c("control","antibiotics", "copper", "drought",
                                            "microplastic", "Ndep", "salinity", "temp",
                                            "fungicide", "glyphosate","insecticide","Level 8"
))


# copies per cell, per category general 

d = data %>% 
  filter(sample != 85) %>% 
  group_by(sample,cat) %>% 
  mutate(n  = sum( ncard / mean_n_mar)) %>% 
  dplyr::select(n,sample,cat,Lv,remark,Salinity..9.) %>% 
  unique()

dcontrol = d %>% 
  filter(remark == 'control') %>%
  group_by(cat) %>% 
  mutate(median_control = median(n)) %>% 
  dplyr::select(cat,median_control) %>% 
  ungroup() %>% 
  unique()


d = d %>% 
  left_join(dcontrol,by = "cat")

ggplot(d,
       aes(x = remark, y = n))+
  geom_boxplot(aes(color = as.factor(Lv),fill = as.factor(Lv)),outlier.shape = NA,alpha = 0.3)+
  coord_flip()+
  facet_wrap(~cat,scales = "free_x",ncol = 2)+
  stat_compare_means(label = "p.signif", method = "wilcox.test",
                     ref.group = "control",hide.ns = TRUE, tip.length = 0)+
  theme_classic()+
  geom_jitter(aes(color = as.factor(Lv),fill = as.factor(Lv)),shape=16, position=position_jitter(0.2),alpha = 0.8)+
  ylab("Copies per cell")+
  theme_minimal()+
  geom_hline(aes(yintercept= median_control), linetype="dashed", 
             color = "grey", size=1)+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2")+
  labs(fill = "Number of factors",color = "Number of factors")




# particular genes
dcontrol = data %>% 
  mutate(n  =  ncard / mean_n_mar) %>% 
  filter(remark == 'control') %>%
  group_by(name) %>% 
  mutate(median_control = median(n)) %>% 
  dplyr::select(name,median_control) %>% 
  ungroup() %>% 
  unique()

d = data %>% 
  left_join(dcontrol,by = "name")

ggplot(d %>% 
         filter(name %in% c("BJP-1","RbpA","Erm(37)","efpA","vanR_in_vanO_cl")) %>% 
         filter(sample != 85),aes(x = remark, y = ncard / mean_n_mar))+
  geom_boxplot(aes(color = as.factor(Lv),fill = as.factor(Lv)),outlier.shape = NA,alpha = 0.3)+
  coord_flip()+
  facet_wrap(~name,scales = "free_x",ncol = 2)+
  stat_compare_means(label = "p.signif", method = "wilcox.test",
                     ref.group = "control",hide.ns = TRUE, tip.length = 0)+
  theme_classic()+
  geom_jitter(aes(color = as.factor(Lv),fill = as.factor(Lv)),shape=16, position=position_jitter(0.2),alpha = 0.8)+
  ylab("Copies per cell")+
  theme_minimal()+
  geom_hline(aes(yintercept= median_control), linetype="dashed", 
             color = "grey", size=1)+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2")+
  labs(fill = "Number of factors",color = "Number of factors")



# corr with composition 
# load collapsed per sample
data = read.table("~/analysis/Berlin/soil 2019 metagenomes/gene_prediction/stats_per_sample.1e7_id80_cov75.tab",header = F, sep = '\t')
names(data) = c("sample","n_genes","ncard","ncard_plasmids","div_card","nres","nres_plasmids","div_res","nmob","div_mob")
data$sample = tstrsplit(data$sample, "_")[[2]]
data$sample = str_replace_all(data$sample,"b",'')
head(data)

data = data %>% 
  left_join(metadata,by = "sample") %>% 
  left_join(mabs,by = "sample") %>% 
  left_join(comp,by = "sample") %>%
  filter(sample != 85)

correlation_coef <- data %>%
  ungroup() %>% 
  summarize(correlation = cor(V1, ncard / mean_n_mar,method = "spearman"),
            p_value = cor.test(V1, ncard / mean_n_mar,method = "spearman")$p.value)


correlation_coef
custom_palette <- c("#999999", "#377EB8", "#F781BF","#4DAF4A", "#FF7F00",
                    "#E41A1C", "#A65628", "#984EA3", "#FFFF33", "#66C2A5", "black", "#8EBA42")

data$remark = factor(data$remark,levels = c("control","antibiotics", "copper", "drought",
                                            "microplastic", "Ndep", "salinity", "temp",
                                            "fungicide", "glyphosate","insecticide","Level 8"
))

ggplot(data,aes(y = ncard / mean_n_mar, x = remark))+
  geom_boxplot(aes(color = remark))+
  stat_compare_means(aes(group=remark,y = ncard / mean_n_mar, x = remark),label = "p.signif", method = "wilcox.test",
                     ref.group = "Level 8",hide.ns = TRUE, tip.length = 0)


pcorr = ggplot(data ,aes(x = ncard / mean_n_mar, y = V1))+
  geom_point(alpha = 0.5,size = 3,aes(color = remark))+
  geom_smooth(method = 'lm')+
  theme_minimal()+
  ylab("Bacterial composition")+
  xlab("ARG copy number per cell")+
  scale_color_manual(values = custom_palette)

pcorr

correlation_coef <- data %>%
  ungroup() %>% 
  summarize(correlation = cor(nmob / mean_n_mar, ncard / mean_n_mar),
            p_value = cor.test(nmob / mean_n_mar, ncard / mean_n_mar)$p.value)

correlation_coef
pcorr2 = ggplot(data,aes(x = ncard / mean_n_mar, y = nmob / mean_n_mar))+
  geom_point(alpha = 0.5,size = 3,aes(color = remark))+
  geom_smooth(method = 'lm')+
  theme_minimal()+
  ylab("MGE copy number per cell")+
  xlab("ARG copy number per cell")+
  scale_color_manual(values = custom_palette)
pcorr2
