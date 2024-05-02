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

# load metadata
metadata = read.table("metadata.tab",header = T,sep = '\t')
metadata$sample = paste(metadata$Tube.number,sep  = '')

######
# Bacteria taxonomic analysis
######

# difference in relative abundance of particular lineages
data = read.table("coverM.rel_abs.per_tax.tab",header = F, sep = '\t')
names(data) = c("species","sample","rel")
data$sample = tstrsplit(data$sample, "_")[[2]]
data$sample = str_replace_all(data$sample,"b",'')
data = data %>%
  left_join(metadata,by = "sample") %>% 
  filter(sample != 85)

### c__Actinomycetia
d = data %>%
  filter(grepl("c__",species)) %>%
  filter(!grepl("o__",species)) %>% 
  group_by(sample) %>% 
  mutate(tot = sum(rel))

# level 8 vs control
wilcox.test(d[d$species == 'd__Bacteria;p__Actinobacteriota;c__Actinomycetia' & d$remark == 'control',]$rel,
            d[d$species == 'd__Bacteria;p__Actinobacteriota;c__Actinomycetia' & d$remark == 'Level 8',]$rel)

# salinity vs control
wilcox.test(d[d$species == 'd__Bacteria;p__Actinobacteriota;c__Actinomycetia' & d$remark == 'control',]$rel,
            d[d$species == 'd__Bacteria;p__Actinobacteriota;c__Actinomycetia' & d$remark == 'salinity',]$rel)


### p__Actinobacteriota
d = data %>%
  filter(grepl("p__",species)) %>%
  filter(!grepl("c__",species)) %>% 
  group_by(sample) %>% 
  mutate(tot = sum(rel))

# level 8 vs control
wilcox.test(d[d$species == 'd__Bacteria;p__Actinobacteriota' & d$remark == 'control',]$rel,
            d[d$species == 'd__Bacteria;p__Actinobacteriota' & d$remark == 'Level 8',]$rel)


# level 8 (samples not including salinity) vs control
wilcox.test(d[d$species == 'd__Bacteria;p__Actinobacteriota' & d$remark == 'control',]$rel,
            d[d$species == 'd__Bacteria;p__Actinobacteriota' & d$remark == 'Level 8' & d$sample %in% c(124, 127, 128),]$rel)

# level 8 (samples including salinity) vs control
wilcox.test(d[d$species == 'd__Bacteria;p__Actinobacteriota' & d$remark == 'control',]$rel,
            d[d$species == 'd__Bacteria;p__Actinobacteriota' & d$remark == 'Level 8' & !d$sample %in% c(124, 127, 128),]$rel)


##### p__Proteobacteria
d = data %>%
  filter(grepl("p__",species)) %>%
  filter(!grepl("c__",species)) %>% 
  group_by(sample) %>% 
  mutate(tot = sum(rel))

# level 8 vs control
wilcox.test(d[d$species == 'd__Bacteria;p__Proteobacteria' & d$remark == 'control',]$rel,
            d[d$species == 'd__Bacteria;p__Proteobacteria' & d$remark == 'Level 8',]$rel)

# salinity 8 vs control
wilcox.test(d[d$species == 'd__Bacteria;p__Proteobacteria' & d$remark == 'control',]$rel,
            d[d$species == 'd__Bacteria;p__Proteobacteria' & d$remark == 'salinity',]$rel)

# level 8 (samples not including salinity) vs control
wilcox.test(d[d$species == 'd__Bacteria;p__Proteobacteria' & d$remark == 'control',]$rel,
            d[d$species == 'd__Bacteria;p__Proteobacteria' & d$remark == 'Level 8' & d$sample %in% c(124, 127, 128),]$rel)

# level 8 (samples including salinity) vs control
wilcox.test(d[d$species == 'd__Bacteria;p__Proteobacteria' & d$remark == 'control',]$rel,
            d[d$species == 'd__Bacteria;p__Proteobacteria' & d$remark == 'Level 8' & !d$sample %in% c(124, 127, 128),]$rel)


###### Firmicutes
d = data %>%
  filter(grepl("p__",species)) %>%
  filter(!grepl("c__",species)) %>% 
  group_by(sample) %>% 
  mutate(tot = sum(rel))

# salinity vs control
wilcox.test(d[d$species == 'd__Bacteria;p__Firmicutes' & d$remark == 'control',]$rel,
            d[d$species == 'd__Bacteria;p__Firmicutes' & d$remark == 'salinity',]$rel)


##### Mycobacterium
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



# load reference MAGs relative abundance data
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

relImp

#####
# Viral taxonomic analysis
##### 

# load relative abundance per bin
data = read.table("coverM.rel_abs.tab",header = F, sep = '\t')
names(data) = c("bin","species","sample","rel","rel_tot")
data$sample = tstrsplit(data$sample, "_")[[2]]
data$sample = str_replace_all(data$sample,"b",'')

data = data %>%
  left_join(metadata,by = "sample")

data[data$species=="NULL",]$species <- "Unclassified"

# unknown viral genomes abundance 
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

# viral alpha diversity
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

# beta diversity
data_d = data %>% 
  filter(sample != 85) %>% 
  dplyr::select(bin,rel,sample,remark) %>% 
  mutate(sample = paste(sample,remark)) %>% 
  dplyr::select(-remark) %>% 
  spread(sample,rel)

data_d[is.na(data_d)] = 0
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
# Functional gene analysis
#####

# load marker copy number per samplpe
mabs = read.table("../maker_gene_count/maker_count_per_sample.tab", header = F, sep = "\t")
names(mabs) = c("sample","mean_n_mar","median_n_mar","total_n_mar")
mabs$sample = tstrsplit(mabs$sample, "_")[[2]]
mabs$sample = str_replace_all(mabs$sample,"b",'')
mabs = mabs %>% 
  group_by(sample) %>% 
  dplyr::select(sample,mean_n_mar) %>% 
  unique()

# load ko number of genes per sample
data_c = read.table("KOs/ko_number_per_sample.desc.tab",header = F, sep = '\t',quote = "")
names(data_c) = c("desc","ko","sample","n_genes")
data_c$sample = tstrsplit(data_c$sample, "_")[[2]]
data_c$sample = str_replace_all(data_c$sample,"b",'')

data = merge(data_c,data_rel,by = c("ko","sample"))
data = merge(data,metadata,by = "sample")
data = merge(data,mabs,by = "sample")

##### respiration 
d = data %>% 
  filter(sample!=85) %>% 
  filter(desc %in% c("coxA, ctaD; cytochrome c oxidase subunit I [EC:7.1.1.9]",
                     "coxB, ctaC; cytochrome c oxidase subunit II [EC:7.1.1.9]",
                     "coxC, ctaE; cytochrome c oxidase subunit III [EC:7.1.1.9]",
                     "coxD, ctaF; cytochrome c oxidase subunit IV [EC:7.1.1.9]")) 

dc = d %>% 
  group_by(sample) %>% 
  mutate(tot_copyn = mean(n_genes / mean_n_mar)) %>% 
  dplyr::select(sample,remark,tot_copyn) %>% 
  unique()


# Level 8 vs control
wilcox.test(dc[dc$remark == "control",]$tot_copyn,dc[dc$remark == "Level 8",]$tot_copyn)

# corr
dc = dc %>% 
  left_join(metadata,by = "sample")

cor.test(dc$tot_copyn ,dc$CO2_ppm_ThirdWeek,method = "spearman")


#### flg 
d = data %>% 
  filter(grepl("flg",desc))

dc = d %>% 
  group_by(sample) %>% 
  mutate(tot_copyn = mean(n_genes / mean_n_mar)) %>% 
  dplyr::select(sample,remark,tot_copyn) %>% 
  unique()

# Level 8 vs control
wilcox.test(dc[dc$remark == "control",]$tot_copyn,dc[dc$remark == "Level 8",]$tot_copyn)


# corr
dc = dc %>% 
  left_join(metadata,by = "sample")
cor.test(dc$tot_copyn ,dc$WSA,method = "spearman")


#####
# AMR gene analysis
#####

# load data
data = read.table("CARD/hits.1e7_id80_cov75.n_genes_ass.gene_name.CARD_code.tab",header = F, sep = '\t',quote = "")
names(data) = c("sample","gene","cat","name","n","ncard","ncard_plasmids","nmob")
data$sample = tstrsplit(data$sample, "_")[[2]]
data$sample = str_replace_all(data$sample,"b",'')

data = data %>% 
  left_join(metadata,by = "sample") %>% 
  left_join(mabs,by = "sample")

# load composition per sample
comp = read.table("MAGs/semibin/dRep/coverM/composition.no85.tab",header = T)
comp$sample = tstrsplit(comp$sample, " ")[[1]]
comp$sample = str_replace_all(comp$sample,"b",'')
comp = comp %>% 
  dplyr::select(sample,V1)

data = data %>% 
  left_join(comp,by = "sample")

# copies per cell, per category, general 
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

correlation_coef_composition <- data %>%
  ungroup() %>% 
  summarize(correlation = cor(V1, ncard / mean_n_mar,method = "spearman"),
            p_value = cor.test(V1, ncard / mean_n_mar,method = "spearman")$p.value)


correlation_coef_composition

correlation_coef_MGEs <- data %>%
  ungroup() %>% 
  summarize(correlation = cor(nmob / mean_n_mar, ncard / mean_n_mar),
            p_value = cor.test(nmob / mean_n_mar, ncard / mean_n_mar)$p.value)

correlation_coef_MGEs
