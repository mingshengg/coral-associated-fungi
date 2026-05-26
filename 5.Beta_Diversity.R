# Load packages, data, and customizations ####
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(vegan); packageVersion("vegan")
library(patchwork); packageVersion("patchwork")
library(microbiome); packageVersion("microbiome")
library(broom); packageVersion("broom")
library(purrr); packageVersion("purrr")
library(corncob); packageVersion("corncob")
library(pairwiseAdonis); packageVersion("pairwiseAdonis")
library(ggfortify); packageVersion("ggfortify")
library(metagMisc); packageVersion("metagMisc")

se <- function(x) {
  sd(n)/sqrt(length(n))
}

# read in phyloseq ####
ps <- readRDS('./Output/ps_assigned.RDS') %>%
  subset_taxa(Kingdom == 'k__Fungi') %>%
  subset_samples(Individual != 'NA') %>%
  subset_samples(Compartment == 'Tissue'| Compartment == 'Skeleton') 

ps <- ps %>% subset_taxa(taxa_sums(ps)>0)

# metadata file
full_meta = ps@sam_data %>% data.frame()
mdf_meta <- full_meta[7:11] %>% mutate_all(., function(x) {as.numeric(x)})

# sampling depth, rarefaction
otu_tab <- ps@otu_table %>% as.data.frame()
sampling_depth <- ps@otu_table %>% as.data.frame() %>% rowSums() %>% min()

mdf_com <- avgdist(otu_tab, iterations = 999, sample = sampling_depth, dmethod = 'bray')

#### Exploring overall composition with PCA ####
otu_res <- metaMDS(mdf_com)
otu_res

data.scores = as.data.frame(scores(otu_res, display =  'sites'))
rownames(data.scores) == rownames(full_meta)

data.scores$Site = full_meta$Site
data.scores$Compartment = full_meta$Compartment
data.scores$Species = full_meta$Species

# Plot NMDS
gg = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores, aes(colour = Site, shape = Species), size = 3, alpha = 1) +
  geom_polygon(stat = 'ellipse', aes(fill = Site), alpha = 0.3)  +
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) + 
  labs(title = round(otu_res$stress,5)) +
  xlab('NMDS1') +
  ylab('NMDS2')

gg

#PERMANOVA to investigate difference in composition ####
full_meta <- full_meta[order(rownames(full_meta)),]
mdf_com = as.matrix(mdf_com, labels=TRUE)
rownames(mdf_com) == rownames(full_meta) # ensure that formatting of metadata table and otu table is the same

# full model with bray-curtis distances
full_model <- adonis2(mdf_com ~ full_meta$Species*full_meta$Compartment*full_meta$Site, by = 'terms')
full_model

# betadisper to assess homogeneity of group dispersion
bd <- betadisper(vegdist(mdf_com, 'bray'), group = full_meta$Compartment)
anova(bd)
plot(bd)
boxplot(bd)
permutest(bd)

# which pairs are significantly different
bd_pair <- TukeyHSD(bd, p.adjust = 'fdr')$group
bd_pair %>% as.data.frame() %>% filter("p adj" < 0.05)

# which sites are significantly different?
pairwise.adonis2(mdf_com ~ Site, data = full_meta)

#### Species level compositional analysis ####
#### Pocillopora acuta ####
ps_poci <- readRDS('./Output/ps_assigned.RDS') %>%
  subset_samples(Species != 'Pachyseris_speciosa' & Species != 'Diploastrea_heliopora') %>%
  subset_samples(Species != 'BLANK' & Species != 'Mock') %>%
  subset_taxa(Kingdom == 'k__Fungi')

poci_meta <- ps_poci@sam_data %>% data.frame()
poci_mdf_meta <- poci_meta[7:11] %>% mutate_all(., function(x) {as.numeric(x)})
poci_otu <- ps_poci@otu_table %>% as.data.frame()

poci_com <- avgdist(poci_otu, iterations = 999, sample = poci_otu %>% rowSums() %>% min(), dmethod = 'bray')

otu_res_poci <- metaMDS(poci_com)
otu_res_poci

data.scores = as.data.frame(scores(otu_res_poci, display =  'sites'))
rownames(data.scores) == rownames(poci_meta)

data.scores$Site = poci_meta$Site
data.scores$Compartment = poci_meta$Compartment

# Plot NMDS
gg_poci = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores, aes(colour = Compartment, shape = Site), size = 3, alpha = 1) +
  geom_polygon(stat = 'ellipse', aes(fill = Compartment), alpha = 0.3)  +
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) + 
  labs(title = paste0("Pocillopora acuta; ", round(otu_res_poci$stress,4))) +
  xlab('NMDS1') +
  ylab('NMDS2')

gg_poci

# poci model with bray-curtis distances
poci_model <- adonis2(poci_com ~ poci_meta$Compartment*poci_meta$Site, by = 'terms')
poci_model

bd <- betadisper(vegdist(poci_com, 'bray'), group = poci_meta$Compartment)
anova(bd)
boxplot(bd)
TukeyHSD(bd)

#### Pachyseris speciosa ####
ps_pachy <- readRDS('./Output/ps_assigned.RDS') %>%
  subset_samples(Species != 'Pocillopora_acuta' & Species != 'Diploastrea_heliopora') %>%
  subset_samples(Species != 'BLANK' & Species != 'Mock') %>%
  subset_taxa(Kingdom == 'k__Fungi')

pachy_meta <- ps_pachy@sam_data %>% data.frame()
pachy_mdf_meta <- pachy_meta[7:11] %>% mutate_all(., function(x) {as.numeric(x)})
pachy_otu <- ps_pachy@otu_table %>% as.data.frame()

pachy_com <- avgdist(pachy_otu, iterations = 999, sample = pachy_otu %>% rowSums() %>% min(), dmethod = 'bray')

otu_res_pachy <- metaMDS(pachy_com)
otu_res_pachy

data.scores = as.data.frame(scores(otu_res_pachy, display =  'sites'))
rownames(data.scores) == rownames(pachy_meta)

data.scores$Site = pachy_meta$Site
data.scores$Compartment = pachy_meta$Compartment

# Plot NMDS
gg_pachy = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores, aes(colour = Compartment, shape = Site), size = 3, alpha = 1) +
  geom_polygon(stat = 'ellipse', aes(fill = Compartment), alpha = 0.3)  +
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) + 
  labs(title = paste0("Pachyseris speciosa; ", round(otu_res_pachy$stress,4))) +
  xlab('NMDS1') +
  ylab('NMDS2')

gg_pachy

# pachy model with bray-curtis distances
pachy_model <- adonis2(pachy_com ~ pachy_meta$Compartment*pachy_meta$Site, by = 'terms')
pachy_model

bd <- betadisper(vegdist(pachy_com, 'bray'), group = pachy_meta$Compartment)
anova(bd)
boxplot(bd)
TukeyHSD(bd)
plot(bd)

#### Diploastrea heliopora ####
ps_diplo <- readRDS('./Output/ps_assigned.RDS') %>%
  subset_samples(Species != 'Pachyseris_speciosa' & Species != 'Pocillopora_acuta') %>%
  subset_samples(Species != 'BLANK' & Species != 'Mock') %>%
  subset_taxa(Kingdom == 'k__Fungi')

diplo_meta <- ps_diplo@sam_data %>% data.frame()
diplo_mdf_meta <- diplo_meta[7:11] %>% mutate_all(., function(x) {as.numeric(x)})
diplo_otu <- ps_diplo@otu_table %>% as.data.frame()

diplo_com <- avgdist(diplo_otu, iterations = 999, sample = diplo_otu %>% rowSums() %>% min(), dmethod = 'bray')

otu_res_diplo <- metaMDS(diplo_com)
otu_res_diplo

data.scores = as.data.frame(scores(otu_res_diplo, display =  'sites'))
rownames(data.scores) == rownames(diplo_meta)

data.scores$Site = diplo_meta$Site
data.scores$Compartment = diplo_meta$Compartment

# Plot NMDS
gg_diplo = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores, aes(colour = Compartment, shape = Site), size = 3, alpha = 1) +
  geom_polygon(stat = 'ellipse', aes(fill = Compartment), alpha = 0.3)  +
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) + 
  labs(title = paste0('Diploastrea heliopora; ', round(otu_res_diplo$stress,4))) +
  xlab('NMDS1') +
  ylab('NMDS2')

gg_diplo

# diplo model with bray-curtis distances
diplo_model <- adonis2(diplo_com ~ diplo_meta$Compartment*diplo_meta$Site, by = 'terms')
diplo_model

bd <- betadisper(vegdist(diplo_com, 'bray'), group = diplo_meta$Compartment)
anova(bd)
boxplot(bd)
TukeyHSD(bd)
plot(bd)


pairwise.adonis2(diplo_com ~ Compartment, data = diplo_meta, strata = 'Site')
pairwise.adonis2(pachy_com ~ Compartment, data = pachy_meta, strata = 'Site')
pairwise.adonis2(poci_com ~ Compartment, data = poci_meta, strata = 'Site')