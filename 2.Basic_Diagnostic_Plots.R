# Load packages ####
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(vegan); packageVersion("vegan")
library(speedyseq); packageVersion("speedyseq")

# Summary plots ####

# Plot of taxon-level assignment efficiency 
ps_sp <- readRDS("./Output/ps_fungi_rare_tree.RDS") %>% subset_taxa(Kingdom == 'k__Fungi')
metadata <- sample_data(ps_sp) %>% data.frame()

sample_data(ps_sp)$Sequenced <- NULL

phy <- tax_table(ps_sp)[,2] %>% data.frame() != 'NA'
cla <- tax_table(ps_sp)[,3] %>% data.frame() != 'NA'
ord <- tax_table(ps_sp)[,4] %>% data.frame() != 'NA'
fam <- tax_table(ps_sp)[,5] %>% data.frame() != 'NA'
gen <- tax_table(ps_sp)[,6] %>% data.frame() != 'NA'
spp <- tax_table(ps_sp)[,7] %>% data.frame() != 'NA'
assignments <- data.frame(Phylum=phy, Class=cla,Order=ord,Family=fam,Genus=gen,Species=spp)

assignments %>% pivot_longer(1:6) %>% mutate(name=factor(name,levels = c("Phylum","Class","Order","Family","Genus","Species"))) %>%
  ggplot(aes(x=name,fill=value)) + geom_bar() + scale_fill_manual(values=c("Gray","Black")) +
  labs(x="Taxonomic level",y="Count",fill="Unambiguous\nassignment")

rm(king,phy,cla,ord,fam,gen,spp,assignments)

# Composition plots
ps_sp <- ps_sp %>% subset_samples(Site != 'BLANK') %>% subset_samples(Site != 'Mock')
ps_merge <- ps_sp %>% merge_samples2('Species', fun_otu = mean, funs = list(mean)) # Can be changed to another factor
ps_merge_comp <- ps_merge %>% microbiome::transform('compositional')

ps_plot <- ps_merge_comp %>%
  psmelt() %>%
  ggplot(aes(y=Abundance, x= sample_Species, fill = Class)) + # Change to another factor, or another taxonomic level
  geom_bar(position = 'stack', stat = 'identity')

ps_plot

ps_merge_facet <- ps_merge_comp %>% psmelt() %>%
  filter(Compartment == 'Skeleton' | Compartment == 'Tissue') %>%
  group_by(sample_Species, Compartment, Class) %>% # change to another taxonomic level
  reframe(sum_Abundance = sum(Abundance))

ps_facet <- ps_merge_facet %>%
  ggplot(aes(y=sum_Abundance, x = sample_Species, fill = Class)) +
  geom_bar(position = 'stack', stat = 'identity') +
  facet_wrap(vars(Compartment)) +
  theme_bw()

ps_facet
