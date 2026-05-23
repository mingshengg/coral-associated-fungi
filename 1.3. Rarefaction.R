# Load packages ####
library(vegan); packageVersion('vegan')
library(reshape2); packageVersion('reshape2')
library(tidyr); packageVersion('tidyr')
library(dplyr); packageVersion('dplyr')
library(ggplot2); packageVersion('ggplot2')
library(phyloseq); packageVersion('phyloseq')

ps <- readRDS('./Output/ps_assigned.RDS')
ps_fungi <- ps %>% subset_taxa(Kingdom == 'k__Fungi') %>% subset_samples(Compartment != 'NA')

fungi_otu <- otu_table(ps_fungi) %>% as.data.frame()
fungi_otu$n_seqs <- rowSums(fungi_otu)

### rarefaction curves
rare_curve <- rarecurve(fungi_otu[1:20,], step = 100, label = F, tidy = TRUE)
levels = rare_curve %>% group_by(Site) %>%
  summarise(Max = max(Species)) %>%
  arrange(desc(Max)) %>%
  .$Site

rare_curve$Site <- factor(rare_curve$Site, levels = levels)

ggplot(data = rare_curve, aes(x = Sample, y = Species, colour = Site)) + 
  geom_line(linewidth = 1) +
  theme_bw() +
  theme(legend.position = 'none')


ggsave('./Plots/rarecurve_all.pdf')

# for bottom 20 samples with the least reads
rare_curve20 <- rarecurve(fungi_otu[1:20,], step = 20, label = F, tidy = TRUE)
levels = rare_curve20 %>% group_by(Site) %>%
  summarise(Max = max(Species)) %>%
  arrange(desc(Max)) %>%
  .$Site

rare_curve20$Site <- factor(rare_curve20$Site, levels = levels)

ggplot(data = rare_curve20, aes(x = Sample, y = Species, colour = Site)) + 
  geom_line(linewidth = 1) +
  theme_bw()


ggsave('./Plots/rarecurve_bottom20.pdf')

### visualizing reads distribution
fungi_otu %>% ggplot(aes(x = n_seqs)) +
  geom_histogram(binwidth = 200) +
  coord_cartesian(xlim=c(0,5000))

meta <- microbiome::meta(ps_fungi)
fungi_otu$SampleType = meta$Compartment

fungi_otu %>%
  ggplot(aes(x=1, y=n_seqs, colour = SampleType)) +
  geom_jitter(size = 2) +
  scale_y_log10() +
  theme_bw() +
  labs(y = 'Number of sequences') +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

ggsave('./Plots/number_of_seqs_jitter.pdf')

# Are there drastic differences in reads?
fungi_otu %>%
  arrange(n_seqs) %>%
  ggplot(aes(x=1:nrow(.), y = n_seqs)) +
  geom_line() +
  coord_cartesian(xlim = c(0,500), ylim=c(0,10000))

fungi_otu %>%
  arrange(n_seqs) %>%
  select(n_seqs) # minimum number of reads = 598

# Investigate goods coverage based on different thresholds
fungi_otu$Sample <- rownames(fungi_otu) 
goods_stats <- fungi_otu %>%
  select(-c(n_seqs,SampleType)) %>%
  pivot_longer(-Sample) %>%
  group_by(Sample) %>%
  summarize(n_seqs = sum(value),
            n_sing = sum(value == 2),
            Goods = 1- n_sing/n_seqs) %>%
  filter(n_seqs > 0)

goods_stats$SampleType = fungi_otu$SampleType

ggplot(goods_stats, aes(x = n_seqs, y = Goods, colour = SampleType)) + 
  geom_point(size = 2) +
  coord_cartesian(xlim=c(0,5000), ylim=c(0,1)) +
  coord_cartesian(ylim=c(0.98,1)) +
  theme_bw()

ggsave('./Plots/goods_coverage_doubletons.pdf')

richness_reads <- fungi_otu %>%
  select(-c(n_seqs,SampleType)) %>%
  pivot_longer(-Sample) %>% 
  group_by(Sample) %>% 
  summarize(richness = sum(value > 0),
            n_seqs = sum(value)) 

richness_reads$SampleType <- fungi_otu$SampleType

ggplot(richness_reads, aes(x = n_seqs, y = richness)) + 
  geom_point(aes(colour = SampleType), size = 2) +
  geom_line(stat = 'smooth') +
  theme_bw() +
  coord_cartesian(xlim=c(0,50000)) +
  labs(x='Number of sequences', y='ASV richness')

ggsave('./Plots/richness_w_sequences.pdf')
