# Load packages ####
library(tidyverse); packageVersion('tidyverse')
library(phyloseq); packageVersion('phyloseq')
library(speedyseq); packageVersion('speedyseq')
library(iCAMP); packageVersion('iCAMP')
library(ggplot2); packageVersion('ggplot2')

ps <- readRDS('./Output/ps_assigned.RDS') %>% subset_taxa(Kingdom == 'k__Fungi')
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

meta <- microbiome::meta(ps)

ps_coral <- ps %>% subset_samples(Compartment == 'Skeleton'| Compartment == 'Tissue') %>%
  prune_taxa(taxa_sums(.) > 0, .)

ps_sw <- ps %>% subset_samples(Compartment == 'Seawater') %>%
  prune_taxa(taxa_sums(.) > 0, .)

ps_sed <- ps %>% subset_samples(Compartment == 'Sediment') %>%
  prune_taxa(taxa_sums(.) > 0, .)

#### coral vs sw ####
sloan_coralsw <- snm.comm(ps_coral@otu_table %>% data.frame(), meta.com = ps_sw@otu_table %>% data.frame(),
                           treat = microbiome::meta(ps_coral) %>% select(Compartment))

saveRDS(sloan_coralsw, './Output/Sloan/sloan_coralsw.RDS')
sloan_coralsw <- readRDS('./Output/Sloan/sloan_coralsw.RDS')

sloan_coralsw_plot <- sloan_coralsw$plot.detail %>% filter(treatment.id == 'Tissue') # change here

sloan_coralsw$stats
sloan_coralsw$pvalues

ggplot(data = sloan_coralsw_plot) +
  geom_point(aes(x = log(p), y = freq, colour = type), size = 2) + 
  geom_line(aes(x = log(p), y = pred.lwr), linetype = 'dashed', linewidth = 1.1) +
  geom_line(aes(x = log(p),  y = pred.upr), linetype = 'dashed', linewidth = 1.1) +
  geom_line(aes(x = log(p), y = freq.pred), linewidth = 1.1, colour = 'grey') +
  ggtitle("Tissue vs SW; R2 = -0.3418") + # change accordingly
  theme_bw() +
  labs(y = 'Occurrence frequency', x = 'log(Mean Relative Abundance)') +
  theme(text=element_text(size=14))

sloan_coral_above <- sloan_coralsw$plot.detail %>% filter(treatment.id == 'Tissue') %>% 
  filter(type == 'Above') %>% arrange(desc(freq)) %>% .$OTU
ps_coral@tax_table[sloan_coral_above]

sloan_coral_below <- sloan_coralsw$plot.detail %>% filter(treatment.id == 'Tissue') %>% 
  filter(type == 'Below') %>% arrange(desc(freq)) %>% .$OTU
ps_coral@tax_table[sloan_coral_below]

sloan_coralsw$plot.detail %>% filter(treatment.id == 'Tissue') %>% count(type) %>% mutate(prop = n/sum(n) * 100)

# piechart, plot and later combined 
typecount <- sloan_coralsw_plot %>% count(type)
ggplot(typecount,aes(x='', y=n, fill=type)) +
  geom_col() +
  coord_polar(theta='y')


#### coral vs sed ####
sloan_coralsed <- snm.comm(ps_coral@otu_table %>% data.frame(), meta.com = ps_sed@otu_table %>% data.frame(),
                          treat = microbiome::meta(ps_coral) %>% select(Compartment))

saveRDS(sloan_coralsed, './Output/Sloan/sloan_coralsed.RDS')
sloan_coralsed <- readRDS('./Output/Sloan/sloan_coralsed.RDS')

sloan_coralsed_plot <- sloan_coralsed$plot.detail %>% filter(treatment.id == 'Tissue')

sloan_coralsed$stats
sloan_coralsed$pvalues

ggplot(data = sloan_coralsed_plot) +
  geom_point(aes(x = log(p), y = freq, colour = type), size = 2) + 
  geom_line(aes(x = log(p), y = pred.lwr), linetype = 'dashed', linewidth = 1.1) +
  geom_line(aes(x = log(p),  y = pred.upr), linetype = 'dashed', linewidth = 1.1) +
  geom_line(aes(x = log(p), y = freq.pred), linewidth = 1.1, colour = 'grey') +
  ggtitle("Tissue vs Sed; R2 = 0.0658") +
  theme_bw() +
  labs(y = 'Occurrence frequency', x = 'log(Mean Relative Abundance)') +
  theme(text=element_text(size=14))

sloan_coral_above <- sloan_coralsed$plot.detail %>% filter(treatment.id == 'Tissue') %>% 
  filter(type == 'Above') %>% arrange(desc(freq)) %>% .$OTU
ps_coral@tax_table[sloan_coral_above]

sloan_coral_below <- sloan_coralsed$plot.detail %>% filter(treatment.id == 'Tissue') %>% 
  filter(type == 'Below') %>% arrange(desc(freq)) %>% .$OTU
ps_coral@tax_table[sloan_coral_below]

sloan_coralsed$plot.detail %>% filter(treatment.id == 'Tissue') %>% count(type) %>% mutate(prop = n/sum(n) * 100)

# piechart, plot and later combined 
typecount <- sloan_coralsed_plot %>% count(type)
ggplot(typecount,aes(x='', y=n, fill=type)) +
  geom_col() +
  coord_polar(theta='y')


#### coral vs coral ####
sloan_coralcomp <- snm.comm(ps_coral@otu_table %>% data.frame(), meta.com = ps_coral@otu_table %>% data.frame(),
                           treat = microbiome::meta(ps_coral) %>% select(Compartment))

saveRDS(sloan_coralcomp, './Output/Sloan/sloan_coralcomp.RDS')

sloan_coralcomp <- readRDS('./Output/Sloan/sloan_coralcomp.RDS')
sloan_coralcomp_plot <- sloan_coralcomp$plot.detail %>% filter(treatment.id == 'Skeleton')

sloan_coralcomp$stats
sloan_coralcomp$pvalues

ggplot(data = sloan_coralcomp_plot) +
  geom_point(aes(x = log(p), y = freq, colour = type), size = 2) + 
  geom_line(aes(x = log(p), y = pred.lwr), linetype = 'dashed', linewidth = 1.1) +
  geom_line(aes(x = log(p),  y = pred.upr), linetype = 'dashed', linewidth = 1.1) +
  geom_line(aes(x = log(p), y = freq.pred), linewidth = 1.1, colour = 'grey') +
  ggtitle("Skeleton vs Coral; R2 = 0.7864") +
  theme_bw() +
  labs(y = 'Occurrence frequency', x = 'log(Mean Relative Abundance)') +
  theme(text=element_text(size=14))

sloan_coral_above <- sloan_coralcomp$plot.detail %>% filter(treatment.id == 'Skeleton') %>% 
  filter(type == 'Above') %>% arrange(desc(freq)) %>% .$OTU
ps_coral@tax_table[sloan_coral_above]

sloan_coral_below <- sloan_coralcomp$plot.detail %>% filter(treatment.id == 'Skeleton') %>% 
  filter(type == 'Below') %>% arrange(desc(freq)) %>% .$OTU
ps_coral@tax_table[sloan_coral_below]

# piechart
typecount <- sloan_coralcomp_plot %>% count(type)
ggplot(typecount,aes(x='', y=n, fill=type)) +
  geom_col() +
  coord_polar(theta='y')

sloan_coralcomp$plot.detail %>% filter(treatment.id == 'Skeleton') %>% count(type) %>% mutate(prop = n/sum(n) * 100)

