# Load packages ####
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(vegan); packageVersion("vegan")
library(VennDiagram); packageVersion("VennDiagram")
library(patchwork); packageVersion("patchwork")
library(ggplot2); packageVersion("ggplot2")
library(speedyseq); packageVersion("speedyseq")

`%notin%` <- Negate(`%in%`)
se <- function(x){
  sd(x)/sqrt(length(x))
}

# custom palette
pal <- c("#d98416","#25802d","#664c13")
# names(pal) <- c("f","l","p","s")

# Load ps object ####
full_ps <- readRDS("./Output/DADA2/ps_fungi_rare.RDS") %>% subset_taxa(Kingdom == 'k__Fungi')

# agglomerate taxa ####
full_ps <- microbiome::transform(full_ps, 'compositional')
full_ps <- full_ps %>% subset_samples(Temperature != "NA")
full_ps <- full_ps %>% subset_taxa(taxa_sums(full_ps)>0)

# Phyla distribution
full_ps_phylum <- full_ps %>%
  subset_samples(Compartment == 'Tissue') %>% # change for other compartments
  tax_glom('Genus')

otu <- as.data.frame(otu_table(full_ps_phylum))
colnames(otu) <- as.data.frame(tax_table(full_ps_phylum))$Genus
otu %>% apply(2, function(x) mean(x)*100) %>% round(2) %>% order()
otu %>% apply(2, function(x) se(x)*100) %>% round(2)

# OTUs without order
full_ps@tax_table %>% as.data.frame() %>% count(Genus)

full_ps_phylum <- prune_taxa(taxa_sums(full_ps_phylum) > 0, full_ps_phylum)
full_ps_phylum@tax_table %>% as.data.frame() %>% count(Genus)

# Most common OTUs (relative abundance)
otu_tab <- otu_table(full_ps_phylum %>% subset_samples(Compartment == 'Tissue')) 
top10 <- otu_tab %>% colMeans() %>% sort(decreasing = T) %>% .[1:10] %>% names()
otu_tab %>% as.data.frame() %>% .[top10] %>% apply(2, se)
full_ps_phylum@tax_table[top10] %>% as.data.frame()

# Most occurring OTUs (found in highest number of samples)
otu_tab <- otu_table(full_ps %>% subset_samples(Species == 'Diploastrea_heliopora' & Compartment == 'Skeleton'))
otu_tab_o <- (otu_tab > 0) * 1
tax_occ <- (colSums(otu_tab_o)/nrow(otu_tab_o)) %>% sort(decreasing = TRUE)
top <- tax_occ %>% .[1:5] %>% names()
full_ps@tax_table[top] %>% as.data.frame() %>% select(Species)

Biostrings::writeXStringSet(full_ps@refseq['ASV16'], './Output/ASVs/ASV16.fa')

# VennDiagram of overlap for compartments ####
ps_venn <- full_ps #Different taxonomic levels
A <- ps_venn %>% subset_samples(Compartment == "Skeleton") %>% 
  prune_taxa(taxa_sums(.)>0, .) %>% tax_table() %>% row.names()
B <- ps_venn %>% subset_samples(Compartment == "Tissue") %>% 
  prune_taxa(taxa_sums(.)>0, .) %>% tax_table() %>% row.names()
C <- ps_venn %>% subset_samples(Compartment == "Seawater") %>% 
  prune_taxa(taxa_sums(.)>0, .) %>% tax_table() %>% row.names()
D <- ps_venn %>% subset_samples(Compartment == "Sediment") %>% 
  prune_taxa(taxa_sums(.)>0, .) %>% tax_table() %>% row.names()

VennDiagram::venn.diagram(list(A,B,C,D),
                          category.names = c("Skeleton","Tissue","Seawater","Sediment"),
                          filename = './Plots/compartment_venn_OTU.tiff',
                          output=FALSE)


A <- ps_venn %>% subset_samples(Species == "Diploastrea_heliopora") %>% 
  subset_taxa(taxa_sums(full_ps_genus %>% subset_samples(Species == "Diploastrea_heliopora"))>0) %>% tax_table() %>% row.names()
B <- ps_venn %>% subset_samples(Species == "Pachyseris_speciosa") %>% 
  subset_taxa(taxa_sums(full_ps_genus %>% subset_samples(Species == "Pachyseris_speciosa"))>0) %>% tax_table() %>% row.names()
C <- ps_venn %>% subset_samples(Species == "Pocillopora_acuta") %>% 
  subset_taxa(taxa_sums(full_ps_genus %>% subset_samples(Species == "Pocillopora_acuta"))>0) %>% tax_table() %>% row.names()

VennDiagram::venn.diagram(list(A,B,C),
                          category.names = c("Diploastrea_heliopora","Pachyseris_speciosa","Pocillopora_acuta"),
                          filename = './Plots/coral_species_venn_OTU.tiff',
                          output=TRUE)

dev.off()

A <- full_ps_genus %>% subset_samples(Compartment == "Tissue") %>% 
  subset_taxa(taxa_sums(full_ps_genus %>% subset_samples(Compartment == "Tissue"))>0) %>% tax_table() %>% row.names()
B <- full_ps_genus %>% subset_samples(Compartment == "Skeleton") %>% 
  subset_taxa(taxa_sums(full_ps_genus %>% subset_samples(Compartment == "Skeleton"))>0) %>% tax_table() %>% row.names()
C <- full_ps_genus %>% subset_samples(Compartment == "Sediment") %>% 
  subset_taxa(taxa_sums(full_ps_genus %>% subset_samples(Compartment == "Sediment"))>0) %>% tax_table() %>% row.names()
D <- full_ps_genus %>% subset_samples(Compartment == "Seawater") %>% 
  subset_taxa(taxa_sums(full_ps_genus %>% subset_samples(Compartment == "Seawater"))>0) %>% tax_table() %>% row.names()

# for triple venn; dev.off() to close previous venn diagram
full <- unique(c(A,B,C))
n12 <- sum(full %in% unique(A) & full %in% unique(B))
n13 <- sum(full %in% unique(A) & full %in% unique(C))
n23 <- sum(full %in% unique(B) & full %in% unique(C))

n123 <- sum(full %in% unique(A) & full %in% unique(B) & full %in% unique(C))

venn.plot1 <- VennDiagram::draw.triple.venn(length((A)), length((B)), length((C)),
                                            n12, n23, n13,
                                            n123,
                                            category = c('Tissue', 'Seawater', 'Sediment'))

## for quad venn
full <- unique(c(A,B,C,D))
n12 <- sum(full %in% unique(A) & full %in% unique(B))
n13 <- sum(full %in% unique(A) & full %in% unique(C))
n14 <- sum(full %in% unique(A) & full %in% unique(D))
n23 <- sum(full %in% unique(B) & full %in% unique(C))
n24 <- sum(full %in% unique(B) & full %in% unique(D))
n34 <- sum(full %in% unique(C) & full %in% unique(D))

n123 <- sum(full %in% unique(A) & full %in% unique(B) & full %in% unique(C))
n124 <- sum(full %in% unique(A) & full %in% unique(B) & full %in% unique(D))
n134 <- sum(full %in% unique(A) & full %in% unique(C) & full %in% unique(D))
n234 <- sum(full %in% unique(B) & full %in% unique(C) & full %in% unique(D))

n1234 <- sum(full %in% unique(A) & full %in% unique(B) & full %in% unique(C) & full %in% unique(D))

venn.plot1 <- VennDiagram::draw.quad.venn(length((A)),length((B)),length((C)), length((D)),
                                          n12,n13,n14,n23,n24,n34,
                                          n123,n124,n134,n234,n1234,
                                          category = c('Tissue','Seawater',
                                                       'Sediment', 'Skeleton'))

## for dual venn
n12 <- sum(full %in% unique(A) & full %in% unique(B))
VennDiagram::draw.pairwise.venn(length((A)), length((B)),
                                n12,
                                category = c("Tissue", "Skeleton"))

# common in all four compartments
common_four <- full[full %in% unique(A) & full %in% unique(B) & full %in% unique(C) & full %in% unique(D)]
full_ps_genus@tax_table[common_four]

full_ps_genus@otu_table[,common_four] %>% as.data.frame() %>% rowSums()

# Exclusively found
# Tissue
exclusive_tissue <- setdiff(B, union(A, union(C,D)))
full_ps_genus@tax_table[exclusive_tissue]
top_excl_tissue <- full_ps %>% subset_samples(Compartment == "Skeleton") %>% 
  otu_table() %>% .[exclusive_tissue,] %>% rowMeans() %>% 
  sort(decreasing = TRUE) %>% .[1:20] %>% names()

full_ps@tax_table[top_excl_tissue] %>% as.data.frame()

# Compositional plots
full_ps_genus <- full_ps %>% subset_samples(Species != 'BLANK' & Species != 'Mock')
meta <- sample_data(full_ps_genus)
meta_new <- meta %>% data.frame() %>% select(-'Sequenced')
full_ps_genus <- phyloseq(sample_data(meta_new), otu_table(full_ps_genus), tax_table(full_ps_genus))

ps_merge <- full_ps_genus %>% merge_samples2('Species', fun_otu = mean, funs = list(mean))
ps_merge_comp <- ps_merge %>% microbiome::transform('compositional')

tax_tab <- data.frame(ps_merge_comp@tax_table)
tax_tab$Class <- ifelse(str_detect(tax_tab$Class, 'Incertae_sedis'), 'NA', tax_tab$Class)
ps_merge_comp@tax_table <- tax_table(as.matrix(tax_tab))

compositional_plot <- ps_merge_comp %>%
  psmelt() %>% 
  ggplot(aes(y = Abundance, x = Sample, fill = Phylum)) + #change fill = Taxonomic rank
  geom_bar(stat = 'identity') +
  theme_minimal()

compositional_plot

ps_melt <- full_ps_genus %>%
  subset_samples(Compartment == 'Tissue' | Compartment == 'Skeleton') %>%
  microbiome::transform('compositional') %>%
  psmelt()

# select top 10 Class with the highest relative abundance
ps_melt$Class <- ifelse(str_detect(ps_melt$Class, 'Incertae_sedis'), 'NA', ps_melt$Class)
class_names <- ps_melt %>% group_by(Class) %>% summarise(mean_rel = mean(Abundance)) %>% 
  arrange(desc(mean_rel)) %>% .[1:13,] %>% pull(Class)

ps_melt$Class10 <- ifelse(ps_melt$Class %in% class_names, ps_melt$Class, 'Others')

ps_melt %>% ggplot(aes(y = Abundance, x = Sample, fill = Class10)) +
  geom_bar(stat = 'identity') +
  theme_minimal() +
  facet_wrap(sample_Species ~ Compartment, scales ="free_x", ncol = 3, nrow = 3) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank())

# select top 10 Phylum with the highest relative abundance
ps_melt$Phylum <- ifelse(str_detect(ps_melt$Phylum, 'Incertae_sedis'), 'NA', ps_melt$Phylum)

ps_melt %>% ggplot(aes(y = Abundance, x = Sample, fill = Phylum)) +
  geom_bar(stat = 'identity') +
  theme_minimal() +
  facet_wrap(sample_Species ~ Compartment, scales ="free_x", ncol = 3, nrow = 3) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank())


# investigating relative abundances
ps_melt_coral <- full_ps_genus %>% 
  subset_samples(Compartment == 'Seawater') %>% # change to other sample types
  microbiome::transform('compositional') %>%
  psmelt() 

ps_melt_coral %>% 
  group_by(Phylum) %>% 
  summarise(mean_rel = sum(Abundance)/9 * 100, se = se(Abundance)/126*100)
