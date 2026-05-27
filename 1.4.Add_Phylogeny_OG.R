# Ran on High-performance computing cluster, National University of Singapore
setwd('/home/svu/e0321475/Chapter2/')

library(tidyverse)
library(Biostrings)
library(phyloseq)
library(DECIPHER)
library(phangorn)

# Grab sequences as a DNAStringSet Object
sequences <- readRDS('./Output/ps_assigned.fa')

# Multiple sequence alignment
alignment <- AlignSeqs(sequences, 
                       refinements = 3,
                       iterations = 2,
                       gapOpening = c(-18, -16),
                       gapExtension = c(-2, -1),
                       processors = NULL)

saveRDS(alignment, './Output/alignment_decipher.RDS')

alignment <- readRDS('./Output/alignment_decipher.RDS')

alignment_staggered <- StaggerAlignment(alignment)
saveRDS(alignment_staggered, './Output/alignment_staggered_decipher.RDS') 

# Convert to phangorn format
phang.align = phyDat(as(alignment_staggered, 'matrix'), type = 'DNA')
write.phyDat(phang.align, './Output/alignment_staggered_decipher.nex', format = 'nexus')

# Distance max likelihood
dm <- dist.ml(phang.align)
saveRDS(dm, './Output/ML_Distance.RDS')

# Initial neighbour-joining tree
treeNJ <- NJ(dm)
saveRDS(treeNJ, './Output/treeNJ.RDS')

# Estimate model parameters
fit <- pml(treeNJ, data = phang.align)
saveRDS(fit, './Output/fit_treeNJ.RDS')

# Likelihood of tree
fitJC <- optim.pml(fit, TRUE)
saveRDS(fitJC, './Output/tree_fitJC.RDS')

fitJC$tree$tip.label <- names(sequences)
saveRDS(fitJC, './Output/tree_fitJC_named.RDS')

identical(fitJC$tree$tip.label, taxa_names(ps))

# add tree to phyloseq object ####
ps2 <- phyloseq(tax_table(tax_table(ps)),
                otu_table(otu_table(ps)),
                sample_data(sample_data(ps)),
                phy_tree(fitJC$tree))

saveRDS('./Output/ps_assigned_tree.RDS')
