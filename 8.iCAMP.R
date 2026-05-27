# Run stegen's unified framework
# Load packages ####
library(iCAMP); packageVersion('iCAMP')
library(phyloseq); packageVersion('phyloseq')
library(dplyr); packageVersion('dplyr')
library(ggplot2); packageVersion('ggplot2')

# Read data
ps_tree <- readRDS('./Output/ps_fungi_rare_tree.RDS') %>% subset_taxa(Kingdom == 'k__Fungi')
ps_tree_comp <- microbiome::transform(ps_tree, transform = 'compositional')

# ps tree sanity checks #
length(phy_tree(ps_tree)$tip.label)
all(sort(taxa_names(ps_tree)) == sort(phy_tree(ps_tree)$tip.label))
summary(phy_tree(ps_tree)$edge.length)
plot_tree(ps_tree, color = "Phylum")

#### Coral Tissues #### 
ps_tissue <- ps_tree_comp %>% subset_samples(Compartment == 'Tissue')
ps_tissue <- prune_taxa(taxa_sums(ps_tissue) > 0, ps_tissue)

bnti_icamp_tissue <- bNTIn.p(ps_tissue@otu_table, cophenetic(ps_tissue@phy_tree),
                               weighted = TRUE) # conduct both weighted and unweighted analysis

saveRDS(bnti_icamp_tissue, './Output/bnti_icamp_tissue_weighted.RDS')
bnti_icamp <- readRDS('./Output/bnti_icamp_tissue_weighted.RDS')

bnti_icamp_tissue_val <- bnti_icamp_tissue$index[upper.tri(bnti_icamp_tissue$index)]
total_comparisons <- length(bnti_icamp_tissue_val) # no of comparisons

# Calculate proportions
variable_selection <- sum(bnti_icamp_tissue_val > 1.96) / total_comparisons * 100
homogeneous_selection <- sum(bnti_icamp_tissue_val < -1.96) / total_comparisons * 100
stochastic_processes <- sum(bnti_icamp_tissue_val >= -1.96 & bnti_icamp_tissue_val <= 1.96) / total_comparisons * 100

# Print results
cat("Variable selection:", variable_selection, "%\n")
cat("Homogeneous selection:", homogeneous_selection, "%\n")
cat("Stochastic processes:", stochastic_processes, "%\n")

rcbray_icamp_tissue <- RC.pc(ps_tissue@otu_table %>% t(),
                               weighted = TRUE)

saveRDS(rcbray_icamp_tissue, './Output/rcbray_tissue_weighted.RDS')

stochastic_pairs <- which(bnti_icamp_tissue_val >= -2 & bnti_icamp_tissue_val <= 2)
rc_bc_stochastic <- rcbray_icamp_tissue$index[upper.tri(rcbray_icamp_tissue$index)][stochastic_pairs]

# Classify dispersal processes
dispersal_limitation <- sum(rc_bc_stochastic > 0.95) / length(rc_bc_stochastic) * 100
homogenizing_dispersal <- sum(rc_bc_stochastic < -0.95) / length(rc_bc_stochastic) * 100
ecological_drift <- sum(rc_bc_stochastic >= -0.95 & rc_bc_stochastic <= 0.95) / length(rc_bc_stochastic) * 100

# Print results
cat("Dispersal limitation:", dispersal_limitation, "%\n")
cat("Homogenizing dispersal:", homogenizing_dispersal, "%\n")
cat("Ecological drift:", ecological_drift, "%\n")

#### Coral skeletons ####
ps_skeleton <- ps_tree %>% subset_samples(Compartment == 'Skeleton')

bnti_icamp_skeleton <- bNTIn.p(ps_skeleton@otu_table, cophenetic(ps_skeleton@phy_tree),
                             weighted = FALSE) # conduct both weighted and unweighted analysis

saveRDS(bnti_icamp, './Output/bnti_icamp_skeleton_unweighted.RDS')
bnti_icamp <- readRDS('./Output/bnti_icamp_skeleton_unweighted.RDS')

bnti_icamp_skeleton_val <- bnti_icamp_skeleton$index[upper.tri(bnti_icamp_skeleton$index)]
total_comparisons <- length(bnti_icamp_skeleton_val) # no of comparisons

# Calculate proportions
variable_selection <- sum(bnti_icamp_skeleton_val > 1.96) / total_comparisons * 100
homogeneous_selection <- sum(bnti_icamp_skeleton_val < -1.96) / total_comparisons * 100
stochastic_processes <- sum(bnti_icamp_skeleton_val >= -1.96 & bnti_icamp_skeleton_val <= 1.96) / total_comparisons * 100

# Print results
cat("Variable selection:", variable_selection, "%\n")
cat("Homogeneous selection:", homogeneous_selection, "%\n")
cat("Stochastic processes:", stochastic_processes, "%\n")

rcbray_icamp_skeleton <- RC.pc(ps_skeleton@otu_table,
                             weighted = FALSE)

saveRDS(rcbray_icamp_skeleton, './Output/rcbray_skeleton.RDS')

stochastic_pairs <- which(bnti_icamp_skeleton_val >= -2 & bnti_icamp_skeleton_val <= 2)
rc_bc_stochastic <- rcbray_icamp_skeleton$index[upper.tri(rcbray_icamp_skeleton$index)][stochastic_pairs]

# Classify dispersal processes
dispersal_limitation <- sum(rc_bc_stochastic > 0.95) / length(rc_bc_stochastic) * 100
homogenizing_dispersal <- sum(rc_bc_stochastic < -0.95) / length(rc_bc_stochastic) * 100
ecological_drift <- sum(rc_bc_stochastic >= -0.95 & rc_bc_stochastic <= 0.95) / length(rc_bc_stochastic) * 100

# Print results
cat("Dispersal limitation:", dispersal_limitation, "%\n")
cat("Homogenizing dispersal:", homogenizing_dispersal, "%\n")
cat("Ecological drift:", ecological_drift, "%\n")


## Plot venn diagrams
icamp_val <- readxl::read_xlsx('./Output/DADA2/icamp.xlsx')

icamp_val_plot <- icamp_val %>% filter(Compartment == 'Skeleton' & Weighted == 'Unweighted')

ggplot(icamp_val_plot,aes(x='', y=Percentage, fill=Type)) +
  geom_col() +
  coord_polar(theta='y')
