# Differential abundance using ANCOM-BC2 ####
# Load packages ####
library(ANCOMBC); packageVersion('ANCOMBC')
library(ggrepel); packageVersion('ggrepel')
library(ggplot2); packageVersion('ggplot2')
library(dplyr); packageVersion('dplyr')
library(phyloseq); packageVersion('phyloseq')

# ANCOM-BC2 takes count data
ps <- readRDS('./Output/DADA2/ps_assigned.RDS') %>% 
  subset_taxa(Kingdom == 'k__Fungi')

# Reformatting ps-object
ps_da <- ps %>% subset_samples(Compartment == 'Skeleton' | Compartment == 'Sediment') %>%
  prune_taxa(taxa_sums(.) > 0, .)

tax <- as.data.frame(tax_table(ps_da))

ancom_result <- ancombc2(data=ps_da,
                                fix_formula = "Compartment",
                                p_adj_method = "BH",
                                prv_cut = 0, #prevalence filter
                                lib_cut = 0, #filtering samples based on library sizes
                                group = "Compartment", #for detection of structural zeros and performing global tests
                                struc_zero = FALSE, #to detect structural zeros
                                neg_lb = FALSE) #whether to classify taxon as a structural zero using its asymptotic lower bound, recommended for n>30


saveRDS(ancom_result, './Output/ancom_result_skeleton-sed.RDS')

ancom_result <- readRDS('./Output/ancom_result_compartment.RDS')

res_df <- ancom_result$res

# volcano plot
ggplot(res_df, aes(x = `lfc_CompartmentTissue`, y = -log10(`q_CompartmentTissue`), color = `diff_CompartmentTissue`)) +
  geom_point() +
  scale_color_manual(values = c("grey","red")) +
  theme_minimal() +
  labs(title = "Skeleton vs Tissue",
       x = "Log2 Fold Change",
       y = '-log10 Adjusted p-value') +
  geom_text_repel(data=subset(res_df, `q_CompartmentTissue` < 0.0001),
                  aes(label = taxon),
                  size = 3, max.overlaps = 20)


# What is different?
tissue_taxa <- ps@tax_table[res_df %>% filter(diff_CompartmentTissue == TRUE & lfc_CompartmentTissue > 0) %>% pull(taxon)] %>% 
  data.frame() 

skeleton_taxa <- ps_da@tax_table[res_df %>% filter(diff_CompartmentTissue == TRUE & lfc_CompartmentTissue < 0) %>% pull(taxon)] %>% 
  data.frame()

# Compare with Sloan Neutral Model
sloan_res <- readRDS('./Output/Sloan/sloan_tissue.RDS')
ancom_res <- readRDS('./Output/ancom_result_tissue-sw.RDS')

sloan <- sloan_res$plot.detail %>% filter(treatment.id == 'Tissue')
sloan$taxon <- sloan$OTU

ancomres <- ancom_res$res
ancomres <- ancomres %>%
  left_join(sloan %>% select(taxon, type), by = 'taxon')

# volcano plot
ggplot(subset(ancomres, `q_CompartmentTissue` < 0.1) , aes(x = `lfc_CompartmentTissue`, y = -log10(`q_CompartmentTissue`), color = `type`, shape = `diff_CompartmentTissue`)) +
  geom_point() +
  theme_minimal() +
  labs(title = "Tissue vs Seawater",
       x = "Log2 Fold Change",
       y = '-log10 Adjusted p-value') +
  geom_text_repel(data=subset(ancomres, `q_CompartmentTissue` < 0.0001),
                  aes(label = taxon),
                  size = 3, max.overlaps = 20)

tax_list <- ancomres %>% filter(diff_CompartmentTissue == TRUE & lfc_CompartmentTissue > 0) %>%
  pull(taxon)
