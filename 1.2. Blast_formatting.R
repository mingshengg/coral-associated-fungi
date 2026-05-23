# Load packages ####
library(dplyr); packageVersion('dplyr')

### taxonomy update from unassigned ASVs
taxa <- readxl::read_xlsx('./Output/BLAST results/processed/fungi_only_ASVs.xlsx') %>% as.data.frame()
taxa <- readxl::read_xlsx('./Output/BLAST results/processed/unassigned_ASVs.xlsx') %>% as.data.frame() #repeat twice for unassigned and those only with kingdom fungi assignment

taxa_dat <- taxa[,2]
n_tab <- length(taxa_dat)
taxa_table <- matrix(, nrow = n_tab, ncol = 8)

for (i in 1:n_tab){
  temp = strsplit(taxa_dat[i], ';')
  taxa_table[i,] = c(temp[[1]], rep('NA', 8 - length(temp[[1]])))
}

rownames(taxa_table) <- taxa[,1]

taxa_table

ps <- readRDS('./Output/ps_assigned.RDS') #first time should be ps.RDS
tax <- ps@tax_table %>% as.data.frame()

tax %>% filter(rownames(tax) %in% rownames(taxa_table))

taxa_table -> tax[match(rownames(taxa_table), rownames(tax)), ]

## checking to ensure that taxonomy has been updated
tax %>% filter(rownames(tax) %in% rownames(taxa_table))

ps_new <- ps
ps_new@tax_table <- tax_table(as.matrix(tax))
temp <- ps_new@tax_table %>% as.data.frame()

saveRDS(ps_new, './Output/ps_assigned.RDS')

# final track reads
ps_assigned <- readRDS('./Output/ps_assigned.RDS') %>% subset_taxa(Kingdom == 'k__Fungi')
track_reads <- read.csv('./Output/track_reads_version2.csv')

fung_reads <- otu_table(ps_assigned) %>% rowSums()
track_reads$X == names(fung_reads)

track_reads$Fungal_reads <- fung_reads

write.csv(track_reads, './Output/track_reads_final.csv')
