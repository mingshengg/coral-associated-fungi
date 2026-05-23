# Load packages ####
library(dplyr); packageVersion('dplyr')

# Output from QIIME2
taxa <- read.csv('./Output/taxonomy.tsv', sep ='\t')

# QIIME2 taxonomy output file: 1st col: DNA sequences; 2nd col: assigned taxonomy; 3rd col: confidence
taxa_dat <- taxa[,2]
n_tab <- length(taxa_dat)
taxa_table <- matrix(, nrow = n_tab, ncol = 8)

for (i in 1:n_tab){
  temp = strsplit(taxa_dat[i], ';')
  taxa_table[i,] = c(temp[[1]], rep('NA', 8 - length(temp[[1]])))
}

colnames(taxa_table) <- c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Strain')
rownames(taxa_table) <- taxa$Feature.ID

saveRDS(taxa_table, './Output/taxa.RDS')

### taxonomy update from unassigned ASVs, output from BLAST outfmt 6
taxa <- readxl::read_xlsx('./Output/BLAST results/unassigned_ASVs.xlsx') %>% as.data.frame()
taxa_dat <- taxa[,2]
n_tab <- length(taxa_dat)
taxa_table <- matrix(, nrow = n_tab, ncol = 8)

for (i in 1:n_tab){
  temp = strsplit(taxa_dat[i], ';')
  taxa_table[i,] = c(temp[[1]], rep('NA', 8 - length(temp[[1]])))
}

rownames(taxa_table) <- taxa[,1]

ps <- readRDS('./Output/ps_assigned.RDS')
tax <- ps@tax_table %>% as.data.frame()

tax %>% filter(rownames(tax) %in% rownames(taxa_table))

# update taxonomy table
taxa_table -> tax[match(rownames(taxa_table), rownames(tax)), ]

## checking to ensure that taxonomy has been updated
tax %>% filter(rownames(tax) %in% rownames(taxa_table))

ps_new <- ps
ps_new@tax_table <- tax_table(as.matrix(tax))
temp <- ps_new@tax_table %>% as.data.frame()

saveRDS(ps_new, './Output/ps_assigned.RDS')
