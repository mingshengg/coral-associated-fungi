# Load packages ####
library(dada2); packageVersion("dada2")
library(purrr); packageVersion("purrr")
library(tidyverse); packageVersion("tidyverse")
library(readxl); packageVersion("readxl")
library(ShortRead); packageVersion("ShortRead")

# Load metadata ####
full_meta <- readxl::read_xlsx("./Metadata.xlsx")

# Find raw fastq files and prepare workspace ####
path <- "./Data"
fqs <- list.files(path, full.names = TRUE, recursive = FALSE, pattern = "_R1.fastq.gz$|_R2.fastq.gz$")

# Parse fwd and rev reads
fnFs <- fqs[grep("_R1.fastq.gz",fqs)]

# Get Sample Names
sample.names <- str_remove(basename(fnFs),"_R1.fastq.gz")

# subset metadata to this run's samples
`%notin%` = Negate(`%in%`)
full_meta %>% filter(SampleID %notin% sample.names) # should be 0, to indicate all samples are accounted for
meta <- full_meta %>% filter(SampleID %in% sample.names)

# Make filtered outfile names
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))

# make new directory for filtered files
if(!dir.exists(file.path(path,"filtered"))){
  dir.create(file.path(path,"filtered"))
}

# check for duplicated sample names
sum(duplicated(sample.names))

# prefilter reads
out <- filterAndTrim(fnFs, filtFs, verbose = TRUE)

saveRDS(out,"./Output/out.RDS") # save object

out <- readRDS("./Output/out.RDS") # reload point

filtpath <- file.path(path,"filtered")

## get sample names again
sample.names <- str_remove(basename(fnFs),"_R1.fastq.gz")

# reassign filts for any potentially lost samples
filtFs <- list.files(filtpath, pattern = "_F_filt", full.names = TRUE)

# Learn error rates ####
errF <- learnErrors(filtFs, multithread=TRUE)
saveRDS(errF,"./Output/errF.RDS") # save object

errF <- readRDS("./Output/errF.RDS") # reload point
out <- readRDS("./Output/out.RDS")

# plot error rates for sanity
plotErrors(errF, nominalQ=TRUE)

# Infer sequence variants ####

# add names to filts
names(filtFs) <- sample.names

# DADA function ####
dadaFs <- dada(filtFs, err = errF, multithread = TRUE)

saveRDS(dadaFs,"./Output/dadaFs.RDS")

# Inspect the sequence variants
head(dadaFs)

# Construct sequence table ####
seqtab <- makeSequenceTable(dadaFs)
saveRDS(seqtab,"./Output/seqtab.RDS") # save object
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

# Remove chimeras ####
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

# Save progress
saveRDS(seqtab.nochim,"./Output/seqtab.nochim.RDS")
seqtab.nochim = readRDS("./Output/seqtab.nochim.RDS")

#Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))

# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "merged", "nonchim")

rownames(track) <- sample.names
head(track)
write.csv(track, './Output/track_reads.csv')

# rarefaction curve
library(vegan)
rare_curve <- rarecurve(seqtab.nochim, step = 50, label = F)

##################################################################################################################
# Assign taxonomy
# export to qiime
uniquesToFasta(seqtab.nochim, fout='./Output/rep-seqs.fna', ids=colnames(seqtab.nochim))

# save as phyloseq object
library(phyloseq)

seqtab.nochim <- readRDS('./Output/seqtab.nochim.RDS')

taxa <- readRDS('./Output/taxa.RDS')

seqtab.df <- as.data.frame(seqtab.nochim)
row.names(seqtab.df)

otu <- otu_table(seqtab.nochim,taxa_are_rows = FALSE)
met <- sample_data(full_meta)
tax <- tax_table(taxa)

sample_names(met) <- met$SampleID

names <- rownames(otu)
rownames(otu) <- names

sample_names(met) <- met$SampleID
taxa_names(tax) <- colnames(otu)

ps <- phyloseq(otu,met,tax)

sequences <- Biostrings::DNAStringSet(taxa_names(ps))
names(sequences) <- taxa_names(ps)
ps <- merge_phyloseq(ps, sequences)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

saveRDS(ps, './Output/ps.RDS')

### decontam ####
library(decontam)

blanks <- which(ps@sam_data$Species == 'BLANK')
contamdf.freq <- isContaminant(ps, neg=blanks)
table(contamdf.freq$contaminant)

ps.noncontam <- prune_taxa(!contamdf.freq$contaminant, ps) #no ASV detected as contaminants

#### Export to BLAST for further taxonomic assignment ####
ps <- readRDS('./Output/ps.RDS')
ps_unassigned <- ps %>% subset_taxa(Kingdom == "Unassigned")
ps_unassigned@refseq %>% Biostrings::writeXStringSet('./Output/unassigned_ASVs.fa')

ps_k_fungi <- ps %>% subset_taxa(Kingdom == 'k__Fungi' & Phylum == 'NA')
ps_k_fungi@refseq %>% Biostrings::writeXStringSet('./Output/fungi_only_ASVs.fa')