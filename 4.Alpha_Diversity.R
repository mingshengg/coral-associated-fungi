# Load packages, data, and customizations ####
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(vegan); packageVersion("vegan")
library(patchwork); packageVersion("patchwork")
library(microbiome); packageVersion("microbiome")
library(broom); packageVersion("broom")
library(lme4); packageVersion("lme4")
library(lmerTest); packageVersion("lmerTest")
library(emmeans); packageVersion("emmeans")
library(metagMisc); packageVersion("metagMisc")

se <- function(x) sd(x)/sqrt(length(x))

# Load ps objects
ps <- readRDS('./Output/ps_fungi_rare.RDS') %>%
  subset_samples(Individual !='NA') %>%
  subset_taxa(Kingdom == 'k__Fungi')
ps <-  prune_taxa(taxa_sums(ps) > 0, ps)
fungi_tax <- ps@tax_table %>% data.frame()

# ps object only with corals
ps_coral <- ps %>% subset_samples(Compartment == 'Skeleton' | Compartment == 'Tissue')
ps_coral <- prune_taxa(taxa_sums(ps_coral) > 0, ps_coral)
ps_coral@tax_table %>% as.data.frame() %>% count(Phylum)

# Model alpha diversity ####
meta <- microbiome::meta(ps)                
meta$Shannon <- vegan::diversity(otu_table(ps),index = "shannon")
meta$Richness <- vegan::specnumber(otu_table(ps))
meta$Evenness <- microbiome::evenness(otu_table(ps), index = "simpson") %>% .$simpson

# add to ps object
ps@sam_data$Richness <- meta$Richness
ps@sam_data$Shannon <- meta$Shannon
ps@sam_data$Evenness <- meta$Evenness

# Explore shannon, richness, and evenness with violin plots
shannon_plot <- meta %>% 
  ggplot(aes(x = Species, y = Shannon, fill = Compartment)) + 
  geom_boxplot() +
  #geom_violin(scale = 'width', trim = FALSE) +
  geom_jitter(width = 0.3, alpha = 0.3) + 
  theme_classic() +
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title = element_text(size=20),
        legend.position = "none")

shannon_plot

ggsave("./Plots/shannon.pdf",dpi=400,width = 12,height = 6)

richness_plot <- meta %>% 
  ggplot(aes(x = Species, y = Richness, fill = Compartment)) + 
  geom_boxplot() +
  geom_jitter(alpha = 0.3) + 
  theme_classic() +
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title = element_text(size=20),
        legend.position = "none")

richness_plot

ggsave("./Plots/richness.pdf",dpi=400,width = 12,height = 6)

evenness_plot <- meta %>%
  ggplot(aes(x = Species, y = Evenness, fill = Compartment)) + 
  geom_boxplot() +
  geom_jitter(alpha = 0.3) +
  ylim(c(0,1.0)) +
  theme_classic() +
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title = element_text(size= 20),
        legend.position = "none")

evenness_plot

ggsave("./Plots/evenness.pdf",dpi=400,width = 12,height = 6)
  
# Overall ANOVA model
shannon_overall <- aov((Shannon) ~ Compartment*Site, data = meta)
summary(shannon_overall)
plot(shannon_overall)
shapiro.test(resid(shannon_overall))
car::leveneTest((Shannon) ~ Compartment*Site, data = meta) # homogeneous variance
emmeans(shannon_overall, pairwise ~ Compartment)$contrasts %>% write.table('./Output/shannon-paired.txt', sep='\t')

richness_overall <- aov(log(Richness) ~ Compartment*Site, data = meta)
summary(richness_overall)
plot(richness_overall)
shapiro.test(resid(richness_overall))
car::leveneTest(log(Richness) ~ Compartment*Site, data = meta) 
emmeans(richness_overall, pairwise ~ Compartment)$contrasts %>% write.table('./Output/richness-paired.txt', sep='\t')
emmeans(richness_overall, pairwise ~ Site)
meta %>% group_by(Site) %>% summarise(Mean = mean(Richness), Se = se(Richness))

evenness_overall <- aov(log(Evenness) ~ Compartment+Site, data = meta)
summary(evenness_overall)
plot(evenness_overall)
shapiro.test(resid(evenness_overall))
car::leveneTest(log(Evenness) ~ Compartment*Site, data = meta) 
emmeans(evenness_overall, pairwise ~ Compartment)$contrasts %>% write.table('./Output/evenness-paired.txt', sep='\t')
emmeans(evenness_overall, pairwise ~ Site)
meta %>% group_by(Site) %>% summarise(Mean = mean(Evenness), Se = se(Evenness))

# average across compartments
meta %>% group_by(Compartment) %>%
  summarise(Shannon_mean = mean(Shannon), Shannon_se = se(Shannon),
            Richness_mean = mean(Richness), Richness_se = se(Richness),
            Evenness_mean = mean(Evenness), Evenness_se = se(Evenness))

meta %>% filter(Compartment == 'Tissue' | Compartment == 'Skeleton') %>%
  summarise(Shannon_mean = mean(Shannon), Shannon_se = se(Shannon),
            Richness_mean = mean(Richness), Richness_se = se(Richness),
            Evenness_mean = mean(Evenness), Evenness_se = se(Evenness))

meta %>% filter(Compartment == 'Sediment' | Compartment == 'Seawater') %>%
  summarise(Shannon_mean = mean(Shannon), Shannon_se = se(Shannon),
            Richness_mean = mean(Richness), Richness_se = se(Richness),
            Evenness_mean = mean(Evenness), Evenness_se = se(Evenness))

# Separating just coral data
coral_meta <- meta %>% filter(Compartment == 'Tissue'| Compartment == 'Skeleton')

shannon_coral <- aov(Shannon ~ Species*Compartment + Site, data = coral_meta)
summary(shannon_coral)

plot(shannon_coral)
shapiro.test(resid(shannon_coral))

emmeans(shannon_coral, pairwise ~ Species|Compartment)
emmeans(shannon_coral, pairwise ~ Compartment|Species)

richness_coral <- aov(log(Richness) ~ Species*Compartment + Site, data = coral_meta)
summary(richness_coral)

plot(richness_coral)
shapiro.test(resid(richness_coral))

emmeans(richness_coral, pairwise ~ Species|Compartment)
emmeans(richness_coral, pairwise ~ Compartment|Species)

evenness_coral <- aov(Evenness ~ Species*Compartment + Site, data = coral_meta)
summary(evenness_coral)

plot(evenness_coral)
shapiro.test(resid(evenness_coral))