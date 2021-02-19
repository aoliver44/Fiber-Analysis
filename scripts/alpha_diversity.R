#################################
######## Alpha Diversity ########
#################################

# path should be set to data directory
path <- getwd()
setwd(path)
source("../scripts/Generate_basic_env.R")

library(vegan)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(nlme)
library(reshape2)
library(ggpubr)
library(janitor)

# calculate alpha diveristy
shannon_div <- diversity(alpha_rare, "shannon") 
richness <- specnumber(alpha_rare)
Evenness <- shannon_div/log(richness)
all_diversity <- cbind(Evenness, richness, shannon_div)

# merge with Nutritional Metadata
alpha_diversity_data <- merge(all_diversity, metadata, by.x = "row.names", by.y = "Metagenome")
alpha_diversity_data_subset <- select(alpha_diversity_data, Individual, Intervention,Cluster, Evenness, richness, shannon_div)

# This is for ggplot. Change the measure vars for what y-axis you want to plot
alpha_diversity_data_subset_melt <- melt(data = alpha_diversity_data_subset, id.vars = c("Individual","Intervention", "Cluster"), measure.vars = c("Evenness", "richness", "shannon_div"))

# removing pseudo-replication by taking mean of the alpha diversity
alpha_diversity_data_means <- alpha_diversity_data_subset_melt %>% group_by(., Individual, Intervention, variable) %>% mutate(., diversity_means = mean(value))

# make single alpha diversity box plots
richness <- ggplot(data = alpha_diversity_data) +
  aes(x = alpha_diversity_data$Intervention, y = alpha_diversity_data$richness, fill = Intervention) +
  geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge(), alpha = 0.3) +
  #geom_jitter(width = 0.15, alpha = 0.2) +
  scale_fill_brewer(palette = "Blues") +
  labs(x = '',
       y = 'Number of Species') +
  theme_bw(base_size = 14) + theme(legend.position = "none")

Evenness <- ggplot(data = alpha_diversity_data) +
  aes(x = alpha_diversity_data$Intervention, y = alpha_diversity_data$Evenness, fill = Intervention) +
  geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge(), alpha = 0.3) +
  #geom_jitter(width = 0.15, alpha = 0.2) +
  scale_fill_brewer(palette = "Greens") +
  labs(x = 'Intervention',
       y = 'Evenness') +
  theme_bw(base_size = 14) + theme(legend.position = "none")

ggplot(data = alpha_diversity_data) +
  aes(x = alpha_diversity_data$Intervention, y = alpha_diversity_data$shannon_div, fill = Intervention) +
  geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge(), alpha = 0.3) +
  #geom_jitter(width = 0.15, alpha = 0.2) +
  scale_fill_brewer(palette = "Blues") +
  labs(x = 'Intervention',
       y = 'Shannon Diversity') +
  theme_bw(base_size = 18) + theme(legend.position = "none")

plot_grid(richness, Evenness, ncol = 1, labels = c("A", "B"))

shapiro.test(alpha_diversity_data$shannon_div)
summary(aov(shannon_div ~ Intervention + Error(as.factor(Individual)), data = alpha_diversity_data))
anova(lme(shannon_div ~ Intervention, random = (~1 | Individual), method = "REML", correlation = corAR1(), data = alpha_diversity_data))
alpha_diversity_data %>%
  group_by(Intervention) %>% 
  summarise_each(funs(mean))

shapiro.test(wilcox_data$value)
wilcox_data <- alpha_diversity_data_means %>% filter(., variable == "shannon_div")
wilcox.test(value ~ Intervention, data = wilcox_data)

########## CAT/BAT ############3

cat_bat_raw <- readRDS(file = "otu_melted_rar14000_filt.rds")
cat_bat_species <- cat_bat_raw %>% select(., 1:101, 107) %>% group_by(., genus) %>% summarise_all(funs(sum)) %>% drop_na() %>%
  column_to_rownames(., var = "genus") %>% t() 
cat_bat_species <- as.data.frame(cat_bat_species) %>% clean_names()

# calculate alpha diveristy
shannon_div <- diversity(cat_bat_species, "shannon") 
richness <- specnumber(cat_bat_species)
Evenness <- shannon_div/log(richness)
all_diversity <- cbind(Evenness, richness, shannon_div)

# merge with Nutritional Metadata
alpha_diversity_data <- merge(all_diversity, metadata, by.x = "row.names", by.y = "Metagenome")
ggplot(data = alpha_diversity_data) +
  aes(x = alpha_diversity_data$Intervention, y = alpha_diversity_data$shannon_div, fill = Intervention) +
  geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge(), alpha = 0.3) +
  #geom_jitter(width = 0.15, alpha = 0.2) +
  scale_fill_brewer(palette = "Blues") +
  labs(x = 'Intervention',
       y = 'Shannon Diversity') +
  theme_bw(base_size = 18) + theme(legend.position = "none")
