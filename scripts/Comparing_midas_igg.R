# comparing IGG and MIDAS

# path should be set to data directory
path <- getwd()
setwd(path)
source("../scripts/Generate_basic_env.R")

##################################
###### Richness comparisons ######
##################################
lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)

library(EcolUtils)
library(readxl)
library(vegan)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(reshape2)

set.seed(999)
lenient <- read.csv("lenient_IGG.txt", sep = "\t", row.names = 8)
igg <- read.csv("default_IGG.txt", sep = "\t", row.names = 8)
midas <- read.csv("default_MIDAS.txt", sep = "\t", row.names = 1)
lenient[,1:7] <- NULL
igg[,1:7] <- NULL
# Metadata
metadata <- read_excel("~/Google Drive File Stream/My Drive/Github/Fiber-Analysis/data/Nutritional_Data.xlsx", skip = 1)
metadata$Intervention <- factor(metadata$Intervention,
                                levels = c('pre','post'),ordered = TRUE)

# check raw data
median(sort(rowSums(t(midas))))
sort(rowSums(t(igg)))

# rarefy
lenient_alpha_rare <- rrarefy.perm(t(lenient), sample = 6700, n = 100, round.out = T)
lenient_rarified <- lenient_alpha_rare[rowSums(lenient_alpha_rare) >= 6700-(6700*.1), 
                                 colSums(lenient_alpha_rare) >= 1]

igg_alpha_rare <- rrarefy.perm(t(igg), sample = 6700, n = 100, round.out = T)
igg_rarified <- igg_alpha_rare[rowSums(igg_alpha_rare) >= 6700-(6700*.1), 
                                       colSums(igg_alpha_rare) >= 1]

midas_alpha_rare <- rrarefy.perm(t(midas), sample = 900, n = 100, round.out = T)
midas_rarified <- midas_alpha_rare[rowSums(midas_alpha_rare) >= 900-(900*.1), 
                               colSums(midas_alpha_rare) >= 1]

# calculate alpha diversity & plot
lenient_richness <- as.data.frame(specnumber(lenient_rarified))
lenient_shannon_div <- diversity(lenient_rarified, "shannon")
lenient_Evenness <- as.data.frame(lenient_shannon_div/log(lenient_richness))

igg_richness <- as.data.frame(specnumber(igg_rarified))
igg_shannon_div <- diversity(igg_rarified, "shannon")
igg_Evenness <- as.data.frame(igg_shannon_div/log(igg_richness))

midas_richness <- as.data.frame(specnumber(midas_rarified))
midas_shannon_div <- diversity(midas_rarified, "shannon")
midas_Evenness <- as.data.frame(midas_shannon_div/log(midas_richness))
midas_all_alpha <- cbind(midas_richness, midas_Evenness)

all_diversity <- cbind(lenient_richness,lenient_Evenness,igg_richness,igg_Evenness)
all_diversity1 <- merge(all_diversity, midas_all_alpha, by.x = "row.names", by.y = "row.names")
colnames(all_diversity1) <- c("sample","lenient_rich","lenient_even","igg_rich","igg_even","midas_rich","midas_even")

all_diversity_meta <- merge(all_diversity1, metadata, by.x = "sample", by.y = "...1")

richness_diversity_meta <- all_diversity_meta %>% select(., 1,2,4,6,15,16) 
richness_diversity_meta_melted <- melt(richness_diversity_meta, id.vars = c("sample", "Intervention", "individual"))
richness_diversity_meta_melted$variable <- gsub(pattern = "lenient_rich", replacement = "IGG-Lenient", x = richness_diversity_meta_melted$variable)
richness_diversity_meta_melted$variable <- gsub(pattern ="midas_rich", replacement = "MIDAS", x = richness_diversity_meta_melted$variable)

evenness_diversity_meta <- all_diversity_meta %>% select(., 1,3,5,7,15,16) 
evenness_diversity_meta_melted <- melt(evenness_diversity_meta, id.vars = c("sample", "Intervention", "individual"))
evenness_diversity_meta_melted$variable <- gsub(pattern = "lenient_even", replacement = "IGG-Lenient", x = evenness_diversity_meta_melted$variable)
evenness_diversity_meta_melted$variable <- gsub(pattern ="midas_even", replacement = "MIDAS", x = evenness_diversity_meta_melted$variable)

richness <- ggplot(data = richness_diversity_meta_melted) + 
  aes(x = Intervention, y = value, fill = variable) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.2) + 
  scale_fill_brewer(palette = "Blues") + labs(x = '', y = 'Number of Species') + 
  theme_bw(base_size = 14) + theme(legend.title = element_blank())

evenness <- ggplot(data = evenness_diversity_meta_melted) + 
  aes(x = Intervention, y = value, fill = variable) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.2) + 
  scale_fill_brewer(palette = "Greens") + labs(x = '', y = 'Pielou J Evenness') + 
  theme_bw(base_size = 14) + theme(legend.title = element_blank())

plot_grid(richness, evenness, ncol = 1, labels = "AUTO")
