#################################
######### Taxa Barplots #########
#################################

source("~/Google Drive File Stream/My Drive/Github/Fiber-Analysis/scripts/Generate_basic_env.R")

library(reshape2)
library(ggplot2)
library(tidyverse)

# now do a little stuff outside of R, for the final taxonomy stuff
out_melted <- read.csv("~/Google Drive File Stream/My Drive/Github/Fiber-Analysis/data/taxonomy_melted_rel_abund.tsv", sep = "\t", header = T)
Genus_melted <- out_melted %>% filter(., L4 != "") %>% group_by(L4, variable) %>% summarise(genus_sum = sum(value))
taxonomy_melt_meta <- merge(Genus_melted, metadata, by.x = "variable", by.y = "Metagenome")

cb7 <- c("#f842b5", "#856800", "#8c159c", "#ff8c48", "#0057b0", "#ff4566", "#7289bd")

cb27 <- c("#dd7f60", "#746dd8", "#7daf3d", "#c864c5", "#6cb558", "#4f388a", "#b9b337", "#5083de", "#d38a32", "#668cd1", "#c2562e", "#43c8ac", "#ca427f", "#5bc47e", "#852c78", "#b7b457", "#b682d3", "#3d7831", "#d975b5", "#777425", "#892c5a", "#c2914b", "#94273d",
          "#de6b87", "#87361a", "#d5474f", "#d06163")

ggplot(data = taxonomy_melt_meta, aes(x = Day, y = genus_sum, fill = L4)) +
  geom_area() +
  theme_bw() + facet_grid(. ~ Individual, space = "free", scales = "free") + 
  scale_fill_manual(values = cb27) +
  theme(panel.spacing = unit(0.1, "lines")) +   
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# L6 genus level
L6 <- read.csv("~/Google Drive File Stream/My Drive/Github/Fiber-Analysis/data/L6.tsv", sep = "\t", header = T, check.names = F)
colnames(L6)[1] <- "sampleID"
L6_melt <- melt(L6, id.vars = "sampleID")
L6_merge_melt <- merge(L6_melt, metadata, by.y = "Metagenome", by.x = "sampleID")
L6_merge_melt$Individual <- as.factor(L6_merge_melt$Individual)

L6_taxabarplot <- ggplot(data = L6_merge_melt) +
  aes(x = as.factor(Day), fill = variable, weight = value) +
  geom_bar() +
  theme_bw() +
  facet_grid(cols = vars(Individual), space = "free", scales = "free") + scale_fill_manual(values = cb27) +
  guides(fill=guide_legend(title="Bacterial Genus")) +
  labs(x = 'Day',
       y = 'Relative Abundance') + theme(axis.ticks.x=element_blank())
