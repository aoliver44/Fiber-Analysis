############################
##### Humann3 analysis #####
############################

# path should be set to data directory
path <- getwd()
setwd(path)
source("../scripts/Generate_basic_env.R")

library(reshape2)
library(ggplot2)
library(tidyverse)


humann3_raw <- read.csv("scfa_humann_subset.txt", sep = "\t")
humann3_meta <- read.csv("humann3_path_meta.txt", sep = "\t")
humann3_sample_meta <- read.csv("humann_sample_meta.txt", sep = "\t")

humann3_melt <- melt(data = humann3_raw)
humann3_melt_meta <- merge(humann3_melt, humann3_meta, by.x = "pathway", by.y = "pathway")
humann3_melt_meta <- merge(humann3_melt_meta, humann3_sample_meta, by.x = "variable", by.y = "sample")
humann3_melt_meta$intervention <- factor(humann3_melt_meta$intervention, levels=c("pre", "post"))
humann3_melt_meta$`SCFA.1` <- factor(humann3_melt_meta$`SCFA.1`, levels=c("Acetate", "Butyrate", "Lactate", "Propionate", "Multiple_pwys"))

humann3_melt_meta_nozero <- subset(humann3_melt_meta, humann3_melt_meta$value > 0)
humann3_melt_meta_nozero <- subset(humann3_melt_meta_nozero, humann3_melt_meta_nozero$Organism != "unclassified")

# # outside of project, where package ggh4x installed correctly
# # renv apparently didnt like it
# library(ggh4x)
# library(ggpubr)
# ggplot(data = humann3_melt_meta_nozero, aes(y = log(value), x = metacyc, fill = Organism)) +
#   geom_boxplot(outlier.alpha = 0) + facet_nested(. ~ `SCFA.1` + intervention + Organism, scales = "free") + theme_bw() +
#   geom_point(position = position_jitter(width = 0.1), alpha = 0.2) +
#   theme(panel.spacing=unit(0,"lines"),
#         strip.background=element_rect(color="grey30", fill="grey90"),
#         panel.border=element_rect(color="grey90")) +
#   theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))

humann3_melt_meta_nozero %>%
  filter(., Organism == "Total") %>% 
  group_by(., metacyc, intervention) %>%
  tally()

ggplot(data = subset(humann3_melt_meta_nozero, humann3_melt_meta_nozero$Organism == "Total"), aes(y = log10(value), fill = intervention)) +
  geom_boxplot(outlier.alpha = 0) + facet_grid(. ~ `SCFA.1` + metacyc, scales = "free") + theme_bw() +
  #geom_point(position = position_jitter(width = 0.1), alpha = 0.2) +
  theme(panel.spacing=unit(0.1,"lines"),
        strip.background=element_rect(color="grey30", fill="grey90"),
        panel.border=element_rect(color="grey90"), panel.grid.minor.x = element_blank(), 
        panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(axis.text.x = element_blank(), axis.ticks.x.bottom = element_blank()) + scale_fill_brewer(palette = "Blues")

tmp <- humann3_melt_meta_nozero %>%
  filter(., metacyc == "hexitol fermentation to lactate, formate, ethanol and acetate") %>%
  filter(., Organism == "g__Anaerostipes.s__Anaerostipes_hadrus") %>%
  group_by(., individual, intervention) %>%
  tally()
  


# anova assumptions
anova_subset <- subset(humann3_melt_meta_nozero, humann3_melt_meta_nozero$metacyc == "hexitol fermentation to lactate, formate, ethanol and acetate")
anova_subset <- subset(anova_subset, anova_subset$Organism == "g__Anaerostipes.s__Anaerostipes_hadrus")
shapiro.test(anova_subset$value)
anova_subset <- anova_subset %>% group_by(., individual, intervention) %>% summarise(., mean = mean(value))
shapiro.test(anova_subset$mean)
wilcox.test(anova_subset$mean ~ anova_subset$intervention)

