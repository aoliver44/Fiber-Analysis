#######################################
######### Functional Analysis #########
#######################################


## CAZY ###

source("~/Google Drive File Stream/My Drive/Github/Fiber-Analysis/scripts/Generate_basic_env.R")
library(reshape2)
library(tidyverse)
library(nlme)
library(ggplot2)
library(ggpubr)
library(knitr)
library(kableExtra)

# Import the cazy results and make into an OTU table
cazy <- read.csv("CAZY_classified_final.txt", check.names = F, header = F, sep = "\t")
cazy_meta <- merge(metadata, cazy, by.x = "Metagenome", by.y = "V1")
cazy_meta <- cazy_meta %>% separate(., V2, sep = "_", extra = "drop", into = c("subclass", "sub-sub"), remove = F)
cazy_meta <- cazy_meta %>% mutate(., V4 = gsub("_", "\n", V4))

## boxplots of major CAzy classes faceted by individual
ggplot(data = cazy_meta) + 
  aes(x = Intervention, y = log2(V3)) +
  geom_boxplot(aes(fill = Intervention)) +
  geom_point(aes(color = as.factor(Individual)), size = 0.4, position = position_jitterdodge()) +
  labs(title = 'Carbohydrate Active Enzymes',
       x = 'Individual',
       y = 'Normalized Gene Counts (log2 transformed)') +
  theme_bw() + theme(axis.text.x = element_text(angle = 90)) +
  facet_grid(V4 ~ Individual) + scale_fill_manual(values=c("deepskyblue3", "orange"))

# Or:
cazy_meta %>% group_by(., Individual, V4, Intervention, Metagenome) %>% summarise(., sum = sum(V3)) %>% 
  group_by(V4, Individual, Intervention) %>% summarise(mean = mean(sum)) %>% 
  ggplot() + 
  aes(x = Intervention, y = mean) +
  geom_boxplot(aes(fill = Intervention), outlier.shape = NA) +
  geom_point(aes(color = as.factor(Individual)), size = 0.8) +
  geom_line(aes(group = Individual), alpha = 0.3) +
  labs(title = 'Carbohydrate Active Enzymes',
       x = '',
       y = 'Mean normalized reads mapping to CAzy') +
  theme_bw() + theme(axis.text.x = element_text(angle = 90)) +
  facet_grid(. ~ V4) + scale_fill_manual(values=c("deepskyblue3", "orange"))

## Table of the amount of enzymes and reads that were measured by class
cazy_meta %>% group_by(., Intervention, V4) %>% 
  rename(., Enzyme_Class = V4) %>%
  summarise(., read_count = round(sum(V3)), 
            distinct_enzymes = n_distinct(V2), 
            distinct_Families = n_distinct(subclass)) %>% 
  kable() %>%
  kable_styling(bootstrap_options = "striped", full_width = F, position = "left")

# Richness
cazy_meta %>% 
  filter(., V4 == "Glycoside\nhydrolase") %>%
  group_by(., Intervention, Individual) %>% 
  summarise(., Distinct_enzymes = n_distinct(subclass)) %>%
  #summarise(., mean(Distinct_enzymes))
  ggplot() + aes(x = Intervention, y = Distinct_enzymes) + 
  geom_boxplot(aes(fill = Intervention)) + 
  geom_point(aes(color = as.factor(Individual)), position = position_jitterdodge(dodge.width = 0.3)) + 
  theme_bw(base_size = 14) +
  scale_fill_brewer(palette = "Blues")

## Distribution of enzymes within Major class
ggplot(data = subset(cazy_meta, V4 == "Glycoside\nhydrolase")) +
  aes(reorder(subclass, V3), V3, fill = Intervention, color = Intervention) +
  geom_point(size = 0.2) +
  theme_bw(base_size = 7) + coord_flip() + 
  labs(x = "Glycoside hydrolase", y = "Normalized Read Counts") +
  scale_color_manual(values=c("deepskyblue3", "orange"))

#########################################################################
# multivariate stats on CAzy

cazy_otu <- cazy %>% 
  dcast(., formula = V1 ~ V2, value.var = "V3") %>%
  column_to_rownames(., "V1") %>%
  replace(., is.na(.), 0)
cazy_otu <- cazy_otu[, colSums(cazy_otu) >= 10]
# Relative abundance:
#cazy_otu <- decostand(cazy_otu, method = "total")

# Random Forest

library(rfPermute)
set.seed(999)
Fiber_rfp_data <- merge(metadata, cazy_otu, by.x = "Metagenome", by.y = "row.names")
RFP_data_coarse <- Fiber_rfp_data %>% dplyr::select(3, (13:NCOL(Fiber_rfp_data)))
Fiber_RFP <- rfPermute(as.factor(Intervention) ~., data = RFP_data_coarse, proximity = TRUE, importance = TRUE, corr.bias = TRUE, mtry = 60, ntree = 901, num.cores = 6)
# plot the plots
proximityPlot(Fiber_RFP)
var_imp <- varImpPlot(Fiber_RFP, type = 1, n.var = 10)
plotConfMat(Fiber_RFP)
heatmap <- impHeatmap(Fiber_RFP, alpha = 0.05, ranks = F, n = 6)
heatmap + theme(axis.text.x = element_text(angle = 90, hjust = 1))
plot(Fiber_RFP, alpha = 0.05)

# check the raw data

library(ggpubr)
ggplot(data = Fiber_rfp_data) +
  aes(x = Intervention, y = Fiber_rfp_data$GH120) +
  geom_boxplot() + geom_point(position = "jitter", width = 0.2, aes(colour = as.factor(Fiber_rfp_data$Individual))) +
  theme_bw() + stat_compare_means(method = "anova") + 
  #facet_grid(. ~ Intervention) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Permanova

library(vegan)

intervention_permanova <- adonis(Fiber_rfp_data[,13:281] ~ Individual*Intervention, data = Fiber_rfp_data, permutations = 999, parallel = 4, method = "euclidean")
intervention_permanova

coef <- coefficients(intervention_permanova)["Intervention.L",]
top.coef <- coef[rev(order(abs(coef)))[1:20]]

par(mar=c(3,9,3,2) + 0.1)
barplot(sort(top.coef), horiz = T, las = 1, main = "Top taxa", col = "deepskyblue3", cex.axis=1, cex.names=0.8, xlim = c(-25,250))
dev.off()


# KW on broad cazy groups by Cazy "richness"
GH_lme <- cazy_meta %>% 
  filter(., V4 == "Glycoside\nhydrolase") %>% 
  group_by(., Intervention, Individual) %>% 
  summarise(., Distinct_enzymes = n_distinct(subclass)) 
shapiro.test(GH_lme$Distinct_enzymes)
anova(aov(as.numeric(Distinct_enzymes) ~ Intervention, data = GH_lme))

# LME on broad cazy groups by Cazy "abundance"
GH_lme <- cazy_meta %>% 
  filter(., V4 == "Carbohydrate\nesterase") %>%
  group_by(., Individual, V4, Intervention, Metagenome) %>% summarise(., sum = sum(V3)) %>% 
  group_by(V4, Individual, Intervention) %>% summarise(Abundance = mean(sum))

wilcox.test(as.numeric(Abundance) ~ Intervention, data = GH_lme)

###############################################################################################################

## MIDAS COMPOUND SEARCH

midas_scfa <- read.csv("midas_butyrate.txt", sep = "\t", header = T)
midas_scfa_meta <- merge(midas_scfa, metadata, by.x = "sample_id", by.y = "Metagenome")
midas_scfa_meta <- midas_scfa_meta %>% mutate(., species_id = gsub("_", "\n", species_id))

butyrate <- read.csv("midas_butyrate.txt", sep = "\t", header = T)
butyrate_meta <- merge(butyrate, metadata, by.x = "sample_id", by.y = "Metagenome")
butyrate_meta <- butyrate_meta %>% mutate(., SCFA = "Butyrate")
acetate <- read.csv("midas_acetate.txt", sep = "\t", header = T)
acetate_meta <- merge(acetate, metadata, by.x = "sample_id", by.y = "Metagenome")
acetate_meta <- acetate_meta %>% mutate(., SCFA = "Acetate")
propionate <- read.csv("midas_propionate.txt", sep = "\t", header = T)
propionate_meta <- merge(propionate, metadata, by.x = "sample_id", by.y = "Metagenome")
propionate_meta <- propionate_meta %>% mutate(., SCFA = "Propionate")
lactate <- read.csv("midas_lactate.txt", sep = "\t", header = T)
lactate_meta <- merge(lactate, metadata, by.x = "sample_id", by.y = "Metagenome")
lactate_meta <- lactate_meta %>% mutate(., SCFA = "Lactate")
ethanol <- read.csv("midas_ethanol.txt", sep = "\t", header = T)
ethanol_meta <- merge(ethanol, metadata, by.x = "sample_id", by.y = "Metagenome")
ethanol_meta <- ethanol_meta %>% mutate(., SCFA = "Ethanol")

all_scfa_melt <- rbind(butyrate_meta, acetate_meta, propionate_meta, lactate_meta, ethanol_meta)

# Figure of pathway reads
# normalize to reads
decon_reads <- read.csv("decon_reads.txt", sep = " ")
decon_reads$total_reads <- decon_reads$decon_reads * 2
decon_reads$norm_factor <- as.numeric(10673396 / decon_reads$total_reads)

all_scfa_melt <- merge(all_scfa_melt, decon_reads, by.x = "sample_id", by.y = "FILE_qc_decon_pe.1.fastq.gz")
all_scfa_melt$norm_reads <- all_scfa_melt$count_reads * all_scfa_melt$norm_factor
all_scfa_melt$Intervention <- factor(all_scfa_melt$Intervention, levels = c("pre", "post"), ordered = T)
all_scfa_melt$rel_abund <- all_scfa_melt$count_reads / all_scfa_melt$total_reads

give.n <- function(x){
  return(c(y = median(x)*1.05, label = length(x))) 
  # experiment with the multiplier to find the perfect position
}

all_scfa_melt %>%
  group_by(., Individual, Intervention, SCFA) %>%
  summarise(., total_norm_reads = mean(rel_abund)) %>% 
  filter(., SCFA %in% c("Butyrate", "Propionate", "Acetate")) %>%
  ggplot() + aes(Intervention, total_norm_reads) +
  geom_boxplot(aes(fill = Intervention), outlier.alpha = 0) +
  geom_jitter(width = 0.2, alpha = 0.4) +
  scale_fill_brewer(palette = "Blues") +
  facet_wrap(~SCFA, scales = "free") +
  theme_bw() + labs(y = "Relative Abundance") #+ 
  stat_summary(fun.data = give.n, geom = "text", fun.y = median, 
                                                            position = position_dodge(width = 0.75))


# LME on different compounds:
pathway_lme <- all_scfa_melt %>%
  group_by(., Individual, Intervention, SCFA) %>%
  summarise(., total_norm_reads = mean(rel_abund)) %>% 
  filter(., SCFA == "Butyrate")

wilcox.test(total_norm_reads ~ Intervention, data = pathway_lme)

# Test significance of intervention on rel abund of specific butyrate genes
butyrate_gene <- all_scfa_melt %>% 
  filter(., enzyme_id == "2.8.3.8") %>% 
  group_by(., Individual, Intervention, SCFA) %>% 
  summarise(., total_norm_reads = mean(rel_abund))
shapiro.test(butyrate_gene$total_norm_reads)
wilcox.test(total_norm_reads ~ Intervention, butyrate_gene)

# Midas reads by pathway, by organism
all_scfa_melt %>% filter(., SCFA == "Butyrate") %>%
  #group_by(., species_id, Intervention) %>%
  #summarise(., total_reads = sum(norm_reads)) %>%
  #filter(str_detect(species_id, "Bifido")) %>% 
ggplot() +
  aes(x = Intervention, y = norm_reads, fill = Intervention) +
  geom_boxplot() +
  theme_bw(base_size = 10) + theme(axis.text.x = element_text(angle = 90)) +
  facet_wrap(~enzyme_id + species_id, scales = "free_y") + 
  scale_fill_manual(values=c("deepskyblue3", "orange"))

# Midas table: which genes are repersented above and how many reads?
midas_scfa_meta %>% group_by(enzyme_id) %>% 
  summarise(., read_count = sum(count_reads)) %>% 
  kable() %>% 
  kable_styling(bootstrap_options = "striped", full_width = F, position = "left") %>%
  add_header_above(c("butyrate", ""))


# Permanova with midas:
permanova_scfa <- all_scfa_melt %>%
  dcast(sample_id ~ enzyme_id, value.var = "norm_reads") 
permanova_scfa_meta <- merge(metadata, permanova_scfa, by.x = "Metagenome", by.y = "sample_id")

adonis(permanova_scfa_meta[,13:36] ~ Individual*Intervention, data = permanova_scfa_meta, permutations = 999, parallel = 4, method = "bray")
