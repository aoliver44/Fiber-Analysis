#################################
######## Bifido Barplots ########
#################################

# path should be set to data directory
path <- getwd()
setwd(path)
source("../scripts/Generate_basic_env.R")

library(ggplot2)
library(corrr)
library(tidyverse)
library(reshape2)
library(ggpubr)
library(rstatix)

# Relative abundance of Bifido only (bifido sum to 1) use this file: Condensed_otu_table_bifido.tsv
# now do a little stuff outside of R, for the final taxonomy stuff
taxonomy_melted <- read.table("~/Google Drive File Stream/My Drive/Github/Fiber-Analysis/data/taxonomy_melted.tsv", sep = "\t", header = T)
taxonomy_melted_bifido <- subset(taxonomy_melted, L6 == "Bifidobacterium")
taxonomy_melt_meta_bifido <- merge(taxonomy_melted_bifido, metadata, by.x = "variable", by.y = "Metagenome")
taxonomy_melt_meta_bifido$Individual <- as.factor(taxonomy_melt_meta_bifido$Individual)
taxonomy_melt_meta_bifido$Day <- as.factor(taxonomy_melt_meta_bifido$Day)
taxonomy_melt_meta_bifido$L7 <- gsub(pattern = "Bifidobacterium Bifido", replacement = "Bifido", x = taxonomy_melt_meta_bifido$L7)
cb12 <- c("#015888","#ffb33c","#6c40c5","#007019", "#a90ca0","#4fdad2","#cf005d", "#02bcf2", "#ed48c2","#704d00","#b79eff","#feafd6","#673e95")

# # Geom-Area plot
# ggplot(data = taxonomy_melt_meta_bifido) + aes(x = as.numeric(Day), y = value, fill = L7) +
#   geom_area() +
#   #geom_vline(aes(xintercept = 3), color = "grey", linetype = "dashed", size = 0.75) +
#   theme_bw() + facet_grid(. ~ Individual, space = "free", scales = "free") + 
#   scale_fill_manual(values = cb12) +
#   theme(panel.spacing = unit(0.1, "lines")) +   
#   theme(axis.text.x=element_blank(),
#         axis.ticks.x=element_blank())

# Geom-bar plot
ggplot(data = taxonomy_melt_meta_bifido) +
  aes(x = as.factor(Day), fill = L7, weight = value) +
  geom_bar() +
  theme_bw() + facet_grid(. ~ Individual, space = "free", scales = "free") + 
  scale_fill_manual(values = cb12) +
  labs(x = "Day", y = 'Relative Abundance (MIDAS)') +
  theme(panel.spacing = unit(0.1, "lines")) +   
  theme(axis.ticks.x=element_blank()) + guides(fill=guide_legend(title="Bifidobacterium"))


### Geom-bar plot: top 5 with alex
taxonomy_melt_meta_bifido_5 <- taxonomy_melt_meta_bifido
taxonomy_melt_meta_bifido_5 %>% group_by(., L7) %>% summarise(., mean = mean(value)) %>% arrange(., mean)
low_abundance <- c("Bifidobacterium biavatii", "Bifidobacterium pseudolongum", "Bifidobacterium dentium", "Bifidobacterium moukalabense", "Bifidobacterium animalis subsp. lactis", "Bifidobacterium ruminantium", "Bifidobacterium breve")
taxonomy_melt_meta_bifido_5$L7[taxonomy_melt_meta_bifido_5$L7 %in% low_abundance == "TRUE"] <- "other"

ggplot(data = taxonomy_melt_meta_bifido_5) +
  aes(x = as.factor(Day), fill = L7, weight = value) +
  geom_bar() +
  theme_bw() + facet_grid(. ~ Individual, space = "free", scales = "free") + 
  scale_fill_manual(values = cb12) +
  labs(x = "Day", y = 'Relative Abundance (MIDAS)') +
  theme(panel.spacing = unit(0.1, "lines")) +   
  theme(axis.ticks.x=element_blank()) + guides(fill=guide_legend(title="Bifidobacterium"))

###########################################################################################

## correlations with bifido + other taxa
taxonomy_melted_meta_cor <- merge(metadata, taxonomy_melted, by.x = "Metagenome", by.y = "variable")
taxonomy_melted_meta_cor %>% group_by(., L6) %>% summarise(., sd = sd(value)) %>% arrange(., sd)
taxonomy_melted_meta_cor <- taxonomy_melted_meta_cor %>% select(., -matches("Bacillus"), -matches("Listeria"))
taxonomy_melted_meta_cor %>%
  mutate_all(na_if,"") %>%
  drop_na() %>% 
  group_by(., L6, Metagenome) %>% 
  summarise(., L6_amount = sum(value)) %>% 
  dcast(Metagenome ~ L6, value.var = "L6_amount") %>% 
  column_to_rownames(var = "Metagenome") %>% 
  corrr::correlate(method = "spearman") %>% 
  corrr::focus(Bifidobacterium) %>% 
  # cut off of 0.318 determined below...keeps all padjusted sig stuff in
  filter(abs(Bifidobacterium) > 0.318) %>% 
  mutate(rowname = factor(rowname, levels = rowname[order(Bifidobacterium)])) %>%
  ggplot(aes(x = rowname, y = Bifidobacterium)) +
  geom_bar(stat = "identity") +
  ylab("Correlation with Bifidobacterium\n (Spearman)") +
  xlab("Genus") + theme_bw() + theme(axis.text.x = element_text(angle = 90)) +
  coord_flip()

# pvalues
Bifido_p_values <- taxonomy_melted_meta_cor %>%
  mutate_all(na_if,"") %>%
  drop_na() %>% 
  group_by(., L6, Metagenome) %>% 
  summarise(., L6_amount = sum(value)) %>% 
  dcast(Metagenome ~ L6, value.var = "L6_amount") %>% 
  column_to_rownames(var = "Metagenome") %>% 
  cor_test(method = "spearman", vars = "Bifidobacterium")
Bifido_p_values <- Bifido_p_values %>% filter(., var1 == "Bifidobacterium") 
Bifido_p_values2 <- as.data.frame(p.adjust(Bifido_p_values$p, method = "fdr"))
Bifido_p_values3 <- cbind(Bifido_p_values2, Bifido_p_values) %>% arrange(., cor)
colnames(Bifido_p_values3)[1] <- "padjust"


## check with raw data
raw_genus_check <- taxonomy_melted_meta_cor %>%
  mutate_all(na_if,"") %>%
  drop_na() %>% 
  group_by(., L6, Metagenome) %>% 
  summarise(., L6_amount = sum(value)) %>% 
  dcast(Metagenome ~ L6, value.var = "L6_amount") %>% 
  column_to_rownames(var = "Metagenome")

microbe2check <- "Eubacterium"

ggscatter(raw_genus_check, x = microbe2check, y = "Bifidobacterium", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = paste0(microbe2check, " read count"), ylab = "Fiber Intake (grams)")
####################################################################################

## correlations with fiber
fiber_corr <- merge(metadata, alpha_rare, by.x = "Metagenome", by.y = "row.names")
fiber_corr %>%
  select(., Metagenome, Fiber, 13:391) %>%
  column_to_rownames(var = "Metagenome") %>%
  corrr::correlate(method = "spearman") %>% 
  corrr::focus(Fiber) %>% 
  filter(abs(Fiber) > 0.2) %>% 
  mutate(rowname = factor(rowname, levels = rowname[order(Fiber)])) %>% 
  ggplot(aes(x = rowname, y = Fiber)) +
  geom_bar(stat = "identity") +
  ylab("Correlation (R) with Fiber\n (Spearman)") +
  xlab("Species") + theme_bw() + theme(axis.text.x = element_text(angle = 90)) +
  coord_flip()

# pvalues - everythin is not significant when corrected
Fiber_p_values <- fiber_corr %>%
  select(., Metagenome, Fiber, 13:391) %>%
  column_to_rownames(var = "Metagenome") %>%
  cor_test(method = "spearman", vars = "Fiber")
Fiber_p_values <- Fiber_p_values %>% filter(., var1 == "Fiber") 
Fiber_p_values2 <- as.data.frame(p.adjust(Fiber_p_values$p, method = "fdr"))
Fiber_p_values3 <- cbind(Fiber_p_values2, Fiber_p_values) %>% arrange(., cor)
colnames(Fiber_p_values3)[1] <- "padjust"

## check with raw data
microbecheck <- "Bifidobacterium_psuedocatenulatum_57754"

fiber_corr$Individual <- as.factor(fiber_corr$Individual)
ggscatter(data = fiber_corr, x = microbecheck, y = "Fiber",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = paste0(microbe2check, " read count"), ylab = "Fiber Intake (grams)")


# supplemental figure on genus Bifido
taxonomy_melt_meta_bifido %>% 
  group_by(., Individual, Intervention) %>% 
  summarize(., total = mean(value)) %>% 
  ggplot() + aes(x = Intervention, y = as.numeric(total)) + geom_boxplot() + 
  geom_jitter(width = 0.1, aes(color = Individual)) + 
  theme_classic(base_size = 14) + labs(y = "Midas Read Count (mean)")

# supplemental figure on genus Bifido
taxonomy_melt_meta_bifido %>% 
  group_by(., L7, Intervention, Individual) %>% 
  summarize(., total = mean(value)) %>% 
  ggplot() + aes(x = Intervention, y = log2(as.numeric(total))) + geom_boxplot() + facet_wrap(.~L7, scales = "free_y") +
  theme_classic() + labs(y = "Midas Read Count \n(mean, log2 transformed)")


taxonomy_melt_meta_bifido %>% 
  group_by(., L7, Intervention) %>% 
  summarize(., total = mean(value)) %>% View()

