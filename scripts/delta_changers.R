##################################
######### Delta Changers #########
##################################

source("~/Google Drive File Stream/My Drive/Github/Fiber-Analysis/scripts/Generate_basic_env.R")

library(ggplot2)
library(corrr)
library(matrixStats)
library(tidyverse)
library(reshape2)
library(ggpubr)
library(rstatix)
library("viridis")
library(ggpubr)

metadata_otu <- merge(metadata, alpha_rare, by.x = "Metagenome", by.y = "row.names")
metadata_otu$Individual_level <- paste0("X",metadata_otu$Individual)

# pvalues
delta_abundances <- metadata_otu %>%
  select(., Metagenome, Individual_level, Intervention, 13:391) %>%
  melt(., id.vars = c("Metagenome", "Individual_level", "Intervention")) %>%
  group_by(., Individual_level, Intervention, variable) %>%
  summarise(., mean_abund = mean(value)) %>% 
  spread(., Intervention, mean_abund) %>%
  mutate(., diff = post - pre) %>% 
  drop_na() %>%
  dcast(Individual_level ~ variable, value.var = "diff") %>% 
  na_if(., 0) %>%
  #be present in at least 25% of samples
  discard(~sum(is.na(.x))/length(.x)* 100 >=75) %>% 
  column_to_rownames(var = "Individual_level") %>% 
  replace(is.na(.), 0) %>%
  cor_test(method = "spearman")
delta_p_values2 <- as.data.frame(p.adjust(delta_abundances$p, method = "fdr"))
delta_p_values3 <- cbind(delta_p_values2, delta_abundances) %>% arrange(., cor)
colnames(delta_p_values3)[1] <- "padjust"

# Average abundances
avg_abund <- metadata_otu %>%
  filter(., Intervention == "pre") %>%
  select(., Metagenome, 14:382) %>%
  column_to_rownames(., "Metagenome") %>%
  summarise_all(funs(median)) %>% melt()

change_in_abund <- merge(delta_p_values3, avg_abund, by.x = "var1", by.y = "variable")
colnames(change_in_abund)[8] <- "avg_pre_abund"

# filter by avg_pre_abund and FDR and correlation value

hm.palette <- colorRampPalette(rev(brewer.pal(9, 'YlOrRd')), space='Lab')

change_in_abund %>%
  filter(., avg_pre_abund > 1, padjust < 0.2, abs(cor) != 1, abs(cor) > 0.2) %>% 
  arrange(desc(cor)) %>% 
  ggplot() + aes(x = var1, y = var2, fill = cor) +
  geom_tile() + coord_flip() + theme(axis.text.x = element_text(angle = 90)) + 
  labs(x = "", y= "") +
  scale_fill_viridis(option = "D") + ggtitle(label = "Spearman correlation between change in abundances")

raw_melt <- metadata_otu %>%
    select(., Metagenome, Individual_level, Intervention, 13:391) %>%
    melt(., id.vars = c("Metagenome", "Individual_level", "Intervention")) %>%
    group_by(., Individual_level, Intervention, variable) %>%
    summarise(., mean_abund = (mean(value)/900)) %>% 
    spread(., Intervention, mean_abund) %>%
    mutate(., diff = post - pre) %>% 
    select(., diff, Individual_level, variable) %>%
    #filter(., variable == "Bifidobacterium_longum_57796" | variable == "Lactobacillus_mucosae_57556") %>%
    spread(., variable, diff) %>% drop_na()

# check raw data
ggplot(raw_melt, 
       aes(x=Bifidobacterium_pseudocatenulatum_57754, y=Ruminococcus_torques_62045)) + 
  geom_point() + geom_smooth(method=lm) + stat_cor(method = "spearman") + ggtitle("Change in abundance in response to Fiber Intervention")


# correlate changers butyrate producers
butyrate_corr <- delta_abundances %>%
  filter(str_detect(var1, c("_praus", "_rectale", "faecis", "_intestinalis", "_inulin", "halii", "hadrus")))
butyrate_p_values2 <- as.data.frame(p.adjust(butyrate_corr$p, method = "fdr"))
butyrate_p_values3 <- cbind(butyrate_p_values2, butyrate_corr) %>% arrange(., cor)
colnames(butyrate_p_values3)[1] <- "padjust"

but_change_in_abund <- merge(butyrate_p_values3, avg_abund, by.x = "var1", by.y = "variable")
colnames(but_change_in_abund)[8] <- "avg_pre_abund"

# fun fact, like nothing is significant when corrected for, hence no p.adjust parameter
but_change_in_abund %>%
  filter(., avg_pre_abund > 1, abs(cor) != 1, abs(cor) > 0.3) %>% 
  arrange(desc(cor)) %>% 
  ggplot() + aes(x = var1, y = var2, fill = cor) +
  labs(x = "Butyrate Producers", y= "") +
  geom_tile() + coord_flip() + theme(axis.text.x = element_text(angle = 90)) + 
  scale_fill_viridis(option = "D")

## correlate with butyrate itself 
diff_otu <- metadata_otu %>%
  select(., Metagenome, Individual_level, Intervention, 13:391) %>%
  melt(., id.vars = c("Metagenome", "Individual_level", "Intervention")) %>%
  group_by(., Individual_level, Intervention, variable) %>%
  summarise(., mean_abund = mean(value)) %>% 
  spread(., Intervention, mean_abund) %>%
  mutate(., diff = post - pre) %>%
  drop_na() %>%
  dcast(Individual_level ~ variable, value.var = "diff") %>% View()
  #be present in at least 25% of samples
  discard(~sum(is.na(.x))/length(.x)* 100 >=75) %>% 
  column_to_rownames(var = "Individual_level") 
