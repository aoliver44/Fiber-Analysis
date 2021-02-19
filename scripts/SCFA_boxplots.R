###################################
########## SCFA Analysis ##########
###################################

# path should be set to data directory
path <- getwd()
setwd(path)
source("../scripts/Generate_basic_env.R")

library(readxl)
library(reshape2)
library(tidyverse)
library(ggplot2)
library(ggpubr)

SCFA <- read_excel("~/Google Drive File Stream/My Drive/Github/Fiber-Analysis/data/fiber_diet_gcfid_data_eric.xlsx")
test_scfa <- SCFA %>% select(., ID, Intervention, norm_acetate, norm_propionate, norm_butyrate, norm_valerate) 
plot_data <- test_scfa %>%
  group_by(as.factor(ID), Intervention) %>% 
  summarise_each(funs(mean))
plot_data_melt <- melt(plot_data, measure.vars = c("norm_acetate", "norm_propionate", "norm_butyrate", "norm_valerate"))
plot_data_melt$Intervention <- factor(plot_data_melt$Intervention, levels = c("pre", "post"), ordered = T)
plot_data_melt$ID <- plot_data_melt$`as.factor(ID)`
plot_data_melt %>% group_by(Intervention, variable) %>% drop_na() %>% summarise(., mean = mean(value), sd = sd(value), cv = ((sd(value)/mean(value))*100))

ggplot(data = plot_data_melt, aes(x = Intervention, y = value)) +
  geom_boxplot(aes(x = Intervention, y = value, fill = Intervention)) + geom_point(aes(color= plot_data_melt$ID)) +
  theme_classic(base_size = 20) + scale_fill_brewer(palette = "Blues") + geom_line(aes(group = plot_data_melt$ID), alpha = 0.2, colour = "black", data = plot_data_melt) + facet_grid(facets = plot_data_melt$variable) + labs(color = "Individual", fill = "Diet Intervention", x = "", y = "SCFA concentration (mg/L)")

##############################################################################

# Statistics
library(nlme)
SCFA$tech_rep <- as.factor(SCFA$`Technical Replicate`)
SCFA$`Ethylbutyrate [ng/ul]` <- as.numeric(as.character(SCFA$`Ethylbutyrate [ng/ul]`))

anova(lme(norm_propionate ~ Intervention*tech_rep, data = SCFA, random = ~1|ID, na.action = na.omit))

# tech rep...better
temp_scfa <- SCFA %>% select(., ID, Sample, `Technical Replicate`, Intervention, norm_acetate, norm_propionate, norm_butyrate, norm_valerate, `Ethylbutyrate [ng/ul]`)
temp_scfa$ethyl_butyrate <- temp_scfa$`Ethylbutyrate [ng/ul]`
temp_scfa$`Ethylbutyrate [ng/ul]` <- NULL
replicated_ids <- c("1","4","7","8","9","11","12","13")
temp_scfa_1 <- subset(temp_scfa, ID %in% replicated_ids)
temp_scfa_2 <- melt(temp_scfa_1, measure.vars = c("norm_acetate", "norm_propionate", "norm_butyrate", "norm_valerate", "ethyl_butyrate"))
temp_scfa_2$ID <- as.factor(temp_scfa_2$ID)
temp_scfa_2$value <- as.numeric(as.character(temp_scfa_2$value))
temp_scfa_2$Sample <- as.factor(temp_scfa_2$Sample)
temp_scfa_2$tech_rep <- as.factor(temp_scfa_2$`Technical Replicate`)

ggplot(data = temp_scfa_2, 
       aes(x = Intervention, y = value, color = tech_rep, fill = tech_rep)) +
  geom_boxplot(position = "dodge") +  
  theme_classic() +
  facet_wrap(. ~ ID + variable, scales = "free_y")

########## ratio analysis ###########
# variation in tech rep, on average
ratio_analysis1 <- temp_scfa_2 %>% group_by(ID, Intervention, variable, tech_rep) %>% summarise_at(vars(value), mean, na.rm = T)
ratio_analysis2 <- ratio_analysis1 %>% group_by(tech_rep) %>% mutate(rn = row_number()) %>% spread(tech_rep, value)
ratio_analysis3 <- ratio_analysis2 %>% group_by(ID, variable) %>% 
  mutate(diff = abs(`1` - `2`)/(abs(`1` + `2`)/2)) %>% View()
  group_by(ID,variable) %>% summarise(mean_diff=mean(diff))
colnames(ratio_analysis3)[3] <- "Technical_diff"
# variation in biological signal, on averate
mean_analysis1 <- temp_scfa_2 %>% group_by(ID, Intervention, variable) %>% summarise_at(vars(value), mean, na.rm = T) 
mean_analysis2 <- mean_analysis1 %>% spread(Intervention, value)
mean_analysis3 <- mean_analysis2 %>% group_by(ID, variable) %>% 
  mutate(diff = (abs(post - pre)/(abs(post + pre)/2))) %>% 
  group_by(ID,variable) %>% 
  summarise(mean_diff=mean(diff))
colnames(mean_analysis3)[3] <- "Biological_diff"

total_diff <- merge(mean_analysis3, ratio_analysis3, by = c("ID", "variable"))
total_diff <- total_diff %>% mutate(diff = (Biological_diff - Technical_diff))
total_diff$color[total_diff$diff < 0 ] <- "Technical noise higher than biological signal"
total_diff$color[total_diff$diff > 0 ] <- "Biological signal higher than technical noise"

ggplot(data = subset(total_diff, total_diff$variable != "ethyl_butyrate")) +
  aes(x = variable, weight = diff, fill = color) +
  geom_bar(position = "dodge") +
  theme_bw() +
  labs(y = '') +
  facet_grid(cols = vars(ID)) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values = c("#39A9AB","#FF7C43"))

#########################################################################

# ethyl butyrate subset
ggplot(data = subset(temp_scfa_2, temp_scfa_2$variable == "ethyl_butyrate"), 
       aes(x = Intervention, y = value, fill = tech_rep)) +
  geom_boxplot(position = "dodge") + geom_point(position = position_jitterdodge(jitter.width = .1)) +
  theme_classic() + 
  facet_wrap(. ~ ID + variable) #+ stat_compare_means(method = "t.test")
#+ theme(axis.text.y = element_blank()) 


