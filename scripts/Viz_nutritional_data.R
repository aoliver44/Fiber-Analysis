##################################
###### Viz Nutritional data ######
##################################

# path should be set to data directory
path <- getwd()
setwd(path)
source("../scripts/Generate_basic_env.R")

# Displaying the Nutritional data over time period of the 3 weeks. 
# The resulting lmes are done over the 3 week period as well.

library(ggplot2)
library(cowplot)
library(nlme)
library(lme4)
library(tidyverse)
# plot the data
fiber <- ggplot(data = Nutritional_Data) +
  aes(x = as.factor(week), y = Fiber, fill = Intervention) +
  geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge(jitter.width = 0.6), alpha = 0.07) +
  #geom_jitter(width = 0.15, alpha = 0.2) +
  scale_fill_brewer(palette = "Blues") +
  labs(x = '',
       y = 'Fiber (grams/day)') +
  theme_bw(base_size = 14) + theme(legend.position = "none")

carbs <- ggplot(data = Nutritional_Data) +
  aes(x = as.factor(week), y = Carbohydrates, fill = Intervention) +
  geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge(jitter.width = 0.6), alpha = 0.07) +
  #geom_jitter(width = 0.15, alpha = 0.2) +
  scale_fill_brewer(palette = "Purples") +
  labs(x = '',
       y = 'Carbohydrates (grams/day)') +
  theme_bw(base_size = 14) + theme(legend.position = "none")

protein <- ggplot(data = Nutritional_Data) +
  aes(x = as.factor(week), y = Protein, fill = Intervention) +
  geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge(jitter.width = 0.6), alpha = 0.07) +
  #geom_jitter(width = 0.15, alpha = 0.2) +
  scale_fill_brewer(palette = "Greens") +
  labs(x = '',
       y = 'Protein (grams/day)') +
  theme_bw(base_size = 14) + theme(legend.position = "none")

fat <- ggplot(data = Nutritional_Data) +
  aes(x = as.factor(week), y = Fat, fill = Intervention) +
  geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge(jitter.width = 0.6), alpha = 0.07) +
  #geom_jitter(width = 0.15, alpha = 0.2) +
  scale_fill_brewer(palette = "Reds") +
  labs(x = '',
       y = 'Fat (grams/day)') +
  theme_bw(base_size = 14) + theme(legend.position = "none")

calories <- ggplot(data = Nutritional_Data) +
  aes(x = as.factor(week), y = Calculated_calories, fill = Intervention) +
  geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge(jitter.width = 0.6), alpha = 0.07) +
  #geom_jitter(width = 0.15, alpha = 0.2) +
  scale_fill_brewer(palette = "Greys") +
  labs(x = '',
       y = 'Calories (kcal/day)') +
  theme_bw(base_size = 14) + theme(legend.position = "none")

plot_grid(fiber, labels = "A")
plot_grid(carbs, protein, fat, calories, labels= c("B","C", "D", "E"), ncol = 2)


############### STATISTICS #########################

# There are probs a lot of ways to get after this info, but here are LMEs on
# week ~ Nutrition, taking into acccount repeated measures.
nut.mod.rand <- lme(Fiber ~ as.factor(week), random = (~1|as.factor(Individual)), data = Nutritional_Data, method = "REML")
anova(nut.mod.rand)
shapiro.test(nut.mod.rand$residuals)
# none of the nutritional data is normally distributed
library(car)
qqPlot(Nutritional_Data$Fiber)

# Average out with the KW?
Nutritional_kw <- Nutritional_Data %>%
  select(., Individual, week, Carbohydrates, Intervention) %>%
  group_by(., Intervention, Individual) %>%
  summarise(., Average = mean(Carbohydrates))
shapiro.test(Nutritional_Data$Protein)
kruskal.test(Average ~ Intervention, data = Nutritional_kw)
wilcox.test(Average ~ Intervention, data = Nutritional_kw)
