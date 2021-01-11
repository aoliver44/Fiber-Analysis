##################################
###### Generate Basic Files ######
##################################

setwd("~/Google Drive File Stream/My Drive/Github/Fiber-Analysis/data")
set.seed(999)

## Read in sequence data, nutritional data/metadata, metabolite data
library(readxl)
# Nutritional Data
Nutritional_Data <- read_excel("All_nutritional_data.xlsx")
Nutritional_Data$Intervention <- factor(Nutritional_Data$Intervention,
                                        levels = c('pre','post'),ordered = TRUE)
Nutritional_Data <- na.omit(Nutritional_Data)

# Sequence Data
midas <- t(read.csv("count_reads.txt", 
                    row.names = 1, check.names = FALSE, sep = "\t"))


# Metabolite data
SCFA <- read_xlsx("fiber_diet_gcfid_data_eric.xlsx")
SCFA$highlow <- ordered(SCFA$Intervention, levels = c("pre","post"))
SCFA$ID <- as.character(SCFA$ID)

# Metadata
metadata <- read.csv("Metadata.txt", sep = "\t")
metadata$Intervention <- factor(metadata$Intervention,
                                levels = c('pre','post'),ordered = TRUE)


## Calculating mean nutritions for pre-post intervention and assigning it to metadata

lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)
library(tidyverse)
library(fpc)
library(cluster)

Nutritional_means <- Nutritional_Data %>%
  group_by(., Individual, Intervention) %>%
  summarise_at(vars(Fiber:Calculated_calories), mean, na.rm = TRUE)


Nutritional_means$diff <- unlist(tapply(Nutritional_means$Fiber, INDEX = Nutritional_means$Individual,
                                        FUN = function(x) c(0, diff(as.numeric(x)))))

post_nutitional_means <- subset(Nutritional_means, Nutritional_means$Intervention == "post")

# pam clustering to get number of clusters of differences
PAMK <- pamk(data = post_nutitional_means$diff)
#plot(PAMK$crit)

# were gonna go with 3 groups, slightly lower than optimum, crit score.
PAM <- pam(post_nutitional_means$diff, k = 3)

# add clusters to dataframe
PAMClust = rep("NA", length(post_nutitional_means$diff))

PAMClust[PAM$clustering == 1] = "High_changers"
PAMClust[PAM$clustering == 2] = "Medium_changers"
PAMClust[PAM$clustering == 3] = "Low_changers"

post_nutitional_means$Cluster = PAMClust
#Hmisc::describe(post_nutitional_means)

# now we combine with means and fiber groups
means_groups <- merge(Nutritional_means, post_nutitional_means, by.x = "Individual", by.y = "Individual")
metadata_new <- merge(metadata, means_groups, by.x = c("individual", "Intervention"), by.y = c("Individual", "Intervention.x"))
metadata_new_clean <- metadata_new %>% 
  select(X, individual, Intervention, day, Fiber.x, Protein.x, Carbohydrates.x, Fat.x, Calories.x, Calculated_calories.x, diff.x, Cluster)
names(metadata_new_clean) <- c("Metagenome","Individual", "Intervention", "Day", "Fiber", "Protein", "Carbohydrates", "Fat", "Calories", "Calculated_calories", "Diff_fiber", "Cluster")
metadata <- metadata_new_clean
#metadata %>%
#  group_by(Cluster) %>% 
#  summarise(Individuals = n_distinct(Individual), Samples = n_distinct(Metagenome), avg_fiber = mean(Fiber))

# These are the true differences
#post_nutitional_means %>%
#  group_by(Cluster) %>% 
#  summarise(Individuals = n_distinct(Individual), avg_fiber = mean(Fiber), avg_intervention_change = mean(diff))

# more average data (pre - post based on clusters)
#metadata %>% group_by(Cluster, Intervention) %>% summarise(avg_fiber = mean(Fiber))
# clean up:
clean_me <- c("PAM", "means_groups", "metadata_new", "metadata_new_clean", "PAMK", "Nutritional_means", "PAMClust")
rm(list = clean_me)


## Rarefy the midas data to 900 

lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)
library(EcolUtils)
library(vegan)
set.seed(seed = 999)
#barplot(sort(rowSums(midas)))
#abline(h = mean(rowSums(midas)), col = "Red")
#abline(h = median(rowSums(midas)), col = "Blue")
midas_alpha_rare <- rrarefy.perm(midas, sample = 900, n = 100, round.out = T)
alpha_rare <- midas_alpha_rare[rowSums(midas_alpha_rare) >= 900-(900*.1), 
                               colSums(midas_alpha_rare) >= 1]
#barplot(sort(rowSums(alpha_rare)))

lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)
