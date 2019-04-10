---
title: "Fiber_project"
output: html_document
---

```r
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(eval = FALSE)
```

## R Markdown

### Libraries to load:

```r
library(readxl)
library(ggplot2)
library(tidyverse)
library(vegan)
library(EcolUtils)
library(ggpubr)
library(reshape2)
library(gplots)
library(reticulate)
library(devtools)
library(cowplot)
library(nlme)
library(labdsv)
library(knitr)
library(kableExtra)
```

### Functions to load:

```r
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
```

### Python env to load:

```r
# use_python("/Library/Frameworks/Python.framework/Versions/3.5/bin/python3", required = T)
# py_config()
# pd <- import("pandas")
# sns <- import("seaborn")
```


## **NUTRITIONAL DATA**

```r
## Great package for starting a ggplot: esquisse
## esquisse::esquisser(yourdataframe)

## Load in the file
Nutritional_Data <- read_excel("Nutritional_Data.xlsx", skip = 1)

## Make the numeric individual a factor
Nutritional_Data$new_individ <- as.factor(Nutritional_Data$individual)

Nutritional_Data$`high-low fiber` <- factor(Nutritional_Data$`high-low fiber`,
                                            levels = c('low','high'),ordered = TRUE)

#average fiber intake between high and low treatments
group_by(Nutritional_Data, `high-low fiber`) %>% summarise(avgFiber = mean(fiber))

## FOR ALL THE SAMPLES (average across all of them)

A <- ggplot(data = Nutritional_Data) +
  aes(x = Nutritional_Data$`high-low fiber`, y = fiber, fill = `high-low fiber`) +
  geom_boxplot() +
  theme_classic((base_size = 18)) + scale_fill_brewer(palette="Set2") + ggtitle("Fiber") + xlab("") + ylab("Fiber (grams)") + scale_fill_manual(values=c("deepskyblue3", "orange")) + stat_compare_means(method = "t.test")

B <- ggplot(data = Nutritional_Data) +
  aes(x = Nutritional_Data$`high-low fiber`, y = protein, fill = `high-low fiber`) +
  geom_boxplot() +
  theme_classic((base_size = 18)) + scale_fill_brewer(palette="Set2") + ggtitle("Protein") + xlab("") + ylab("Protein (grams)") + scale_fill_manual(values=c("deepskyblue3", "orange")) #+ stat_compare_means(method = "t.test")

C <- ggplot(data = Nutritional_Data) +
  aes(x = Nutritional_Data$`high-low fiber`, y = carbs, fill = `high-low fiber`) +
  geom_boxplot() +
  theme_classic((base_size = 18)) + scale_fill_brewer(palette="Set2") + ggtitle("Carbohydrates") + xlab("") + ylab("Carbs (grams)") + scale_fill_manual(values=c("deepskyblue3", "orange")) #+ stat_compare_means(method = "t.test")

E <- ggplot(data = Nutritional_Data) +
  aes(x = Nutritional_Data$`high-low fiber`, y = fat, fill = `high-low fiber`) +
  geom_boxplot() +
  theme_classic((base_size = 18)) + scale_fill_brewer(palette="Set2") + ggtitle("Fat") + xlab("") + ylab("Fat (grams)") + scale_fill_manual(values=c("deepskyblue3", "orange")) #+ stat_compare_means(method = "t.test")

D <- ggplot(data = Nutritional_Data) +
  aes(x = Nutritional_Data$`high-low fiber`, y = calories, fill = `high-low fiber`) +
  geom_boxplot() +
  theme_classic((base_size = 18)) + scale_fill_brewer(palette="Set2") + ggtitle("Calories") + xlab("") + ylab("Calories (kcal)") + scale_fill_manual(values=c("deepskyblue3", "orange")) #+ stat_compare_means(method = "t.test")

plot_grid(A, B, C, E, D, ncol = 2, labels = "AUTO")
#multiplot(A, B, C, E, D, cols=2)

##############################################################################################

# for each individual

F <- ggplot(data = subset(Nutritional_Data, subset = Nutritional_Data$new_individ == 5)) +
  aes(x = `high-low fiber`, y = fiber, fill = `high-low fiber`) +
  geom_boxplot() +
  theme_classic((base_size = 16)) + scale_fill_brewer(palette="Set2") + ggtitle("Fiber") + xlab("") + ylab("Fiber (grams)") 

G <- ggplot(data = subset(Nutritional_Data, subset = Nutritional_Data$new_individ == 5)) +
  aes(x = `high-low fiber`, y = protein, fill = `high-low fiber`) +
  geom_boxplot() +
  theme_classic((base_size = 16)) + scale_fill_brewer(palette="Set2") + ggtitle("Protein") + xlab("") + ylab("Protein (grams)")

H <- ggplot(data = subset(Nutritional_Data, subset = Nutritional_Data$new_individ == 5)) +
  aes(x = `high-low fiber`, y = carbs, fill = `high-low fiber`) +
  geom_boxplot() +
  theme_classic((base_size = 16)) + scale_fill_brewer(palette="Set2") + ggtitle("Carbohydrates") + xlab("") + ylab("Carbs (grams)")

I <- ggplot(data = subset(Nutritional_Data, subset = Nutritional_Data$new_individ == 5)) +
  aes(x = `high-low fiber`, y = fat, fill = `high-low fiber`) +
  geom_boxplot() +
  theme_classic((base_size = 16)) + scale_fill_brewer(palette="Set2") + ggtitle("Fat") + xlab("") + ylab("Fat (grams)")

J <- ggplot(data = subset(Nutritional_Data, subset = Nutritional_Data$new_individ == 5)) +
  aes(x = `high-low fiber`, y = calories, fill = `high-low fiber`) +
  geom_boxplot() +
  theme_classic((base_size = 16)) + scale_fill_brewer(palette="Set2") + ggtitle("Calories") + xlab("") + ylab("Calories (kcal)")


```


## **SIMPER TEST**

```r
## Load in the raw data
simper <- read_excel("SIMPER_high_low.xlsx")

## Melt the data
melted_simper <- melt(simper)

## Subset the melted data to only have the bacteria that you want
new_melted_simper <- subset(melted_simper, 
                            melted_simper$Species == "Bifidobacterium_adolescentis_56815" 
                            | melted_simper$Species == "Eubacterium_rectale_56927" 
                            | melted_simper$Species == "Blautia_wexlerae_56130")

## Subset the new melted data to only have the variables you care about
clean_simper <- new_melted_simper[ which( new_melted_simper$variable == "High_avg_abund" 
                                          | new_melted_simper$variable == "Low_avg_abund") , ]

## Order the variables
clean_simper$variable <- factor(clean_simper$variable, 
                                levels = c('Low_avg_abund','High_avg_abund'),ordered = TRUE)

## Plot the data
ggplot(data = clean_simper) +
  aes(x = Species, fill = variable, weight = value) +
  geom_bar(position = 'dodge') +
  theme_classic((base_size = 18)) + scale_fill_brewer(palette="Set2") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6)) + 
  ylab("Average Abundance")
```

## **Rarefaction**

```r

## Open up your OTUtable. Make sure the samples are rows and the features are columns.
## Transpose if necessary. IF its a tsv, make sure thats right too.

## Bring in the data and check the lowest read count
midas <- t(read.csv("count_reads.txt", 
                    row.names = 1, check.names = FALSE, sep = "\t"))

barplot(sort(rowSums(midas)), ylim = c(0, max(rowSums(midas))),
        xlim = c(0,NROW(midas)), col = "Orange")
sort(rowSums(midas))

## Alpha rarefaction to 980

midas_alpha_rare <- rrarefy.perm(midas, sample = 980, n = 10, round.out = T)

alpha_rare <- midas_alpha_rare[rowSums(midas_alpha_rare) >= 980-(980*.1), colSums(midas_alpha_rare) >= 1]
sort(rowSums(alpha_rare))

barplot(sort(rowSums(alpha_rare)), ylim = c(0, max(rowSums(alpha_rare))), 
        xlim = c(0,NROW(alpha_rare)), col = "Blue")

## Beta rarefaction to 980

midas_beta_rare <- avgdist(midas, sample = 980, iterations = 100, 
                                              meanfun = median, dmethod = "bray")

bray_distance_midas <- as.data.frame(as.matrix(midas_beta_rare))

## write to files
setwd("~/Google Drive File Stream/My Drive/Class/M130L Fiber Study/Metagenomic_seq/midas/midas_take2/")
write.csv(bray_distance_midas, file = "bray_distances.csv")
write.csv(alpha_rare, file = "rarified_otu-table.csv")
```

## Metaphlan

```r

# condense merged_metaphlan.txt into genus
#library(shiny)
#runGitHub("taxonomy_solution",username = "swandro")

metaphlan <- read.csv("Condensed_metaphlan_table.tsv", row.names = 1, check.names = FALSE, sep = "\t")

fungal <- read.csv("fungal_metaphlan.txt", row.names = 1, check.names = FALSE, sep = "\t")

viruses <- read.csv("virus_metaphlan.txt", row.names = 1, check.names = FALSE, sep = "\t")
```



## Alpha Diverisity
```r
h <- diversity(alpha_rare, "shannon") 
richness <- specnumber(alpha_rare)
Evenness <- h/log(richness)

comm_meta <- merge(comm_meta, as.data.frame(h), by.x = "Row.names", by.y = "row.names")
comm_meta <- merge(comm_meta, as.data.frame(richness), by.x = "Row.names", by.y = "row.names")
comm_meta <- merge(comm_meta, as.data.frame(Evenness), by.x = "Row.names", by.y = "row.names")

comm_meta$intervention <- factor(comm_meta$`high-low fiber`, levels = c('low','high'),ordered = TRUE)

shannon <- ggplot(data = comm_meta) +
  aes(x = comm_meta$intervention, y = comm_meta$h, fill = intervention) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2) +
  theme_classic((base_size = 14)) + ggtitle("Alpha Diversity") + xlab("Intervention") + ylab("Shannon Diversity") + scale_fill_manual(values=c("deepskyblue3", "orange")) + stat_compare_means(method = "t.test", label.x.npc = .8)

richness <- ggplot(data = comm_meta) +
  aes(x = comm_meta$intervention, y = comm_meta$richness, fill = intervention) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2) +
  theme_classic((base_size = 14)) + ggtitle("Species Richness") + xlab("Intervention") + ylab("Number of Species") + scale_fill_manual(values=c("deepskyblue3", "orange")) + stat_compare_means(method = "t.test", label.x.npc = .8)

evenness <- ggplot(data = comm_meta) +
  aes(x = comm_meta$intervention, y = comm_meta$Evenness, fill = intervention) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2) +
  theme_classic((base_size = 14)) + ggtitle("Evenness") + xlab("Intervention") + ylab("Evenness (J')") + scale_fill_manual(values=c("deepskyblue3", "orange")) + stat_compare_means(method = "t.test", label.x.npc = .6)

plot_grid(richness, shannon, evenness, ncol = 3, labels = "AUTO")

```


## **Relative abundance Species**

```r
#make sure you generate alpha_rare from above (rarefaction scripts)
dat <- as.data.frame(alpha_rare)
#make everything relative abundance
comm <- decostand(dat, method = "total")
# select top 8 features (basically the ones that show up the most). Why 8? 
# because thats how many colorblind palette colors there are.
comm_order <- comm[, order(-colMeans(comm))]
comm_collapse <- unite(comm_order, "Other", colnames(comm_order[8]):colnames(rev(comm_order)[1]), remove = TRUE)
comm_collapse$Other <- NULL
comm_collapse$low_abundant <- (1 - rowSums(comm_collapse))
bacteria <- as.list(colnames(comm_collapse))

#add metadata
Nutritional_Data <- read_excel("/Volumes/GoogleDrive/My Drive/Class/M130L Fiber Study/Metagenomic_seq/Primer/Nutritional_Data.xlsx", skip = 1)
metadata <- dplyr::select(Nutritional_Data, X__1, individual, day, `high-low fiber`, fiber, carbs, percent_change)
#merge with metadata
comm_meta <- merge(comm_collapse, metadata, by.x = "row.names", by.y = "X__1")

#melt the data
comm_melted <- melt(comm_meta, id=c("individual", "day", "high-low fiber", "fiber", "carbs", "percent_change"), measure.vars = c(bacteria))
comm_melted <- rename(comm_melted, Taxonomy = variable, highlow = `high-low fiber`)
#plot the data: Relative abundance!!!
cbPalette <- rev(c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))
#facet wrap (big box)
p <- ggplot(comm_melted, aes(day, y=value, fill = Taxonomy)) + 
  geom_bar(stat = "identity") + 
  facet_wrap(.~ individual, scales = "free_y") + ggtitle("Taxa Relative Abundance") + 
  xlab("Individual by day") + 
  ylab("Relative abundance") +
  theme_minimal() + 
  scale_fill_manual(values=cbPalette)

p

# Scatter plot of the biggest changers vs percent change in fiber
ggscatter(comm_melted, x = "fiber", y = "value",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson") + 
          facet_wrap(.~ Taxonomy, scales = "free_y") + 
          ylab("Relative abundance") +
          theme_minimal(base_size = 16) + geom_jitter() +
          xlab("Fiber (grams)")

#plot the data: Scatterplot with data!!!!!


#facet grid (long box)
ggplot(comm_melted, aes(day, y=value, fill = Taxonomy, na.rm = TRUE)) + 
    geom_bar(stat = "identity") + 
    facet_grid(~ individual, ) + ggtitle("Taxa Relative Abundance") + 
    xlab("Individual by day") + 
    ylab("Relative abundance") + theme_minimal((base_size = 20)) + scale_fill_manual(values=cbPalette)

# LME on comm_meta
comm_meta$intervention <- comm_meta$`high-low fiber`

Ba.lme <- lme(Bifidobacterium_adolescentis_56815 ~ intervention, data = comm_meta, 
                        random = ~ 1|individual, cor=corAR1())
summary(Ba.lme)
anova(Eu.lme) %>% kable(caption = "LME Health Status", booktabs = T) %>% kable_styling(position = "float_right", font_size = 16, latex_options = c("striped", "hold_position"))

```





## **SCFA Plotting**

```r
SCFA <- read_xlsx("SCFA_concentrations.xlsx")
SCFA$highlow <- ordered(SCFA$`pre-post`, levels = c("pre","post"))
compare_means(norm_butyrate ~ highlow, data = SCFA, method = "t.test")

#esquisse::esquisser(SCFA)
acetate <- ggplot(data = SCFA) +
       aes(x = highlow, y = norm_acetate, fill = highlow) +
       geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2) +
       labs(x = 'Fiber Intervention',
                       y = 'Acetate (ng/ul)') +
       theme_classic(base_size = 14)
A <- acetate + scale_fill_manual(values=c("deepskyblue3", "orange")) + stat_compare_means(method = "t.test")
propionate <- ggplot(data = SCFA) +
         aes(x = highlow, y = norm_propionate, fill = highlow) +
         geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2) +
         labs(x = 'Fiber Intervention', y = 'Propionate (ng/ul)') +
         theme_classic(base_size = 14)
P <- propionate + scale_fill_manual(values=c("deepskyblue3", "orange")) + stat_compare_means(method = "t.test")
isob <- ggplot(data = SCFA) +
         aes(x = highlow, y = norm_isobutyrate, fill = highlow) +
         geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2) +
         labs(x = 'Fiber Intervention',
                         y = 'Isobutyrate (ng/ul)') +
         theme_classic(base_size = 14)
IB <- isob + scale_fill_manual(values=c("deepskyblue3", "orange")) + stat_compare_means(method = "t.test")
but <- ggplot(data = SCFA) + 
         aes(x = highlow, y = norm_butyrate, fill = highlow) +
         geom_boxplot(outlier.shape = NA) +geom_jitter(width = 0.2) +
         labs(x = 'Fiber Intervention',
                         y = 'Butyrate (ng/ul)') +
         theme_classic(base_size = 14)
B <- but + scale_fill_manual(values=c("deepskyblue3", "orange")) + stat_compare_means(method = "t.test")
isoval <- ggplot(data = SCFA) + 
         aes(x = highlow, y = norm_isovalerate, fill = highlow) +
         geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2) +
         labs(x = 'Fiber Intervention',
                         y = 'Isovalerate (ng/ul)') +
         theme_classic(base_size = 14)
IV <- isoval + scale_fill_manual(values=c("deepskyblue3", "orange")) + stat_compare_means(method = "t.test")
val <- ggplot(data = SCFA) +
         aes(x = highlow, y = norm_valerate, fill = highlow) +
         geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2) +
         labs(x = 'Fiber Intervention',
                         y = 'Valerate (ng/ul)') +
         theme_classic(base_size = 14)
V <- val + scale_fill_manual(values=c("deepskyblue3", "orange")) + stat_compare_means(method = "t.test")

pdf("SCFA_all.pdf")
plot_grid(A, P, IB, B, IV, V, ncol = 2, labels = "AUTO")
dev.off()
```




## **PRES/ABSENCE**

```r
# So to get the genes_presabs.txt, you have to run a few things on the cluster using midas.
# First, in the directory: /bio/aoliver2/fiber_metagenome/quality_filtered/midas/midas_take2/ 
# i ran this script:

# /bio/aoliver2/software/MIDAS-1.3.2/scripts/merge_midas.py genes merged_genes/ 
# -i /bio/aoliver2/fiber_metagenome/quality_filtered/midas/midas_take2/ -t dir --min_samples 10
# (the above only worked for 1 genome so far, Eubacterium_rectale)

# then i ran this:
#compare_genes.py  merged_genes/Eubacterium_rectale_56927/ --out merged_genes/Eubacterium_rectale.distances.txt

# Not sure what that distances.txt file is useful for yet, but in the directory merged_genes/Eubacterium_rectale_56927/
# youll find genes_presabs.txt, which is input for below:

core_pan <- read.csv("genes_presabs_E-rectale.txt", sep = '\t', check.names = F)
rownames(core_pan) <- core_pan[,1]
core_pan[,1] <- NULL
core_pan <- core_pan[order(rowSums(core_pan)),]
#write.csv(core_pan, file = "core_pan_ordered.txt", sep = "\t")
core_pan.m <- as.matrix(core_pan)
data.dist <- vegdist(t(core_pan.m), method = "euclidean")
row.clus <- hclust(data.dist, "aver")

 pdf("pres_abs_blau_w")
 heatmap.2(t(core_pan.m),
           Rowv = as.dendrogram(row.clus),
           Colv = F,
           dendrogram = "row",
           #scale = "row",
           col = c("navy", "white"),
           #density.info = "none",
           trace = "none",
           #srtCol = -45,
           #adjCol = c(.1, .5),
           xlab = "Identifier",
           ylab = "Rows",
           rowsep = c(4,7,11,16,20,22,27,30,34,38,44),
           key = FALSE
 )
 dev.off()
```

## **LME on PC**

```r
# LME on PC-1 
tmp.pco <- pco(bray_distance_midas, k =2)
metadata <- merge(as.data.frame(tmp.pco$points), Nutritional_Data, by.x = "row.names", by.y = "X__1")
names(metadata)[2] <- "PCO1"
names(metadata)[3] <- "PCO2"
metadata$intervention <- metadata$`high-low fiber`

fiber.lme <- lme(PCO1 ~ intervention, data = metadata, 
                        random = ~ 1|individual, cor=corAR1())
summary(fiber.lme)
anova(fiber.lme) %>% kable(caption = "LME Intervention", booktabs = T) %>% kable_styling(position = "float_right", font_size = 16, latex_options = c("striped", "hold_position"))
```

## PERMANOVA

```r
permanova.input <- merge(metaphlan, Nutritional_Data, by.x= "row.names", by.y= "X__1")
permanova.input.midas <- merge(alpha_rare, Nutritional_Data, by.x= "row.names", by.y= "X__1")

intervention.perm <- adonis(permanova.input[,2:13] ~ permanova.input$`high-low fiber`, permutations = 9999, method = "bray", by = "terms", strata = as.character(permanova.input$individual)) 

individual.perm <- adonis(permanova.input[,2:13] ~ as.character(permanova.input$individual), permutations = 9999, method = "bray", by = "terms")

## THIS IS THE ONE!! DONT TOUCH. COPY.########################
intervention.perm.midas <- adonis(permanova.input.midas[,2:367] ~ as.character(permanova.input.midas$individual)/permanova.input.midas$`high-low fiber`, permutations = 9999, method = "bray")
##########################################################################

coef <- coefficients(intervention.perm)["permanova.input$`high-low fiber`1",]
top.coef <- coef[rev(order(abs(coef)))[1:20]]
barplot(sort(top.coef), horiz = T, las = 1, main = "Top taxa", col = "deepskyblue3", cex.axis=1.5, cex.names=1.5)

anova(betadisper(d = as.dist(bray_distance_midas), permanova.input$`high-low fiber`))


```


## Beta Diveristy 

```r

fiber.mds <- metaMDS(midas_alpha_rare, distance="bray", k=2)
stressplot(fiber.mds)

sites <- as.data.frame(scores(fiber.mds, display = "sites"))
species <- as.data.frame(scores(fiber.mds, display = "species"))

nmds.sites <- merge(sites, comm_meta, by.x = "row.names", by.y = "Row.names")
nmds.sites <- nmds.sites[order(nmds.sites$day),]

ndms.sites.subset <- subset(nmds.sites, nmds.sites$individual=="1")
ndms.sites.subset <- ndms.sites.subset[order(ndms.sites.subset$day),]

nmds.plot <- ggplot() + 
  geom_point(data = nmds.sites, aes(NMDS1, NMDS2, fill = as.factor(individual), size = fiber, alpha = 0.8), colour="black",pch=21) + 
  geom_path(data=nmds.sites, aes(x=NMDS1,y=NMDS2,group=as.factor(individual),colour=as.factor(individual)), linetype = "dashed",size=0.8, arrow = arrow(), show.legend = FALSE) + 
  theme_classic() + scale_size(range = c(0.5, 10)) + 
  annotate("text", label = "Stress = 0.1546197\nk = 3", x = 1.2, y = 0.7) +
  theme(axis.title = element_text(size = 20), axis.text = element_text(size = 18), 
        legend.text = element_text(size = 20), legend.title = element_text(size = 20)) +
  geom_text(data = nmds.sites, aes(NMDS1, NMDS2, label = individual), check_overlap = TRUE)

nmds.plot.subset <- ggplot() + 
  geom_point(data = ndms.sites.subset, aes(NMDS1, NMDS2, fill = as.factor(individual), size = fiber, alpha = 0.8), colour="black",pch=21) + 
  geom_path(data=ndms.sites.subset, aes(x=NMDS1,y=NMDS2,group=as.factor(individual),colour=as.factor(individual)), linetype = "dashed",size=0.8, arrow = arrow(), show.legend = FALSE) + 
  theme_classic() + scale_size(range = c(0.5, 10)) + 
  annotate("text", label = "Stress = 0.1546197\nk = 3", x = 1.2, y = 0.7) +
  theme(axis.title = element_text(size = 20), axis.text = element_text(size = 18), 
        legend.text = element_text(size = 20), legend.title = element_text(size = 20)) +
  geom_text(data = ndms.sites.subset, aes(NMDS1, NMDS2, label = individual), check_overlap = TRUE)

ggsave(filename = nmds_plot.pdf, plot = nmds.plot)

```


### 16S analysis

```r
## Open up your OTUtable. Make sure the samples are rows and the features are columns.
## Transpose if necessary. IF its a tsv, make sure thats right too.

## Bring in the data and check the lowest read count
OTU_raw <- read.csv("OTU_table.csv", 
                    row.names = 1, check.names = FALSE, sep = ",")

barplot(sort(rowSums(OTU_raw)), ylim = c(0, max(rowSums(OTU_raw))),
        xlim = c(0,NROW(OTU_raw)), col = "Orange")
sort(rowSums(OTU_raw))

## Alpha rarefaction to 980

OTU_alpha_rare <- rrarefy.perm(OTU_raw, sample = 12500, n = 10, round.out = T)

OTU_alpha_rare <- OTU_alpha_rare[rowSums(OTU_alpha_rare) >= 12500-(12500*.1), colSums(OTU_alpha_rare) >= 1]
sort(rowSums(OTU_alpha_rare))

barplot(sort(rowSums(OTU_alpha_rare)), ylim = c(0, max(rowSums(OTU_alpha_rare))), 
        xlim = c(0,NROW(OTU_alpha_rare)), col = "Blue")

## Beta rarefaction to 980

OTU_beta_rare <- avgdist(OTU_raw, sample = 12500, iterations = 100, 
                                              meanfun = median, dmethod = "bray")

bray_distance_OTU <- as.data.frame(as.matrix(OTU_beta_rare))

## write to files
#setwd("~/Google Drive File Stream/My Drive/Class/M130L Fiber Study/Metagenomic_seq/midas/midas_take2/")
#write.csv(bray_distance_midas, file = "bray_distances.csv")
#write.csv(alpha_rare, file = "rarified_otu-table.csv")

otu_with_metadata <- merge(OTU_alpha_rare, Nutritional_Data, by.x = "row.names", by.y = "uniq-id")
row.names(otu_with_metadata) <- otu_with_metadata$Row.names
otu_with_metadata$Row.names <- NULL

otu_shannon <- diversity(otu_with_metadata[,1:2035], "shannon") 
otu_richness <- specnumber(otu_with_metadata[,1:2035])
otu_Evenness <- otu_shannon/log(otu_richness)
h_otu <- as.data.frame(otu_shannon)
otu_with_metadata <- merge(otu_with_metadata, h_otu, by.x = "row.names", by.y = "row.names")
otu_with_metadata <- merge(otu_with_metadata, as.data.frame(otu_richness), by.x = "row.names", by.y = "Row.names")
otu_with_metadata <- merge(otu_with_metadata, as.data.frame(otu_Evenness), by.x = "row.names", by.y = "row.names")


shannon_otu <- ggplot(data = otu_with_metadata) +
  aes(x = otu_with_metadata$`high-low fiber`, y = otu_with_metadata$otu_shannon, fill = otu_with_metadata$`high-low fiber`) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2) +
  theme_classic((base_size = 14)) + ggtitle("Alpha Diversity") + xlab("Intervention") + ylab("Shannon Diversity") + scale_fill_manual(values=c("deepskyblue3", "orange")) + stat_compare_means(method = "t.test", label.x.npc = .8)

richness_otu <- ggplot(data = otu_with_metadata) +
  aes(x = comm_meta$intervention, y = comm_meta$richness, fill = intervention) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2) +
  theme_classic((base_size = 14)) + ggtitle("Species Richness") + xlab("Intervention") + ylab("Number of Species") + scale_fill_manual(values=c("deepskyblue3", "orange")) + stat_compare_means(method = "t.test", label.x.npc = .8)

evenness_otu <- ggplot(data = otu_with_metadata) +
  aes(x = otu_with_metadata$intervention, y = otu_with_metadata$Evenness, fill = intervention) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2) +
  theme_classic((base_size = 14)) + ggtitle("Evenness") + xlab("Intervention") + ylab("Evenness (J')") + scale_fill_manual(values=c("deepskyblue3", "orange")) + stat_compare_means(method = "t.test", label.x.npc = .6)

plot_grid(richness, shannon, evenness, ncol = 3, labels = "AUTO")

```

##FUNGAL Analysis

```r
fungal_counts <- read.csv("all_fungal_counts.txt", check.names = FALSE, sep = "\t")

fungal_counts$intervention <- factor(fungal_counts$intervention, levels = c('low','high'), ordered = TRUE)

library(esquisse)
esquisse::esquisser(fungal_counts)

ggplot(data = fungal_counts) +
  aes(x = species, fill = intervention, weight = log2(normalized_reads), na.rm = TRUE) +
  geom_bar(position = "dodge") +
  labs(y = 'log2(Normalized Read Counts)') +
  theme_classic() +
  facet_wrap(vars(individual)) + scale_fill_manual(values=c("deepskyblue3", "orange")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

```