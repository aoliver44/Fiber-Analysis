##################################
######### Beta Diversity #########
##################################

# path should be set to data directory
path <- getwd()
setwd(path)
source("../scripts/Generate_basic_env.R")

library(vegan)
library(labdsv)
library(ggplot2)
library(grid)
library(dplyr)
library(reshape2)
# Permanovas
set.seed(999)
beta_diversity_data <- merge(metadata, alpha_rare, by.y = "row.names", by.x = "Metagenome")
nested_permanova <- adonis2(beta_diversity_data[,13:NCOL(beta_diversity_data)] ~ as.character(Individual)*Intervention, data = beta_diversity_data, permutations = 999, parallel = 4, method = "bray", by = "margin")
nested_permanova

intervention_permanova <- adonis(beta_diversity_data[,13:NCOL(beta_diversity_data)] ~ Intervention, data = beta_diversity_data, permutations = 999, parallel = 4, method = "bray", by = NULL)
intervention_permanova

cluster_permanova <- adonis(beta_diversity_data[,13:NCOL(beta_diversity_data)] ~ Cluster/as.character(Individual), data = beta_diversity_data, permutations = 999, parallel = 4, method = "bray")
cluster_permanova

############## Coeffecients on L6 (Genus) #################################
out_melted <- read.csv("taxonomy_melted.tsv", sep = "\t", header = T)
Genus_melted <- out_melted %>% filter(., L6 != "") %>% group_by(L6, variable) %>% summarise(genus_sum = sum(value))
Genus_OTU <- dcast(Genus_melted, formula = variable ~ L6)
rownames(Genus_OTU) <- Genus_OTU$variable
Genus_OTU$variable <- NULL
Genus_OTU_merge <- merge(metadata, Genus_OTU, by.x = "Metagenome", by.y = "row.names")

genus_permanova <- adonis(Genus_OTU_merge[,13:NCOL(Genus_OTU_merge)] ~ Intervention, data = Genus_OTU_merge, permutations = 999, parallel = 4, method = "bray")

coef <- coefficients(genus_permanova)["Intervention.L",]
top.coef <- coef[rev(order(abs(coef)))[1:20]]

par(mar=c(3,9,3,2) + 0.1)
barplot(sort(top.coef), horiz = T, las = 1, main = "Top taxa", col = "deepskyblue3", cex.axis=1, cex.names=0.8, xlim = c(-25,25))
dev.off()
# double check with random forest
library(rfPermute)
L6_intervention <- Genus_OTU_merge %>% dplyr::select(3, (13:NCOL(Genus_OTU_merge)))
names(L6_intervention) <- gsub(pattern = " ", replacement = "_", x = names(L6_intervention))
L6_RFP <- rfPermute(as.factor(Intervention) ~., data = L6_intervention, proximity = TRUE, importance = TRUE, corr.bias = TRUE, mtry = 60, ntree = 901, num.cores = 6)

proximityPlot(L6_RFP)
var_imp <- varImpPlot(L6_RFP, type = 1, n.var = 20)
plotConfMat(L6_RFP)
heatmap <- impHeatmap(L6_RFP, alpha = 0.05, ranks = F, n = 20)
heatmap + theme(axis.text.x = element_text(angle = 90, hjust = 1))
plot(L6_RFP)

# Double check with raw data
plot_data_coef <- merge(Genus_melted, metadata, by.x = "variable", by.y = "Metagenome")

ggplot(data = subset(plot_data_coef, plot_data_coef$L6 == "Bifidobacterium"), 
       aes(x = Intervention, y = log2(genus_sum))) +
  geom_boxplot() + geom_point(position = position_jitterdodge(jitter.width = .1), aes(color = as.factor(Individual))) +
  guides(fill=guide_legend(title="Individual")) +
  theme_bw(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Bifidobacterium") +
  labs(y = 'Read Counts',x = 'Intervention') 

###############################################

# LMEs
bray_distance <- vegdist(alpha_rare, method = "bray")
pco <- pco(bray_distance, k =2)
beta_lme <- merge(as.data.frame(pco$points), metadata, by.x = "row.names", by.y = "Metagenome")
names(beta_lme)[2] <- "PCO1"
names(beta_lme)[3] <- "PCO2"
beta_lme$rank_pco <- rank(beta_lme$PCO1)
beta.lme <- lme(rank_pco ~ Intervention, data = beta_lme, 
                random = ~ 1|as.factor(Individual), correlation = corAR1())
anova(beta.lme, beta.lme.nocor)

# vis of MDS of the data
set.seed(seed = 999)
beta.mds <- metaMDS(alpha_rare, distance="bray", k=2)
stressplot(beta.mds)

sites <- as.data.frame(scores(beta.mds, display = "sites"))
species <- as.data.frame(scores(beta.mds, display = "species"))

nmds.sites <- merge(sites, metadata, by.x = "row.names", by.y = "Metagenome")
nmds.sites <- nmds.sites[order(nmds.sites$Day),]

colors_inset <- c("#863636", "#8ad747", "#6e41c8", "#dbcb57", "#cb4bc0", "#6bd183", "#562d6f", "#839243", "#6f79cf", "#cf843b", "#8ab6d6", "#d74b34", "#8dceb6", "#d04a76", "#3c613a", "#ca8abe", "#37242f", "#d4aa97", "#56697f", "#785837")
colors_subset <- c("#863636","#8ad747","#6e41c8","#dbcb57","#cb4bc0","#6bd183","#562d6f","#839243", "#6f79cf","#cf843b","#8ab6d6","#d74b34","#8dceb6","#d04a76","#3c613a","#ca8abe","#37242f","#56697f","#785837")
# remove individual 13
nmds.sites.new <- subset(nmds.sites, nmds.sites$Individual != "13")
all.nmds.plot <- ggplot() + 
  geom_point(data = nmds.sites, aes(NMDS1, NMDS2, fill = "grey"), colour="black",pch=21, alpha = 0.45) + 
  geom_path(data=nmds.sites, aes(x=NMDS1,y=NMDS2,group=as.factor(Individual),colour="grey"), linetype = "dashed",size=0.3, arrow = arrow(), show.legend = FALSE) + 
  theme_classic() + scale_size(range = c(0.2, 6)) + 
  theme(axis.title = element_text(size = 20), axis.text = element_text(size = 18), 
        legend.text = element_text(size = 20), legend.title = element_text(size = 20)) +
  #geom_text(data = nmds.sites, aes(NMDS1, NMDS2, label = Individual), check_overlap = TRUE, size = 2.5) + 
  theme(legend.position = "none", axis.title.x=element_blank(),
        axis.title.y=element_blank(), axis.text.x=element_blank(),
        axis.text.y=element_blank(), axis.ticks=element_blank()) #+
#scale_color_manual(values = colors_inset) + scale_fill_manual(values = colors_inset) 

subset.nmds.plot <- ggplot(data = nmds.sites.new, aes(NMDS1, NMDS2)) + 
  geom_point(aes(color = as.factor(Individual), size = 1, shape = as.factor(Intervention)), alpha = 0.7) + 
  geom_path(aes(group=as.factor(Individual),color=as.factor(Individual)), linetype = "dashed",size=0.6, arrow = arrow(), show.legend = FALSE) + 
  theme_classic() + scale_size(range = c(0.2, 6)) + 
  theme(axis.title = element_text(size = 20), axis.text = element_text(size = 18), 
        legend.text = element_text(size = 20), legend.title = element_text(size = 20)) +
  geom_text(data = nmds.sites.new, aes(label = Individual), check_overlap = TRUE, size = 2.5) + theme(legend.position = "none") +
  scale_color_manual(values = colors_subset) + scale_fill_manual(values = colors_subset) + coord_cartesian(ylim = c(-0.6,0.9))

subset.nmds.plot + annotation_custom(ggplotGrob(all.nmds.plot), xmin = -0.78, xmax = -0.35, 
                                     ymin = 0.49, ymax = 0.95)

# determine what is pushing sample 13 away from everyone
all.nmds.plot <- ggplot() + 
  geom_point(data = nmds.sites, aes(NMDS1, NMDS2, fill = as.factor(Individual)), colour="black",pch=21, alpha = 0.45) + 
  geom_path(data=nmds.sites, aes(x=NMDS1,y=NMDS2,group=as.factor(Individual),colour=as.factor(Individual)), linetype = "dashed",size=0.3, arrow = arrow(), show.legend = FALSE) + 
  geom_segment(data=species,aes(x=0,xend=NMDS1,y=0,yend=NMDS2),
               arrow = arrow(length = unit(0.5, "cm")),colour="grey",inherit_aes=FALSE) +
  geom_text(data=species,aes(x=NMDS1,y=NMDS2,label=rownames(species)),size=3, position=position_jitter(width=0.8,height=0.9)) +
  theme_classic() + scale_size(range = c(0.2, 6)) + 
  theme(axis.title = element_text(size = 20), axis.text = element_text(size = 18), 
        legend.text = element_text(size = 20), legend.title = element_text(size = 20)) +
  #geom_text(data = nmds.sites, aes(NMDS1, NMDS2, label = Individual), check_overlap = TRUE, size = 2.5) + 
  theme(legend.position = "none", axis.title.x=element_blank(),
        axis.title.y=element_blank(), axis.text.x=element_blank(),
        axis.text.y=element_blank(), axis.ticks=element_blank()) +
  scale_color_manual(values = colors) + scale_fill_manual(values = colors) + coord_cartesian(xlim = c(6,14))

# Plot raw abundances
ggplot(data = Fiber_rfp_data) +
  aes(x = as.factor(Individual), y = Fiber_rfp_data$Salmonella_enterica_58266) +
  geom_boxplot() + geom_point(position = "jitter", width = 0.2, aes(colour = as.factor(Fiber_rfp_data$Individual))) +
  theme_bw() + 
  #facet_grid(. ~ Intervention) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_color_manual(values = colors)
