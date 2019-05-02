library(readxl)
library(ggplot2)
library(ggpubr)
library(cowplot)

setwd("~/Google Drive File Stream/My Drive/Github/Fiber-Analysis/")
SCFA <- read_xlsx("SCFA_concentrations.xlsx")
SCFA$highlow <- ordered(SCFA$`pre-post`, levels = c("pre","post"))

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

## TECHNICAL VARIATION

setwd("~/Google Drive File Stream/My Drive/Github/Fiber-Analysis/")
tech_var <- read_xlsx("SCFA_concentrations.xlsx",sheet = 2)
tech_var$highlow <- ordered(tech_var$`pre-post`, levels = c("pre","post"))

dacetate <- ggplot(data = tech_var) +
  aes(x = highlow, y = delta_acetate, fill = highlow) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2) +
  labs(x = 'Fiber Intervention',
       y = 'Acetate (ng/ul)') +
  theme_classic(base_size = 14)
dA <- dacetate + scale_fill_manual(values=c("deepskyblue3", "orange")) + stat_compare_means(method = "t.test")
dpropionate <- ggplot(data = tech_var) +
  aes(x = highlow, y = delta_propionate, fill = highlow) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2) +
  labs(x = 'Fiber Intervention', y = 'Propionate (ng/ul)') +
  theme_classic(base_size = 14)
dP <- dpropionate + scale_fill_manual(values=c("deepskyblue3", "orange")) + stat_compare_means(method = "t.test")
disob <- ggplot(data = tech_var) +
  aes(x = highlow, y = delta_isobutyrate, fill = highlow) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2) +
  labs(x = 'Fiber Intervention',
       y = 'Isobutyrate (ng/ul)') +
  theme_classic(base_size = 14)
dIB <- disob + scale_fill_manual(values=c("deepskyblue3", "orange")) + stat_compare_means(method = "t.test")
dbut <- ggplot(data = tech_var) + 
  aes(x = highlow, y = delta_butyrate, fill = highlow) +
  geom_boxplot(outlier.shape = NA) +geom_jitter(width = 0.2) +
  labs(x = 'Fiber Intervention',
       y = 'Butyrate (ng/ul)') +
  theme_classic(base_size = 14)
dB <- dbut + scale_fill_manual(values=c("deepskyblue3", "orange")) + stat_compare_means(method = "t.test")
disoval <- ggplot(data = tech_var) + 
  aes(x = highlow, y = delta_isovalerate, fill = highlow) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2) +
  labs(x = 'Fiber Intervention',
       y = 'Isovalerate (ng/ul)') +
  theme_classic(base_size = 14)
dIV <- disoval + scale_fill_manual(values=c("deepskyblue3", "orange")) + stat_compare_means(method = "t.test")
dval <- ggplot(data = tech_var) +
  aes(x = highlow, y = delta_valerate, fill = highlow) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2) +
  labs(x = 'Fiber Intervention',
       y = 'Valerate (ng/ul)') +
  theme_classic(base_size = 14)
dV <- dval + scale_fill_manual(values=c("deepskyblue3", "orange")) + stat_compare_means(method = "t.test")

pdf("SCFA_all.pdf")
plot_grid(dA, dP, dIB, dB, dIV, dV, ncol = 2, labels = "AUTO")
dev.off()
