### Global ---------------------------------------------------------------------
source("Data_prep.R")

library(glmmTMB)
library(ggplot2)
library(visreg)
library(corrplot)
library(factoextra)
library(dplyr)
library(scales)
library(glmmTMB)
library(car)
library(performance)

library(MuMIn)
options(na.action = na.omit)

### Descriptive statistics -----------------------------------------------------
## Correlation 
COR_minerals <- round(cor(df_out[,c(11:21)]), 2) # coarse.silt.proportion seems to be a good predictor, as highly correlated with the rest 
COR_explanatory <- round(cor(df_out[,c(6,7,9,10,17,23,24)]), 2)

cor_by_glacier <- by(df_out[,c(6,7,9,10,17,23,24)], df_out$glacier, cor)

ggplot(df_in, aes(log(Pb), log(Cs))) + 
  geom_smooth(method = "lm", aes(log(Pb), log(Cs)), color = "bisque4", se = FALSE, linetype = 8) +
  geom_smooth(method = "lm", aes(color = Glacier), alpha = 0.1) + 
  geom_point(aes(color = Glacier), size = 2) + theme_classic(base_size = 13) + scale_color_brewer(palette = "Dark2") + 
  ylab(expression(Cs-137~(log~activity~concentration))) + theme(plot.title = element_text(hjust = 0.5)) +
  xlab(expression(Pb-210~(log~activity~concentration))) + ggtitle("Relationship of Cs-137 and Pb-210")
ggsave("Output/Paper_1_Cs_Pb.png", height = 12, width = 16, dpi = 300, units = "cm")

## means and SEs
# Cs and Pb 
means_gl <- aggregate(Pb ~ Glacier, FUN = mean, data = df_in) 
means_gl <- cbind(means_gl, aggregate(Cs ~ Glacier, FUN = mean, data = df_in)[,2])
names(means_gl)[3] <- "Cs"

sds_Pb <- tapply(df_in$Pb, df_in$Glacier, sd)
sds_Cs <- tapply(df_in$Cs, df_in$Glacier, sd)

len <- 10
means_gl$Pb_SE <- sds_Pb/sqrt(len)
means_gl$Cs_SE <- sds_Cs/sqrt(len)

# Pu isotopes
Pu_raw <- readxl::read_xlsx("Input/alfa_sediment.xlsx")

means_gl <- cbind(means_gl, aggregate(Pu238 ~ Glacier, FUN = mean, data = Pu_raw)[,2])
means_gl <- cbind(means_gl, aggregate(Pu239_240 ~ Glacier, FUN = mean, data = Pu_raw)[,2])
names(means_gl)[6] <- names(Pu_raw)[4]
names(means_gl)[7] <- names(Pu_raw)[5]

sds_Pu238 <- tapply(Pu_raw$Pu238, Pu_raw$Glacier, sd)
sds_Pu239_240 <- tapply(Pu_raw$Pu239_240, Pu_raw$Glacier, sd)

len_pu <- 3 
means_gl$Pu238_SE <- sds_Pu238/sqrt(len)
means_gl$Pu23_240_SE <- sds_Pu239_240/sqrt(len)

means_gl <- cbind(means_gl, aggregate(Ratio_238to239_240 ~ Glacier, FUN = mean, data = Pu_raw)[,2])
names(means_gl)[10] <- "ratio_238to239_240"

sds_ratioPu <- tapply(Pu_raw$Ratio_238to239_240, Pu_raw$Glacier, sd)
means_gl$PuRatio_SE <- sds_ratioPu/sqrt(len)

means_gl <- means_gl %>% 
  mutate_if(is.numeric, round, digits = 2)

write.csv(means_gl, "Output/summary_statistics_radioactivity.csv")

### Raw data visualisation -----------------------------------------------------
## Correlation
corrplot(COR_minerals, type = "upper", order = "original", 
         tl.col = "black", tl.srt = 45)
corrplot(COR_explanatory, type = "upper", order = "original", 
         tl.col = "black", tl.srt = 45)

png("Output/cor_1_Blanc.png", height = 12, width = 12, units = "cm", res = 200)
corrplot(cor_by_glacier[[1]], type = "upper", order = "original", 
         tl.col = "black", tl.srt = 45)# Blanc
dev.off()

png("Output/cor_2_Ferpecle.png", height = 12, width = 12, units = "cm", res = 200)
corrplot(cor_by_glacier[[2]], type = "upper", order = "original", 
         tl.col = "black", tl.srt = 45)# Ferpecle
dev.off()

png("Output/cor_3_Pasterze.png", height = 12, width = 12, units = "cm", res = 200)
corrplot(cor_by_glacier[[3]], type = "upper", order = "original", 
         tl.col = "black", tl.srt = 45)# Pasterze
dev.off()

png("Output/cor_4_Preda_Rossa.png", height = 12, width = 12, units = "cm", res = 200)
corrplot(cor_by_glacier[[4]], type = "upper", order = "original", 
         tl.col = "black", tl.srt = 45)# Preda rossa
dev.off()

png("Output/cor_5_Ventina.png", height = 12, width = 12, units = "cm", res = 200)
corrplot(cor_by_glacier[[5]], type = "upper", order = "original", 
         tl.col = "black", tl.srt = 45)# Ventina
dev.off()

## Biological components across glaciers
ggplot(df_out, aes(y = mean_volume_mm3, x = Glacier)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(color = "darkgreen") + theme_bw(base_size = 13) +
  ylab(expression(Volume~of~cryoconite~granules~(mm)^{3})) +
  ggtitle("Cryoconite granules") + theme(plot.title = element_text(hjust = 0.5)) 
ggsave("Output/Cryoconite_size_glaciers.png", height = 8, width = 12, dpi = 300, units = "cm")

ggplot(df_out, aes(y = Chlorophyll, x = Glacier)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(color = "darkgreen") + theme_bw(base_size = 13) +
  ylab(expression(Chlorophyll~concentration~(ug~g)^{-1})) +
  ggtitle("Chlorophyll concentration") + theme(plot.title = element_text(hjust = 0.5)) 
ggsave("Output/Chlorophyll_glaciers.png", height = 8, width = 12, dpi = 300, units = "cm")

ggplot(df_out, aes(y = OM, x = Glacier)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(color = "darkgreen") + theme_bw(base_size = 13) +
  ylab(expression(Organic~matter~content~("%"))) +
  ggtitle("Organic matter content") + theme(plot.title = element_text(hjust = 0.5)) 
ggsave("Output/OM_glaciers.png", height = 8, width = 12, dpi = 300, units = "cm")

ggplot(df_out, aes(y = Cyano_to_all_ratio, x = Glacier)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(color = "darkgreen") + theme_bw(base_size = 13) + 
  ylab(expression(Cyanobacteria~ratio~to~all~bacteria)) +
  ggtitle("Cyanobacteria share") + theme(plot.title = element_text(hjust = 0.5)) 
ggsave("Output/Cyano_bacteria_glaciers.png", height = 8, width = 12, dpi = 300, units = "cm")

ggplot(df_out, aes(y = Coarse.Silt.Proportion,  x = Glacier)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(color = "darkgreen") + theme_bw(base_size = 13) + 
  ylab(expression(Silt~proportion~("%"))) +
  ggtitle("Silt proportion") + theme(plot.title = element_text(hjust = 0.5)) 
ggsave("Output/minerals_size_glaciers.png", height = 12, width = 18, dpi = 300, units = "cm")

## Radioisotopes across glaciers
ggplot(df_out, aes(y = Pb, x = Glacier)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(color = "darkgreen") + theme_bw(base_size = 13) + 
  ylab(expression(Activity~concentration~(Bq~kg)^{-1})) +
  ggtitle("Pb-210") + theme(plot.title = element_text(hjust = 0.5)) 
ggsave("Output/Pb_across_glaciers.png", height = 12, width = 18, dpi = 300, units = "cm")

ggplot(df_out, aes(y = Cs, x = Glacier)) + 
  geom_boxplot() + geom_jitter(color = "darkgreen") + theme_bw(base_size = 13) + 
  ylab(expression(Activity~concentration~(Bq~kg)^{-1})) +
  ggtitle("Cs-137") + theme(plot.title = element_text(hjust = 0.5))
ggsave("Output/Cs_across_glaciers.png", height = 12, width = 18, dpi = 300, units = "cm")

### Mean-centered predictors ---------------------------------------------------
## Mean centering
df_out$Pb_LOG <- log(df_out$Pb)
df_out$Cs_LOG <- log(df_out$Cs)

df_pred_d <- df_out[,c(2, 3, 9, 10, 17, 23, 24)]

df_out_delta <- df_pred_d %>%
  group_by(glacier) %>%
  mutate(across(where(is.numeric), ~ . - mean(.)))

## Rescaling 
df_out_delta <- df_out_delta %>%
  group_by(glacier) %>%
  mutate(across(where(is.numeric), ~ rescale(.)))

df_out_delta <- cbind(df_out_delta, df_out[,c(6:8, 25, 26)])

### PCA ------------------------------------------------------------------------
## Radioisotopes 
df_pca_1 <- df_out[, c(2, 6:8)]

PCA_1 <- prcomp(df_pca_1[,-c(1,4)], center = TRUE, scale = TRUE)

PCA_1_res <- get_pca_var(PCA_1)
PCA_1_res$coord; PCA_1_res$contrib

get_eig(PCA_1)
fviz_screeplot(PCA_1, addlabels = TRUE, ylim = c(0, 50))

fviz_pca_var(PCA_1, col.var = "black")

fviz_pca_ind(PCA_1,
             label = "none",
             habillage = df_pca_1$glacier,
             palette = c("#00AFBB", "#E7B800", "#FC4E07", "Black", "Green"),
             addEllipses = TRUE, title = "bla")
## predictors
df_pca_2 <- df_out[, c(2, 9, 10, 23, 24)]

PCA_2 <- prcomp(df_pca_2[,-c(1)], center = TRUE, scale = TRUE)

PCA_2_res <- get_pca_var(PCA_2)
PCA_2_res$coord; PCA_2_res$contrib

get_eig(PCA_2)
fviz_screeplot(PCA_2, addlabels = TRUE, ylim = c(0, 50))

fviz_pca_var(PCA_2, col.var = "black")

fviz_pca_ind(PCA_2,
             label = "none",
             habillage = df_pca_2$glacier,
             palette = c("#00AFBB", "#E7B800", "#FC4E07", "Black", "Green"),
             addEllipses = TRUE) + 
  labs(title = "PCA biplot - variation between glaciers (PC1 and PC2)", color = "Glacier", fill = "Glacier", shape = "Glacier", 
       subtitle = "Predictors: OM, Chlorophyll, Cryoconite Volume, Cyanobacteria") + 
  theme_bw(base_size = 13) + theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5, size = 12)) + 
  scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2")
ggsave("Output/Paper_2_PCA.png", height = 12, width = 16, dpi = 600, units = "cm")

### Models - scaled and within subject mean-centered  data ---------------------
## 210Pb - log 
m1.f_pb <- glmmTMB(Pb_LOG ~ OM + mean_volume_mm3 + Chlorophyll + Cyano_to_all_ratio + Coarse.Silt.Proportion + 
                     Coarse.Silt.Proportion:OM + OM:Chlorophyll + OM:Cyano_to_all_ratio +
                     (1|glacier), 
                   data = df_out_delta); summary(m1.f_pb)

res_pb_1 <- dredge(m1.f_pb, trace = 2)
subset(res_pb_1, delta <= 2, recalc.weights = FALSE)
get.models(res_pb_1, subset = delta <= 2, REML = TRUE)

m1.s_pb <- glmmTMB(Pb_LOG ~ OM + Chlorophyll + (1|glacier), 
                   data = df_out_delta, REML = TRUE); summary(m1.s_pb)

check_collinearity(m1.s_pb)
check_convergence(m1.s_pb)

plot(residuals(m1.s_pb), fitted(m1.s_pb))
qqnorm(residuals(m1.s_pb)); qqline(residuals(m1.s_pb))

Anova(m1.s_pb)

# OM 
mod_pb_out_1 <- visreg(m1.s_pb, "OM", gg = TRUE, line.par = list(col = "darkblue"))
mod_pb_out_1 + theme_classic(base_size = 13) + geom_point(size = 3.5, col = "darkblue") +
  ggtitle("Pb-210 ~ organic matter, partial residuals") + theme(plot.title = element_text(hjust = 0.5)) +
xlab("Scaled organic matter content") + ylab(expression(Pb-210~(log~activity~concentration)))
ggsave("Output/Paper_3_Pb_OM.png", height = 12, width = 16, dpi = 600, units = "cm")

# Chlorophyll
mod_pb_out_3 <- visreg(m1.s_pb, "Chlorophyll", gg = TRUE, line.par = list(col = 'darkblue'))
mod_pb_out_3 + theme_classic(base_size = 13) + geom_point(size = 3.5, col = "darkblue") +
  ggtitle("Pb-210 ~ chlorophyll concentration, partial residuals") + theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Scaled chlorophyll concentration") + ylab(expression(Pb-210~(log~activity~concentration)))
ggsave("Output/Paper_4_Pb_Chlorophyll.png", height = 12, width = 16, dpi = 600, units = "cm")

## 137Cs - log 
m1.f_cs <- glmmTMB(Cs_LOG ~ OM + mean_volume_mm3 + Chlorophyll + Cyano_to_all_ratio + Coarse.Silt.Proportion + 
                     Coarse.Silt.Proportion:OM + OM:Chlorophyll + OM:Cyano_to_all_ratio +
                     (1|glacier), 
                   data = df_out_delta); summary(m1.f_cs)

res_cs_1 <- dredge(m1.f_cs, trace = 2)
subset(res_cs_1, delta <= 2, recalc.weights = FALSE)
get.models(res_cs_1, subset = delta < 2, REML = TRUE)

m1.s_cs <- glmmTMB(Cs_LOG ~  OM + (1|glacier), 
                   data = df_out_delta, REML = TRUE); summary(m1.s_cs)

check_convergence(m1.s_cs)

plot(residuals(m1.s_cs), fitted(m1.s_cs))
qqnorm(residuals(m1.s_cs)); qqline(residuals(m1.s_cs))

Anova(m1.s_cs)

# OM
mod_cs_out_1 <- visreg(m1.s_cs, "OM", gg = TRUE, line.par = list(col = 'darkgreen'))
mod_cs_out_1 + theme_classic(base_size = 13) + geom_point(size = 3.5, color = "darkgreen") +
  ggtitle("Cs-137 ~ organic matter, partial residuals") + theme(plot.title = element_text(hjust = 0.5)) + 
  xlab("Scaled organic matter content") + ylab(expression(Cs-137~(log~activity~concentration)))
ggsave("Output/Paper_5_Cs_OM.png", height = 12, width = 16, dpi = 600, units = "cm")