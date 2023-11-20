### Global ---------------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(glmmTMB)
library(car)

options(na.action = na.omit)

### Uptake ratio ---------------------------------------------------------------
df_up <- readxl::read_xlsx("Input/in_animals.xlsx")
df_up2 <- readxl::read_xlsx("Input/in_animals.xlsx", sheet = 2)

ggplot(df_up2, aes(y = Ratio, x = Isotope)) + 
  geom_boxplot(outlier.shape = NA, aes(fill = Isotope, alpha = 0.4), show.legend = FALSE) + 
  geom_jitter(aes(color = Glacier), size = 3.5) + 
  theme_classic(base_size = 13) + theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom") +
  ggtitle("Uptake ratio of radioisotopes in springtails") + ylab("Uptake ratio") + 
  scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "OrRd")
ggsave("Output/Paper_6_uptake.png", height = 12, width = 16, dpi = 300, units = "cm")

# test 
m_up.1 <- glmmTMB(Ratio ~ Isotope + (1|Glacier),
                  data = df_up2[-c(9:16),]); summary(m_up.1)
res_up1 <- residuals(m_up.1)
qqnorm(res_up1); qqline(res_up1)

plot(res_up1 ~ fitted(m_up.1))
hist(res_up1)
Anova(m_up.1, type = "II")


### Mobility -------------------------------------------------------------------
df_mob <- readxl::read_xlsx("Input/mobility.xlsx", sheet = 1)
df_mob$Form <- factor(df_mob$Form, levels = c("Exchangeable", "Complex", "Specifically adsorbed"))

df_mob2 <- readxl::read_xlsx("Input/mobility.xlsx", sheet = 2)

mean(subset(df_mob2, Isotope == "Cs-137")$Ratio)
sd(subset(df_mob2, Isotope == "Cs-137")$Ratio)/sqrt(5)

## Visualisation 
# loosely
ggplot(subset(df_mob), aes(x = Form, y = Activity, fill = Form)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(aes(shape = Glacier)) + 
  facet_wrap(.~ Isotope, scale = "free") + theme_classic(base_size = 13) + scale_fill_brewer(palette = "Greens", guide="none") + 
  ylab(expression(paste("Activity concentration  ", "(Bq kg" ^-1,")"))) +
  ggtitle("Differences between loosely bound fractions") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "right",
        axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))
ggsave("Output/Paper_8_loosely.png", height = 12, width = 16, dpi = 600, units = "cm")

# ratio 
ggplot(subset(df_mob2), aes(x = Isotope, y = Ratio, fill = Isotope)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(aes(shape = Glacier)) + 
  theme_classic(base_size = 13) + scale_fill_brewer(palette = "Pastel1", guide = "none") + 
  ylab("Ratio") + 
  ggtitle("The ratio of loosely bound to total activity") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "right")
ggsave("Output/Paper_9_ratio_loosely.png", height = 12, width = 16, dpi = 600, units = "cm")

# test loosely (within)
df_mob$log_activity <- log(df_mob$Activity + 0.001)
m_mob.1 <- glmmTMB(log_activity ~ Form * Isotope + (1|Glacier),
                   data = df_mob, REML = TRUE); summary(m_mob.1)

res_mob1 <- residuals(m_mob.1)
plot(res_mob1 ~ fitted(m_mob.1))
hist(res_mob1)
qqnorm(res_mob1); qqline(res_mob1)

Anova(m_mob.1, type = "III")
visreg(m_mob.1, "Form", "Isotope")

# test ratio (loosely to total)
m_mob.2 <- glmmTMB(Ratio ~ Isotope + (1|Glacier),
                   data = df_mob2); summary(m_mob.2)

res_mob2 <- residuals(m_mob.2)
plot(res_mob2 ~ fitted(m_mob.2))
hist(res_mob2)
qqnorm(res_mob2); qqline(res_mob2)

Anova(m_mob.2, type = "II")