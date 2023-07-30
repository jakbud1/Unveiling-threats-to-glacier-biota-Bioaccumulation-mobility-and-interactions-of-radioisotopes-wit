### Global ---------------------------------------------------------------------
library(dplyr)
library(readr)
library(ggplot2)
library(tidyverse)

stock_names <- readxl::read_xlsx("Input/Names.xlsx")
stock_names$sample <- as.character(stock_names$sample)

### Main_data ------------------------------------------------------------------
df_in <- readxl::read_xlsx("Input/main_df.xlsx")

### Minerals size --------------------------------------------------------------
df_minerals <- readxl::read_xlsx("Input/qgrain_results.xlsx")
df_minerals <- data.frame(df_minerals[,c(1,14:24)])

### Granules -------------------------------------------------------------------
df_gran <- readxl::read_xlsx("Input/granules_size.xlsx")

out_gran <- as.data.frame(tapply(df_gran$Length_mm, list(df_gran$Sample_name), 
                                 mean, trim = 0.05))
names(out_gran) <- "mean_diameter"
out_gran$sample_name <- row.names(out_gran) 

out_gran$mean_volume_mm3 <- 4/3 * pi * ((out_gran$mean_diameter/2)^3); rm(df_gran)

### Cyanobacteria relative abundance -------------------------------------------
df_v4 <- read_delim("Output/taxa.silva_V4_all.csv", 
                                delim = ";", escape_double = FALSE, trim_ws = TRUE)

df_v4 <- df_v4[complete.cases(df_v4$Phylum), ]

controls <- c(44:49,73)
df_v4$All_controls <- rowSums(df_v4[,controls])
df_v4 <- df_v4[,-controls]; rm(controls)

df_v4_nc <- subset(df_v4, All_controls == 0 | Phylum == "Cyanobacteria")

dfmol_out <- aggregate(. ~ Phylum, df_v4_nc[,-c(1:2, 67,69:72)], FUN = sum, simplyfy = FALSE)
dfmol_out$Phylum <- as.character(dfmol_out$Phylum)

dfmol_out <- rbind(dfmol_out, c(NA, mapply(sum, dfmol_out[,-1])))

dfmol_out[30,1] <- "total"

dfmol_out <- as.data.frame(t(rbind(dfmol_out, c("Cyano_to_all_ratio", t(dfmol_out[10,-1]/dfmol_out[30,-1])))))

names(dfmol_out) <- t(unname(dfmol_out[1,]))
dfmol_out <- dfmol_out[-1,]

dfmol_out <- mutate_all(dfmol_out, function(x) as.numeric(as.character(x)))
dfmol_out <- dfmol_out[-nrow(dfmol_out),]

# Names 
dfmol_out$run <- ifelse(substring(rownames(dfmol_out), 13) == "1", 1, 2)
dfmol_out$sample <- substr(substring(rownames(dfmol_out), 1), 1,3)
rownames(dfmol_out) <- NULL

# Merge by runs 
df_rat_out <- aggregate(Cyano_to_all_ratio ~ sample, data = dfmol_out, FUN = mean)

# write full data
write.csv(dfmol_out, "Output/NGS_V4_phyla.csv")

### Merge all dfs --------------------------------------------------------------
df_list <- list(stock_names, df_in, df_minerals, out_gran)
df_out <- df_list %>% reduce(full_join, by='sample_name')

df_out <- merge(df_out, df_rat_out, by = "sample")

### Fill the NAs, it's not a lack of data, there were no granules at all ------- 
df_out[11:20,22:23] <- 0