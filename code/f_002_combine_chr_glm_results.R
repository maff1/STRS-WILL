#!/usr/bin/env Rscript
################################################################################
# COMBINE GLM RESULTS ## V001 ## MAFF ##########################################

rm(list=ls())
library(dplyr)

lsRes <- list.files("./results", pattern = ".rds", full.names = TRUE)
names(lsRes) <- gsub("glm_|_pruned_UKB.rds", "", basename(lsRes))
df <- do.call(rbind, lapply(lsRes, readRDS))
df["chr"] <- unlist(lapply(
  strsplit(rownames(df), split = "\\."), function(x) x[1])
)
# p.value adjusted
df <- df %>%
  group_by(chr) %>% 
  mutate(p.value.adj = p.adjust(pvalue, method='BH'))
# ---------------------------------------------------------------------------- #

df <- data.frame(chr = df$chr,
                 strName = df$strName,
                 beta = df$beta,
                 se = df$se,
                 t.value = df$t.value,
                 p.value = df$pvalue,
                 p.value.adj = df$p.value.adj,
                 aic = df$aic,
                 bic = df$bic,
                 r2 = df$r2,
                 sampleSize = df$sample.size,
                 row.names = NULL
)
saveRDS(df, "./results/combined_chr_glm.rds")
# ---------------------------------------------------------------------------- #