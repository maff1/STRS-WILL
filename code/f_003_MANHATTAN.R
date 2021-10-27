#!/usr/bin/env Rscript
################################################################################
# MANHATTAN PLOT ## V001 ## MAFF ###############################################

library(ggplot2)

df <- readRDS("./results/combined_chr_glm.rds")
df_sub <- df[df$strName %in% sample(df$strName, size = 10000),]
pdf(file = "./results/MANHATTAN.pdf")
ggplot(df_sub, aes(beta, -log10(p.value))) +
  geom_point() +
  facet_grid(~ chr, scales = "free_x", switch = "x") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  theme_bw()
dev.off()
# ---------------------------------------------------------------------------- #