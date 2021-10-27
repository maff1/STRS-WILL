#!/usr/bin/env Rscript
################################################################################
# GLM GAUSSIAN MODELS ## V001 ## MAFF ##########################################

rm(list=ls())
library(data.table)
library(parallel)
library(rsq)

# "/t1-data/project/psychdementia/shared/Repeats/WGS_MRI_21Oct2021/MRI_GWAS/chr1_pruned_UKB.txt"
# "/t1-data/project/psychdementia/shared/Repeats/Variables/Covariates_MRI_HIP_1stAttempt"
# ---------------------------------------------------------------------------- #

args <- commandArgs(trailingOnly=TRUE)
datRepeats = args[1]
datCovariates = args[2]
ncores = as.integer(args[3])
# ---------------------------------------------------------------------------- #

print(datRepeats)
print(datCovariates)
# ---------------------------------------------------------------------------- #

cov <- data.table::fread(datCovariates); colnames(cov)[1] <- "ID"
setkey(cov, ID)
chr <- data.table::fread(datRepeats, key = "ID")

# CONVERT TO NUMERIC OR INTEGER
colsToConvert <- colnames(chr)[2:ncol(chr)]
chr[, (colsToConvert) := lapply(.SD, as.numeric), .SDcols = colsToConvert]

# NA's ----------------------------------------------------------------------- #
naSTR <- subset(data.frame(
  varNames = colnames(chr),
  percNA = colSums(is.na(chr))/nrow(chr)*100,
  row.names = NULL
  ), percNA == 0
)$varNames

###
naID <- data.frame(
  varNames = chr$ID,
  percNA = rowSums(is.na(chr))/ncol(chr)*100,
  row.names = NULL
)
##
setDF(chr)
chr_sub <- chr[, colnames(chr) %in% naSTR]
print(paste(ncol(chr) - length(naSTR), "STR's removed", sep = "............."))
anyNA(chr_sub) # FALSE

# COMBINE STR's and COVARIATES ----------------------------------------------- #
dat <- na.omit(merge(cov,chr_sub, by ="ID", all.y = TRUE))

# GLM MODEL PARALLEL --------------------------------------------------------- #

outcome <- "hippocampus_average"
varCov <- colnames(dat)[2:30]
varSTR <- colnames(dat)[34:ncol(dat)]

lsFit <- mclapply(varSTR,
                  mc.cores = ncores,
                  function(sSTR) {
                    sFormula <- as.formula(
                      paste(outcome,
                            paste(c(sSTR, varCov), collapse = "+"), sep = "~")
                    )
                    fit <- glm(formula = sFormula,
                                family = gaussian,
                                data = dat)
                    resCoef <- data.frame(summary(fit)$coefficients, 
                                          AIC(fit), 
                                          BIC(fit), 
                                          rsq::rsq(fit), 
                                          nrow(dat))[2,]
                    colnames(resCoef) <- c("beta", "se", "t.value", "pvalue", 
                                           "aic", "bic", "r2", "sample.size")
                    return(resCoef)
                  } )
names(lsFit) <- varSTR
resFit = do.call(rbind, lsFit)
resFit["strName"] <- rownames(resFit)
rownames(resFit) <- NULL
saveRDS(resFit, 
        paste0(getwd(), "/glm_", gsub("\\.txt", "", basename(datRepeats)), ".rds")
        )
# ---------------------------------------------------------------------------- #