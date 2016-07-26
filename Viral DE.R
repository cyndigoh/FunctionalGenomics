# Viral Differential Expression

# Load packages
library(limma)

# Load volcano plot function
source("../VolcanoPlotFunction.R")

# read in data 
data.exprs <- read.delim("../../Data/GAinS_normalised_265samples.txt", 
                         row.names=1)
# removes first column of probe labels/names to make the analysis easier
data.exprs <- data.exprs[, 2:266]
colnames(data.exprs) <- gsub("X", "", colnames(data.exprs))
# read in the covariate data so you know which individuals are viral/non-viral
s.info <- read.delim("../../Data/covariates_265samples_noNA.txt")

# Define variables to include in model
viral <- as.factor(s.info$Viral)

# Design and fit model
design <- model.matrix( ~ 0 + viral)
fit <- lmFit(data.exprs, design)
contrast.matrix <- makeContrasts(Viral=viral1 - viral0, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Results
# Sort by p-value
allhitsv <- topTable(fit2, coef="Viral", adjust.method="fdr", sort.by="P", 
                     num=nrow(data.exprs))
# relabel with probeID's
allhitsv <- data.frame(cbind("ProbeID"=rownames(allhitsv), allhitsv))
# Table of hits which are significant at p<0.05
hitsv <- topTable(fit2, coef="Viral", adjust.method="fdr", sort.by="P", 
                  p.value=0.05, num=nrow(data.exprs))
hitsv <- data.frame(cbind("ProbeID"=rownames(hitsv), hitsv))
# Table of probes differentially expressed at >1.5 fold change and p<0.05
hits1.5v <- topTable(fit2, coef="Viral", lfc=log2(1.5), p.value=0.05, 
                     sort.by="P", adjust.method="fdr", num=nrow(data.exprs))
hits1.5v <- data.frame(cbind("ProbeID"=row.names(hits1.5v), hits1.5v))

# get genelists from probe IDs
ann.data <- read.delim("../../Data/HumanHT-12_V4_0_R2_15002873_B.txt")
allhitsv$Gene <- ann.data$ILMN_Gene[match(allhitsv$ProbeID, 
                                          ann.data$Array_Address_Id)]
hitsv$Gene <- ann.data$ILMN_Gene[match(hitsv$ProbeID, 
                                          ann.data$Array_Address_Id)]
hits1.5v$Gene <- ann.data$ILMN_Gene[match(hits1.5v$ProbeID, 
                                          ann.data$Array_Address_Id)]

# volcano plot
pdf("N:/jknight/Data/GAinS/Expression/CAP vs FP paper/Figures/Fig 2G Viral VP.pdf", useDingbats = FALSE)
volcano.plot(allhitsv$logFC, allhitsv$adj.P.Val, allhitsv$Gene, n.label=10)
dev.off()

# Restrict to flu versus non-viral

# remove non-flu viral
remove <- which(s.info$Viral == "1" & s.info$H1N1 == "0")
data.exprs <- data.exprs[, -(remove)]
s.info <- s.info[-(remove), ]

# Define variables to include in model
viral <- as.factor(s.info$Viral)

# Design and fit model
design <- model.matrix( ~ 0 + viral)
fit <- lmFit(data.exprs, design)
contrast.matrix <- makeContrasts(Viral=viral1 - viral0, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Results
allhitsv <- topTable(fit2, coef="Viral", adjust.method="fdr", sort.by="P", 
                     num=nrow(data.exprs))
allhitsv <- data.frame(cbind("ProbeID"=rownames(allhitsv), allhitsv))
hitsv <- topTable(fit2, coef="Viral", adjust.method="fdr", sort.by="P", 
                  p.value=0.05, num=nrow(data.exprs))
hitsv <- data.frame(cbind("ProbeID"=rownames(hitsv), hitsv))
hits1.5v <- topTable(fit2, coef="Viral", lfc=log2(1.5), p.value=0.05, 
                     sort.by="P", adjust.method="fdr", num=nrow(data.exprs))
hits1.5v <- data.frame(cbind("ProbeID"=row.names(hits1.5v), hits1.5v))

# get genelists from probe IDs
allhitsv$Gene <- ann.data$ILMN_Gene[match(allhitsv$ProbeID, 
                                          ann.data$Array_Address_Id)]
hitsv$Gene <- ann.data$ILMN_Gene[match(hitsv$ProbeID, 
                                       ann.data$Array_Address_Id)]
hits1.5v$Gene <- ann.data$ILMN_Gene[match(hits1.5v$ProbeID, 
                                          ann.data$Array_Address_Id)]

# volcano plot
pdf("N:/jknight/Data/GAinS/Expression/CAP vs FP paper/Figures/Fig 2G Flu VP.pdf", useDingbats = FALSE)
volcano.plot(allhitsv$logFC, allhitsv$adj.P.Val, allhitsv$Gene, n.label=10)
dev.off()
