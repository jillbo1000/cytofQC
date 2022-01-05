## ---- include = FALSE--------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>", 
  message = FALSE, 
  warning = FALSE
)

## ----setup-------------------------------------------------------------------------------------------------------------
library(cytofQC)
library(ggplot2)

## ----------------------------------------------------------------------------------------------------------------------
f <- "../../data/FlowRepository_FR-FCM-Z29V_files/REP_1_deid.fcs"
tech <- dataPrep(f)
head(tech)

## ----------------------------------------------------------------------------------------------------------------------
labels <- labelQC(tech, types = c("bead", "doublet", "debris"))
head(labels)
table(labels$label)

## ----------------------------------------------------------------------------------------------------------------------
labels <- qcDataFrame(tech)
head(labels)
table(labels$label)

## ----eval = FALSE------------------------------------------------------------------------------------------------------
#  # K-means nearest neighbors
#  knn <- NN(tech, n = 10, type = "knn")
#  head(knn)
#  dim(knn)
#  

## ----------------------------------------------------------------------------------------------------------------------
# Mutual nearest neighbors
mnn <- NN(tech, n = 20, type = "mnn")
head(mnn)
dim(mnn)


## ----------------------------------------------------------------------------------------------------------------------
beads <- initialBead(tech, labels)
head(beads)
sum(beads$init)

## ----------------------------------------------------------------------------------------------------------------------
beadPlots <- NULL
for (i in 1:4) {
  beadPlots[[i]] <- cytofHist(x = tech[, i + 1], group = beads$init, 
                             title = paste("Bead channel", i))
}


## ----------------------------------------------------------------------------------------------------------------------
ggpubr::ggarrange(plotlist = beadPlots, ncol = 2, nrow = 2)

## ----------------------------------------------------------------------------------------------------------------------
beadPlots2 <- NULL
for (i in 1:4) {
  tmp <- data.frame(bead = tech[, i + 1], group = beads$init)
  tmp <- tmp[tmp$bead > 0, ]
  beadPlots2[[i]] <- cytofHist(x = tmp$bead, group = tmp$group, 
                               type = "density", 
                               title = paste("Bead channel", i))
}


## ----------------------------------------------------------------------------------------------------------------------
ggpubr::ggarrange(plotlist = beadPlots2, ncol = 2, nrow = 2)

## ----------------------------------------------------------------------------------------------------------------------
sure <- surematch(mnn, beads$init, threshold = 0)
head(sure)
sum(sure$Mismatch)

## ----------------------------------------------------------------------------------------------------------------------
ind <- modelData(labels, subset = sure$match, init = beads$init, n = 4000)
length(ind)
head(ind)


## ----------------------------------------------------------------------------------------------------------------------
index <- rep(FALSE, nrow(tech))
index[ind] <- TRUE
beadPlots3 <- NULL
for (i in 1:4) {
  beadPlots3[[i]] <- cytofHist(x = tech[, i + 1], group = index, 
                               type = "density", 
                               title = paste("Bead channel", i))
}


## ----------------------------------------------------------------------------------------------------------------------
ggpubr::ggarrange(plotlist = beadPlots3, ncol = 2, nrow = 2)

## ----------------------------------------------------------------------------------------------------------------------
rf.labs <- rfLabel(tech, labels, type = "bead", init = beads$init, index = ind)
head(rf.labs)
table(rf.labs$label)


## ----------------------------------------------------------------------------------------------------------------------
beadPlots <- NULL
for (i in 1:4) {
  tmp <- data.frame(bead = tech[, i + 1], group = rf.labs$label)
  tmp <- tmp[tmp$bead > 0, ]
  beadPlots[[i]] <- cytofHist(tmp$bead, group = tmp$group, type = "count", 
                              title = paste("Bead channel", i))
}


## ----------------------------------------------------------------------------------------------------------------------
ggpubr::ggarrange(plotlist = beadPlots, ncol = 2, nrow = 2)


## ----------------------------------------------------------------------------------------------------------------------
doublets <- initialDoublet(tech, rf.labs, score = 1)
head(doublets)


## ----------------------------------------------------------------------------------------------------------------------
cytofHist(doublets$doubletScore, doublets$doubletGroup, na.rm = TRUE, 
          title = "Doublet score and initial classification")

## ----------------------------------------------------------------------------------------------------------------------
sureDoublet <- surematch(mnn, doublets$init, threshold = 0)
head(sureDoublet)
sum(sureDoublet$match)

## ----------------------------------------------------------------------------------------------------------------------
mmPlots <- NULL
mmPlots[[1]] <- cytofHist(doublets$doubletScore, sureDoublet$match, type = "count", 
          title = "Matched doublets")
mmPlots[[2]] <- cytofHist(doublets$doubletScore, sureDoublet$match, type = "density", 
          title = "Matched doublets")


## ----fig.height = 3----------------------------------------------------------------------------------------------------
ggpubr::ggarrange(plotlist = mmPlots, ncol = 2, nrow = 1)

## ----------------------------------------------------------------------------------------------------------------------
d.rf.index <- modelData(rf.labs, subset = sureDoublet$match, init = doublets$init, n = 4000)


## ----------------------------------------------------------------------------------------------------------------------
index <- rep("All", nrow(tech))
index[d.rf.index] <- "Subset"
cytofHist(x = doublets$doubletScore, group = index, type = "count", 
          title = "Doublet subset")


## ----------------------------------------------------------------------------------------------------------------------
rf.labs <- rfLabel(tech, rf.labs, type = "doublet", init = doublets$init, index = d.rf.index)
head(rf.labs)
table(rf.labs$label)


## ----------------------------------------------------------------------------------------------------------------------
cytofHist(doublets$doubletScore, rf.labs$label, type = "count", 
          title = "Labels on doublet score")

## ----------------------------------------------------------------------------------------------------------------------
debris <- initialDebris(tech, rf.labs, score = 1)
head(debris)


## ----------------------------------------------------------------------------------------------------------------------
cytofHist(debris$debrisScore, debris$debrisGroup, na.rm = TRUE, 
          title = "Initial debris classification")

## ----------------------------------------------------------------------------------------------------------------------
debm <- surematch(mnn, debris$init, threshold = 0)
head(debm)
sum(debm$Mismatch)

## ----------------------------------------------------------------------------------------------------------------------
mmPlots <- NULL
mmPlots[[1]] <- cytofHist(debris$debrisScore, debm$match, type = "count", 
          title = "Mismatched debris")
mmPlots[[2]] <- cytofHist(debris$debrisScore, debm$match, type = "density", 
          title = "Mismatched debris")


## ----fig.height = 3----------------------------------------------------------------------------------------------------
ggpubr::ggarrange(plotlist = mmPlots, ncol = 2, nrow = 1)

## ----------------------------------------------------------------------------------------------------------------------
de.rf.index <- modelData(rf.labs, subset = debm$match, init = debris$init, n = 4000)


## ----------------------------------------------------------------------------------------------------------------------
index <- rep("All", nrow(tech))
index[de.rf.index] <- "Subset"
cytofHist(x = debris$debrisScore, group = index, type = "density", 
          title = "Debris subset")


## ----------------------------------------------------------------------------------------------------------------------
rf.labs <- rfLabel(tech, rf.labs, type = "debris", init = debris$init, index = de.rf.index)
head(rf.labs)
table(rf.labs$label)


## ----------------------------------------------------------------------------------------------------------------------
cytofHist(debris$debrisScore, rf.labs$label, type = "count", 
          title = "Labels on debris score")

## ----------------------------------------------------------------------------------------------------------------------
cytofHist(rowSums(tech[, c("Bead1", "Bead2", "Bead3", "Bead4")]), rf.labs$label, 
          title = "Final labels on Bead 4 measures")
tmp <- data.frame(bead = rowSums(tech[, c("Bead1", "Bead2", "Bead3", "Bead4")]), 
                  group = rf.labs$label)
tmp <- tmp[tmp$bead > 0, ]
cytofHist(tmp$bead, tmp$group, title = "Final labels on Bead measure - zeros removed")
cytofHist(doublets$doubletScore, rf.labs$label, title = "Final labels on doublet score")
cytofHist(debris$debrisScore, rf.labs$label, title = "Final labels on debris score")

## ----------------------------------------------------------------------------------------------------------------------
s.tech <- scale(tech[, -1])
lab.umap <- uwot::umap(s.tech, ret_model = TRUE)
lab.umapD <- data.frame(x = lab.umap$embedding[, 1], y = lab.umap$embedding[, 2])


## ----------------------------------------------------------------------------------------------------------------------
umapData <- data.frame(lab.umapD, labels = rf.labs$label)
umapData$bead <- ifelse(rf.labs$label == "bead", "Bead", "Not-Bead")
umapData$doublet <- ifelse(rf.labs$label == "doublet", "Doublet", "Not-Doublet")
umapData$debris <- ifelse(rf.labs$label == "debris", "Debris", "Not-Debris")
umapData$GDPzero <- ifelse(rf.labs$label == "GDPzero", "GDPzero", "Not-GDPzero")
umapData$cell <- ifelse(rf.labs$label == "cell", "Cells", "Not-Cells")

umapData$labels <- factor(umapData$labels, levels = c("cell", "GDPzero", "bead", "doublet", "debris"))


## ----------------------------------------------------------------------------------------------------------------------
require(RColorBrewer)
cols <- brewer.pal(5, "Set2")
umapG <- ggplot(umapData, aes(x = x, y = y, color = labels)) + 
  geom_point(size = 0.25) + 
  scale_color_manual(values = cols) +
  labs(title = "Random forest", x = "", y = "") +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  theme_bw()

umapG

## ----------------------------------------------------------------------------------------------------------------------
umaps <- NULL
cell_type <- levels(umapData$labels)

for (i in 1:5) {
  tmp <- data.frame(umapData[, 1:2], lab = umapData[, cell_type[i]])
  umaps[[i]] <- ggplot(tmp, aes(x = x, y = y, color = lab)) +
    geom_point(size = 0.25, show.legend = FALSE) +
    scale_color_manual("", values = c(cols[i], "gray80")) +
    labs(title = cell_type[i], x = "", y = "") +
    theme_bw()
}


## ----fig.height = 9----------------------------------------------------------------------------------------------------
library(gridExtra)
grid.arrange(umapG + guides(col = "none"), umaps[[1]], umaps[[2]], 
             umaps[[3]], umaps[[4]], umaps[[5]], nrow = 3, ncol = 2)

