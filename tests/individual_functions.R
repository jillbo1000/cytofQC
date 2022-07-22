library(cytofQC)
library(CATALYST)
library(SingleCellExperiment)

# f_data_all <- "R:/jill/cytof/data/sarah/data_for_prelim_clustering/20210611_U2OS_NT-1_GFP_spike-in_Cellsave_vs_EDTA/0hr_cellsave_no-perm_01_0.fcs"
# f_data_all <- "R:/jill/cytof/data/sarah/data_for_prelim_clustering/20210611_U2OS_NT-1_GFP_spike-in_Cellsave_vs_EDTA/0hr_EDTA_no-perm_01_1.fcs"
# f_data_all <- "R:/jill/cytof/data/sarah/data_for_prelim_clustering/20210611_U2OS_NT-1_GFP_spike-in_Cellsave_vs_EDTA/24hr_cellsave_no-perm_01_1.fcs"
# f_data_all <- "R:/jill/cytof/data/sarah/data_for_prelim_clustering/20210826_U2OS_NT-1_GFP_spike-in_for_Jill_and_Kelly/U2OS_NT-1_GFP_Spike_NORMAL_01_1.fcs"
# f_data_all <- "R:/jill/cytof/data/sarah/data_for_prelim_clustering/20220202_header_fixed/A673_Processed.fcs"
# f_data_all <- "R:/jill/cytof/data/FlowRepository_FR-FCM-Z29V_files/REP_1_deid.fcs"
f_data_all <- "R:/jill/cytof/data/FlowRepository_FR-FCM-Z29V_files/REP_25_deid.fcs"


# x <- CATALYST::prepData(f_data_all)
# names(int_colData(x))
# rownames(x)

## SVM

# sce <- readCytof(f_data_all, 
#                  beads = c("Beads"),
#                  dna = c("DNA1", "DNA2"),
#                  event_length = "Event_length",
#                  viability = "Viability",
#                  gaussian = c("Center", "Offset", "Width", "Residual"))

sce <- readCytof(f_data_all)

sce <- initialBead(sce)
cytofHist(sce$scores$beadScore, sce$initial$beadInitial)
cytofHist(sce$scores$beadScore, sce$initial$beadInitial, type = "density")
sce <- svmLabel(sce, type = "bead")
cytofHist(sce$scores$beadScore, sce$label)
cytofHist(sce$scores$beadScore, sce$label, type = "density")
cytofHist(sce$scores$beadScore, sce$label, type = "density") +
    ggplot2::lims(x = c(6, 30)) 


sce <- initialDebris(sce)
cytofHist(sce$scores$debrisScore, sce$initial$debrisInitial)
cytofHist(sce$scores$debrisScore, sce$initial$debrisInitial, type = "density")
cytofHist(sce$scores$debrisScore, sce$label)
cytofHist(sce$scores$debrisScore, sce$label, type = "density")
sce <- svmLabel(sce, type = "debris")
cytofHist(sce$scores$debrisScore, sce$label)
cytofHist(sce$scores$debrisScore, sce$label, type = "density")

sce <- initialDoublet(sce)
cytofHist(sce$scores$doubletScore, sce$initial$doubletInitial)
cytofHist(sce$scores$doubletScore, sce$initial$doubletInitial, type = "density")
cytofHist(sce$scores$doubletScore, sce$label)
cytofHist(sce$scores$doubletScore, sce$label, type = "density")
sce <- svmLabel(sce, type = "doublet")
cytofHist(sce$scores$doubletScore, sce$label)
cytofHist(sce$scores$doubletScore, sce$label, type = "density")

sce <- initialDead(sce)
cytofHist(sce$scores$deadScore, sce$initial$deadInitial)
cytofHist(sce$scores$deadScore, sce$initial$deadInitial, type = "density")
cytofHist(sce$scores$deadScore, sce$label)
cytofHist(sce$scores$deadScore, sce$label, type = "density")
sce <- svmLabel(sce, type = "dead")
cytofHist(sce$scores$deadScore, sce$label)
cytofHist(sce$scores$deadScore, sce$label, type = "density")


## GBM

# sce <- readCytof(f_data_all, 
#                  beads = c("Beads"),
#                  dna = c("DNA1", "DNA2"),
#                  event_length = "Event_length",
#                  viability = "Viability",
#                  gaussian = c("Center", "Offset", "Width", "Residual"))

sce <- readCytof(f_data_all)

sce <- initialBead(sce)
cytofHist(sce$scores$beadScore, sce$initial$beadInitial)
cytofHist(sce$scores$beadScore, sce$initial$beadInitial, type = "density")
sce <- gbmLabel(sce, type = "bead")
cytofHist(sce$scores$beadScore, sce$label)
cytofHist(sce$scores$beadScore, sce$label, type = "density") +
    ggplot2::lims(x = c(20, 30))

sce <- initialDebris(sce)
cytofHist(sce$scores$debrisScore, sce$initial$debrisInitial)
cytofHist(sce$scores$debrisScore, sce$initial$debrisInitial, type = "density")
cytofHist(sce$scores$debrisScore, sce$label)
cytofHist(sce$scores$debrisScore, sce$label, type = "density")
sce <- gbmLabel(sce, type = "debris")
cytofHist(sce$scores$debrisScore, sce$label)
cytofHist(sce$scores$debrisScore, sce$label, type = "density")

sce <- initialDoublet(sce)
cytofHist(sce$scores$doubletScore, sce$initial$doubletInitial)
cytofHist(sce$scores$doubletScore, sce$initial$doubletInitial, type = "density")
cytofHist(sce$scores$doubletScore, sce$label)
cytofHist(sce$scores$doubletScore, sce$label, type = "density")
sce <- gbmLabel(sce, type = "doublet")
cytofHist(sce$scores$doubletScore, sce$label)
cytofHist(sce$scores$doubletScore, sce$label, type = "density")

sce <- initialDead(sce)
cytofHist(sce$scores$deadScore, sce$initial$deadInitial)
cytofHist(sce$scores$deadScore, sce$initial$deadInitial, type = "density")
cytofHist(sce$scores$deadScore, sce$label)
cytofHist(sce$scores$deadScore, sce$label, type = "density")
sce <- gbmLabel(sce, type = "dead")
cytofHist(sce$scores$deadScore, sce$label)
cytofHist(sce$scores$deadScore, sce$label, type = "density")


## RF

sce <- readCytof(f_data_all, 
                 beads = c("Beads"),
                 dna = c("DNA1", "DNA2"),
                 event_length = "Event_length",
                 viability = "Viability",
                 gaussian = c("Center", "Offset", "Width", "Residual"))

sce <- readCytof(f_data_all) 

sce <- initialBead(sce)
cytofHist(sce$scores$beadScore, sce$initial$beadInitial)
cytofHist(sce$scores$beadScore, sce$initial$beadInitial, type = "density")
sce <- rfLabel(sce, type = "bead")
cytofHist(sce$scores$beadScore, sce$label)
cytofHist(sce$scores$beadScore, sce$label, type = "density") +
    ggplot2::lims(x = c(20, 30))

sce <- initialDebris(sce)
cytofHist(sce$scores$debrisScore, sce$initial$debrisInitial)
cytofHist(sce$scores$debrisScore, sce$initial$debrisInitial, type = "density")
cytofHist(sce$scores$debrisScore, sce$label)
cytofHist(sce$scores$debrisScore, sce$label, type = "density")
sce <- rfLabel(sce, type = "debris")
cytofHist(sce$scores$debrisScore, sce$label)
cytofHist(sce$scores$debrisScore, sce$label, type = "density")

sce <- initialDoublet(sce)
cytofHist(sce$scores$doubletScore, sce$initial$doubletInitial)
cytofHist(sce$scores$doubletScore, sce$initial$doubletInitial, type = "density")
cytofHist(sce$scores$doubletScore, sce$label)
cytofHist(sce$scores$doubletScore, sce$label, type = "density")
sce <- rfLabel(sce, type = "doublet")
cytofHist(sce$scores$doubletScore, sce$label)
cytofHist(sce$scores$doubletScore, sce$label, type = "density")

sce <- initialDead(sce)
cytofHist(sce$scores$deadScore, sce$initial$deadInitial)
cytofHist(sce$scores$deadScore, sce$initial$deadInitial, type = "density")
cytofHist(sce$scores$deadScore, sce$label)
cytofHist(sce$scores$deadScore, sce$label, type = "density")
sce <- rfLabel(sce, type = "dead")
cytofHist(sce$scores$deadScore, sce$label)
cytofHist(sce$scores$deadScore, sce$label, type = "density")


## Semi-supervised SVM

# sce <- readCytof(f_data_all, 
#                  beads = c("Beads"),
#                  dna = c("DNA1", "DNA2"),
#                  event_length = "Event_length",
#                  viability = "Viability",
#                  gaussian = c("Center", "Offset", "Width", "Residual"))

sce <- readCytof(f_data_all)

sce <- initialBead(sce)
cytofHist(sce$scores$beadScore, sce$initial$beadInitial)
cytofHist(sce$scores$beadScore, sce$initial$beadInitial, type = "density")
sce <- s3vmLabel(sce, type = "bead")
cytofHist(sce$scores$beadScore, sce$label)
cytofHist(sce$scores$beadScore, sce$label, type = "density")
cytofHist(sce$scores$beadScore, sce$label, type = "density") +
    ggplot2::lims(x = c(20, 30)) 


sce <- initialDebris(sce)
cytofHist(sce$scores$debrisScore, sce$initial$debrisInitial)
cytofHist(sce$scores$debrisScore, sce$initial$debrisInitial, type = "density")
cytofHist(sce$scores$debrisScore, sce$label)
cytofHist(sce$scores$debrisScore, sce$label, type = "density")
sce <- s3vmLabel(sce, type = "debris")
cytofHist(sce$scores$debrisScore, sce$label)
cytofHist(sce$scores$debrisScore, sce$label, type = "density")

sce <- initialDoublet(sce)
cytofHist(sce$scores$doubletScore, sce$initial$doubletInitial)
cytofHist(sce$scores$doubletScore, sce$initial$doubletInitial, type = "density")
cytofHist(sce$scores$doubletScore, sce$label)
cytofHist(sce$scores$doubletScore, sce$label, type = "density")
sce <- s3vmLabel(sce, type = "doublet")
cytofHist(sce$scores$doubletScore, sce$label)
cytofHist(sce$scores$doubletScore, sce$label, type = "density")

sce <- initialDead(sce)
cytofHist(sce$scores$deadScore, sce$initial$deadInitial)
cytofHist(sce$scores$deadScore, sce$initial$deadInitial, type = "density")
cytofHist(sce$scores$deadScore, sce$label)
cytofHist(sce$scores$deadScore, sce$label, type = "density")
sce <- s3vmLabel(sce, type = "dead")
cytofHist(sce$scores$deadScore, sce$label)
cytofHist(sce$scores$deadScore, sce$label, type = "density")


# labelQC

# sce <- readCytof(f_data_all, 
#                  beads = c("Beads"),
#                  dna = c("DNA1", "DNA2"),
#                  event_length = "Event_length",
#                  viability = "Viability",
#                  gaussian = c("Center", "Offset", "Width", "Residual"))

sce <- readCytof(f_data_all)

sceS <- labelQC(sce, model = "svm")
cytofHist(sceS$scores$beadScore, sceS$label)
cytofHist(sceS$scores$beadScore, sceS$label, type = "density") 
cytofHist(sceS$scores$beadScore, sceS$label, type = "density") +
    ggplot2::lims(x = c(20, 30))
cytofHist(sceS$scores$debrisScore, sceS$label)
cytofHist(sceS$scores$debrisScore, sceS$label, type = "density")
cytofHist(sceS$scores$doubletScore, sceS$label)
cytofHist(sceS$scores$doubletScore, sceS$label, type = "density")
cytofHist(sceS$scores$deadScore, sceS$label)
cytofHist(sceS$scores$deadScore, sceS$label, type = "density")

sceG <- labelQC(sce, model = "gbm")
cytofHist(sceG$scores$beadScore, sceG$label)
cytofHist(sceG$scores$beadScore, sceG$label, type = "density")
cytofHist(sceS$scores$beadScore, sceS$label, type = "density") +
    ggplot2::lims(x = c(6, 30))
cytofHist(sceG$scores$debrisScore, sceG$label)
cytofHist(sceG$scores$debrisScore, sceG$label, type = "density")
cytofHist(sceG$scores$doubletScore, sceG$label)
cytofHist(sceG$scores$doubletScore, sceG$label, type = "density")
cytofHist(sceG$scores$deadScore, sceG$label)
cytofHist(sceG$scores$deadScore, sceG$label, type = "density")
cytofHist(sceG$scores$deadScore, sceG$label) +
    ggplot2::lims(x = c(10, 40))

sceR <- labelQC(sce, model = "rf")
cytofHist(sceR$scores$beadScore, sceR$label)
cytofHist(sceR$scores$beadScore, sceR$label, type = "density")
cytofHist(sceS$scores$beadScore, sceS$label, type = "density") +
    ggplot2::lims(x = c(6, 30))
cytofHist(sceR$scores$debrisScore, sceR$label)
cytofHist(sceR$scores$debrisScore, sceR$label, type = "density")
cytofHist(sceR$scores$doubletScore, sceR$label)
cytofHist(sceR$scores$doubletScore, sceR$label, type = "density")
cytofHist(sceR$scores$deadScore, sceR$label)
cytofHist(sceR$scores$deadScore, sceR$label, type = "density")
cytofHist(sceR$scores$deadScore, sceR$label) +
    ggplot2::lims(x = c(10, 40))


