library(SingleCellExperiment)
library(Seurat)
library(monocle3)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(grid)
library(CellChat)
library(OmnipathR)
library(future)


# population.size: This argument specifies whether to adjust the communication probability based on the size of the cell populations involved in the communication. When set to TRUE, the function will adjust the calculated communication probabilities by considering the relative sizes of the sender and receiver cell populations. This can be important in cases where there are significant differences in the sizes of cell populations, as larger populations might inherently have a higher likelihood of communication simply due to having more cells.
# Setting population.size = FALSE: When you set this parameter to FALSE, it means that the computation of communication probabilities will not take into account the size of the cell populations. This is useful in scenarios where you want to focus purely on the potential for communication based on the expression of signaling molecules, receptors, and ligands, without adjusting for the number of cells in each population. This can provide insights into the intrinsic communication potential that is not influenced by population size.

#### TUMOR ########
file_path <- "/Users/juliat/Desktop/bachelor/daten 7.09.23/RDEB data/cds.RDS"
cds.t <- readRDS(file_path)

# Subset the SingleCellExperiment object based on the 'tissue' column in colData
cds_tumor <- cds.t[, cds.t$tissue == "tumor"]

logcnts.t <- normalizeData(assay(cds_tumor, "counts"), scale.factor = 1e6, do.log = TRUE)

cellchat.t <- createCellChat(object = logcnts.t, meta = data.frame(colData(cds_tumor)), group.by = "DE") 

# Setting the ligand-receptor interaction database 
CellChatDB <- CellChatDB.human
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
cellchat.t@DB <- CellChatDB.use
cellchat.t <- subsetData(cellchat.t) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # Use multisession backend
cellchat.t <- identifyOverExpressedGenes(cellchat.t)
cellchat.t <- identifyOverExpressedInteractions(cellchat.t)

levels(cellchat.t@idents) <- c ("B cell", "Endothel", "Endothel_vasc._31", "Erythrocyte", "Fibroblast_18", "Fibroblast_20", "Fibroblast_25", "Fibroblast_32", "Fibroblast_36", "T gd", "Keratinocyte_26", "Mast_27", "Monocyte", "Monocyte_23", "Neurons_47", "Neutrophil", "NK", "NK T like", "Plasma_15", "Plasma_35", "Plasma_37", "Platelet", "Schwann_43", "T cytotox. CD8", "T fh", "T mem 1", "T mem 2", "T naive CD4", "T reg CD4", "T undef.")
cellchat.t@idents<- droplevels(cellchat.t@idents)
cellchat.t <- computeCommunProb(cellchat.t, type = "triMean", population.size = FALSE)
cellchat.t <- filterCommunication(cellchat.t) # filter communications occurring with less than 10 cells in a cell group
# The cell-cell communication related with the following cell groups are excluded due to the few number of cells:  T naive CD4 

cellchat.t <- aggregateNet(cellchat.t)
df.net.t <- subsetCommunication(cellchat.t) #  returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors
df.net.t

# Compute the communication probability on signaling pathway level by summarizing all related ligands/receptors
cellchat.t <- computeCommunProbPathway(cellchat.t)


# Save data ---------------------------------------------------------------

cellchat.t %>% saveRDS("/Users/juliat/Desktop/bachelor/Final thesis/ccc_cellchat.t.2.rds")
df.net.t %>% saveRDS("/Users/juliat/Desktop/bachelor/Final thesis/ccc_df.net.t.2.rds")

cellchat.t
cellchat.t@netP$pathways

#### SKIN ########
file_path <- "/Users/juliat/Desktop/bachelor/daten 7.09.23/RDEB data/skin_only_cell_data_set.RDS"
cds.s <- readRDS(file_path)

logcnts.s <- normalizeData(assay(cds.s, "counts"), scale.factor = 1e6, do.log = TRUE)

cellchat.s <- createCellChat(object = logcnts.s, meta = data.frame(colData(cds.s)), group.by = "DE") 

# Setting the ligand-receptor interaction database 
CellChatDB <- CellChatDB.human
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
cellchat.s@DB <- CellChatDB.use
cellchat.s <- subsetData(cellchat.s) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # Use multisession backend
# future::plan("multiprocess", workers = 4) # do parallel # No such strategy for futures: ‘multiprocess’
cellchat.s <- identifyOverExpressedGenes(cellchat.s)
cellchat.s <- identifyOverExpressedInteractions(cellchat.s)

levels(cellchat.s@idents) <- c ("B cell", "Endothel", "Endothel_vasc._31", "Erythrocyte", "Fibroblast_18", "Fibroblast_20", "Fibroblast_25", "Fibroblast_32", "Fibroblast_36", "T gd", "Keratinocyte_26", "Mast_27", "Monocyte", "Monocyte_23", "Neurons_47", "Neutrophil", "NK", "NK T like", "Plasma_15", "Plasma_35", "Plasma_37", "Platelet", "Schwann_43", "T cytotox. CD8", "T fh", "T mem 1", "T mem 2", "T naive CD4", "T reg CD4", "T undef.")
cellchat.s@idents<- droplevels(cellchat.s@idents)
cellchat.s <- computeCommunProb(cellchat.s, type = "triMean", population.size = FALSE) 
cellchat.s <- filterCommunication(cellchat.s, min.cells = 10) #? # filter communications occurring with less than 10 cells in a cell group
# The cell-cell communication related with the following cell groups are excluded due to the few number of cells:  T naive CD4 

cellchat.s <- aggregateNet(cellchat.s)
df.net.s <- subsetCommunication(cellchat.s) #  returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors
df.net.s

# Compute the communication probability on signaling pathway level by summarizing all related ligands/receptors
cellchat.s <- computeCommunProbPathway(cellchat.s)

# Save data ---------------------------------------------------------------
cellchat.s %>% saveRDS("/Users/juliat/Desktop/bachelor/Final thesis/ccc_cellchat.s.2.rds")
df.net.s %>% saveRDS("/Users/juliat/Desktop/bachelor/Final thesis/ccc_df.net.s.2.rds")