library(dplyr)
library(Seurat)
library(patchwork)
library(CellChat)
options(stringsAsFactors = FALSE)

#Seurat
#Setup the Seurat object

skin.data <- Read10X(data.dir="/Users/pengyihao/Desktop/R/CCC/Data/Skin/GSE113854_RAW")
#skin.data <- ReadMtx(mtx = "count_matrix.mtx.gz", features = "genes.tsv.gz",cells = "barcodes.tsv.gz")
skin <- CreateSeuratObject(counts = skin.data, project = "skin")


#Pre-processing

#Selection and Filter

# pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")  #Selection of cells

#Visualize metrics
VlnPlot(skin, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)   
Featureplot <- FeatureScatter(skin, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

skin <- subset(skin, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)  #Filter cells

#Normalization
skin <- NormalizeData(skin, normalization.method = "LogNormalize", scale.factor = 10000)

#Identification of highly variable features
skin <- FindVariableFeatures(skin, selection.method = "vst", nfeatures = 2000)

#Scaling 
all.genes <- rownames(skin)
skin <- ScaleData(skin, features = all.genes)

#PCA
skin <- RunPCA(skin, features = VariableFeatures(object = skin))

#Visualization
VizDimLoadings(skin, dims = 1:2, reduction = "pca")
DimPlot(skin, reduction = "pca")
DimHeatmap(skin, dims = 1:2, cells = 500, balanced = TRUE)

#Determine the dimensionality(PC)
skin <- JackStraw(skin,num.replicate = 100)
skin <- ScoreJackStraw(skin, dims= 1:20)

#Visualization
JackStrawPlot(skin, dims = 1:15)
ElbowPlot(skin)

saveRDS(skin, file = "/Users/pengyihao/Desktop/R/CCC/Data/Skin/GSE113854_RAW/Output/skin.rds", version = 3)

#Cluster the Cells

skin <- FindNeighbors(skin, dims = 1:10)
skin <- FindClusters(skin, resolution = 0.5)

#UMAP
skin <- RunUMAP(skin, dims = 1:10)
#Visualization
DimPlot(skin, reduction = "umap")

#Finding DE features(cluster bioomarkers)
cluster0.markers <- FindMarkers(skin, ident.1=0, min.pct =0.25)
skin.markers <- FindAllMarkers(skin, only.pos =  TRUE, min.pct = 0.25, logfc.threshold=0.25)

#Cell Annotation
new.cluster.ids <- c("FIB-A","FIB-B","MYL-A","PC-A","FIB-C","FIB-D","FIB-E","FIB-F","PC-B","TC","PC-C","FIB-G","ENDO-A","ENDO-B","SC","DEN")
#new.cluster.ids <- c("FIB-A","FIB-B","MYL-A","PC-A","MYL-B","FIB-C","FIB-D","FIB-E","FIB-F","FIb-G","PC-B","TC","FIB-H","PC-C","ENDO-A","ENDO-B","SC","DEN","ENDO-C","TC-A","LYME","TC-B")
names(new.cluster.ids) <- levels(skin)
skin <- RenameIdents(skin, new.cluster.ids)
DimPlot(skin, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

saveRDS(skin, file = "/Users/pengyihao/Desktop/R/CCC/Data/Skin/GSE113854_RAW/Output/skinAnnotation1.rds",version = 3)



#CellChat

#Load data from Seurat 
data.input <- GetAssayData(skin, assay = "RNA", slot = "data")
labels <- Idents(skin)
meta <- data.frame(group = labels, row.names = names(labels))
colnames(meta) <- c("labels")

#Create Cellchat object
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")

#Set database
CellChatDB <- CellChatDB.mouse     #mouse or human
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB

#Preprocessing the expression data
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

#Compute communication probability
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)  #Fliter out CCC

#Extract the inferrred CCC
df.net <- subsetCommunication(cellchat)

#Infer CCC at signalling pathway level

cellchat <- computeCommunProbPathway(cellchat)

#Calculate aggregated CCC network

cellchat <- aggregateNet(cellchat)
#Visuaize whole network
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
#Visualize each cell type network
mat <- cellchat@net$weight
par(mfrow = c(2,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

#Visualiztion of signaling pathway with Hierarchy, Circle, Chord Plot
#Hierarchy Plot
pathways.show <- c("ncWNT") 
par(mfrow = c(1,1), xpd=TRUE)
vertex.receiver = seq(1,4)
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)

#Circle Plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

# Chord Plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")

#Contribution of each Ligand-Receptor piar to the sinaling pathwy
netAnalysis_contribution(cellchat, signaling = pathways.show)

#Visualize the cell-cell communication mediated by a single ligand-receptor pair


pairLR.ncWNT <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.ncWNT[1,]
# Hierarchy plot
vertex.receiver = seq(1,4)
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
# Circle plot
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
# Chord diagram
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")

#Automatically save the plots of the all inferred network for quick exploration
# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
vertex.receiver = seq(1,4)
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
}


#Visualize cell-cell communication mediated by multiple ligand-receptors or signaling pathways


# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(cellchat, sources.use = 1, targets.use = c(5:11), remove.isolate = FALSE)
# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), signaling = c("CCL","CXCL"), remove.isolate = FALSE)
# show all the significant interactions (L-R pairs) based on user's input (defined by `pairLR.use`)
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("CCL","CXCL","FGF"))
netVisual_bubble(cellchat, sources.use = c(3,4), targets.use = c(5:8), pairLR.use = pairLR.use, remove.isolate = TRUE)

#Plot the signaling gene expression distribution using violin/dot plot
plotGeneExpression(cellchat, signaling = "ncWNT")


#Systems analysis of cell-cell communication network


#Identify signaling roles (dominant senders, receivers) of cell groups as well as the major contributing signaling

#Compute network centrality score
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
#Vosualize the compute centrality score with heatmap at signalling pathway
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
#Visualize the dominant senders (sources) and receivers (targets) in a 2D space
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
gg1
#Identify signals contributing most to outgoing or incoming signaling of certain cell groups
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2

#Identify global communication patterns

library(NMF)
library(ggalluvial)

#Identify and visualize outgoing communication pattern of secreting cells
selectK(cellchat, pattern = "outgoing")
nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)

# river plot
netAnalysis_river(cellchat, pattern = "outgoing")
# dot plot
netAnalysis_dot(cellchat, pattern = "outgoing")

#Identify and visualize incoming communication pattern of target cells
selectK(cellchat, pattern = "incoming")
nPatterns = 5
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)

# river plot
netAnalysis_river(cellchat, pattern = "incoming")
# dot plot
netAnalysis_dot(cellchat, pattern = "incoming")

#Manifold and classification learning analysis of signaling networks

#Identify signaling groups based on their functional similarity
cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)

#Identify signaling groups based on structure similarity
cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "structural")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "structural", label.size = 3.5)


saveRDS(cellchat, file = "/Users/pengyihao/Desktop/R/CCC/Data/Skin/GSE113854_RAW/Output/CellChat/SkinCCC.rds")
