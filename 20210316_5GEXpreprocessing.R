# Procedure detailed in Supplementary Figure S9A (left side)
# Installation.  ###########################################################################
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# for preprocessing
BiocManager::install("BiocNeighbors")
BiocManager::install(c('DropletUtils', 'SingleCellExperiment', 'scater', 'scran'))
# for batch correction
install.packages('devtools')
devtools::install_github('immunogenomics/harmony')
# for cell type annotation
BiocManager::install('monocle')
BiocManager::install(c('DelayedArray', 'DelayedMatrixStats', 'org.Hs.eg.db', 'org.Mm.eg.db'))
devtools::install_github('cole-trapnell-lab/garnett')


# Libraries. ###########################################################################
library(DropletUtils)
library(Matrix)
library(SingleCellExperiment)
library(scater)
library(scran)
library(BiocSingular)


# Functions. ###########################################################################
namesIDtoSymbol <- function(sceObject, x) {
  out <- rowData(sceObject)[match(x, rowData(sceObject)$ID), 'Symbol']
  x <- out
  x
}


# Paths. ###########################################################################
path.to.files <- './data/'
path.out <- './results/preprocessing/'

dir.create(path.out)


# Load data. ###########################################################################
sample.info <- read.table(paste(path.to.files, 'annotationFiles/samples.txt', sep = ''), header = T)
sce <- read10xCounts(samples = sample.info[sample.info$dataType == '5gex', 'path'],
                     sample.names = sample.info[sample.info$dataType == '5gex', 'patientID'],
                     version = c('3'))
# Add batch info and unique colnames to sce object
sce$batch <- sample.info[match(sce$Sample, sample.info$patientID), 'batch']
colnames(sce) <- paste(sce$Sample, sce$Barcode, sep = '_')


# Quality controls. ###########################################################################
# Identify mitochondrial genes and compute cell QC metrics
is.mito <- grepl("^MT-", rowData(sce)$Symbol)
sce <- addPerCellQC(sce, subsets = list(Mito = is.mito))
rm(is.mito)
# Identify high quality cells
qc.stats <- quickPerCellQC(colData(sce)[3:ncol(colData(sce))], 
                           percent_subsets = 'subsets_Mito_percent',
                           batch = sce$Sample)
sce <- sce[,!qc.stats$discard]
rm(qc.stats)

# Compute gene QC metrics
sce <- addPerFeatureQC(sce)
# Remove non-expressed genes
keep.feature <- nexprs(sce, byrow = TRUE) > 0
sce <- sce[keep.feature,]
rm(keep.feature)
# Remove uninformative genes
# Identify ribosomal genes
ribo.info <- read.table(paste(path.to.files, 
                              'annotationFiles/riboGenes_HGNC_dec2019.txt', sep = ''), header = TRUE, 
                        sep = '\t')
rowData(sce)$is.ribo <- rownames(sce) %in% ribo.info$Ensembl.gene.ID
rm(ribo.info)
# Identify bcr/tcr genes
# biomart files mainly used to remove IG|TR pseudogenes
gene.info <- read.table(paste(path.to.files, 
                              'annotationFiles/mart_export_geneLevel_GRCh38.p12_Jun2020.txt',
                              sep = ''), header = TRUE, sep = '\t')
rowData(sce)$is.var.gene <- rownames(sce) %in% unique(gene.info[grepl("^[IT]{1}[GR]{1}_[VDJ]{1}_+", 
                                                                      gene.info$Gene.type), 
                                                                'Gene.stable.ID']) | grepl("^TR[ABGD]V+|IG[HLK]V+", 
                                                                                           rowData(sce)$Symbol)
rm(gene.info) 
# keep non-ribo and non-tcr/bcr variable genes only
keep.feature <- !(rowData(sce)$is.ribo | rowData(sce)$is.var.gene)
sce <- sce[keep.feature,]
rm(keep.feature)


# Data normalization. ###########################################################################
# Compute size factors
clust.sce.block <- quickCluster(sce, method = 'igraph', use.ranks = FALSE,
                                min.mean = 0.1, BSPARAM = IrlbaParam(),
                                block = sce$Sample)
sce <- computeSumFactors(sce, cluster = clust.sce.block, min.mean = 0.1)
rm(clust.sce.block)
# log normalization (pseudo count = 1)
sce <- logNormCounts(sce)
# cpm
cpm(sce) <- calculateCPM(sce)
# tpm
tpm(sce) <- calculateTPM(sce)


# Feature selection ###########################################################################
dec.sce.sample <- modelGeneVarByPoisson(sce, block = sce$Sample)
hvgs.sample <- getTopHVGs(dec.sce.sample, prop = 0.1)


# Dimensionality reduction.  ###########################################################################
# by PCA
sce <- runPCA(sce, subset_row = hvgs.sample)
dPCs <- getDenoisedPCs(sce, technical = dec.sce.sample, subset.row = hvgs.sample)
dPCs <- ncol(dPCs$components)
reducedDim(sce, 'PCA') <- reducedDim(sce, 'PCA')[,1:dPCs]
rm(dec.sce.sample, hvgs.sample, dPCs)


# Batch correction. ###########################################################################
sce <- harmony::RunHarmony(sce, group.by.vars = c('batch', 'Sample'), 
                           reduction.save = 'harmony')


# Graphical clustering. ###########################################################################
gs <- buildSNNGraph(sce, k = 30, use.dimred = 'harmony', type = c('jaccard'))
clust.id <- igraph::cluster_louvain(gs)$membership
sce$cluster.seurat <- factor(clust.id)
rm(gs, clust.id)


# Data visualization. ###########################################################################
# by UMAP
sce <- runUMAP(sce, ncomponents = 2, dimred = 'harmony', name = 'UMAP')
plotReducedDim(sce, dimred = 'UMAP', colour_by = 'cluster.seurat', text_by = 'cluster.seurat') + theme_bw()


# Cell type annotation with garnett. ################################################################
library(garnett)
library(org.Hs.eg.db)
sce.monocle <- convertTo(sce, type = 'monocle')
pData(sce.monocle)$UMAP1 <- reducedDim(sce, 'UMAP')[,1]
pData(sce.monocle)$UMAP2 <- reducedDim(sce, 'UMAP')[,2]
# create the cluster column required by garnett
pData(sce.monocle)$garnett_cluster <- pData(sce.monocle)$cluster.seurat

# define marker genes
markerFilePath <- paste(path.to.files, 
                        'annotationFiles/markers_ENSEMBL_GRCh38.93.txt', sep = '')
markerCheck <- check_markers(sce.monocle, markerFilePath, org.Hs.eg.db, 
                             cds_gene_id_type = 'ENSEMBL', marker_file_gene_id_type = 'ENSEMBL')
markerCheck$marker_gene <- namesIDtoSymbol(sce, markerCheck$marker_gene)
plot_markers(markerCheck)

# train classifier, longest step
classifier <- train_cell_classifier(cds = sce.monocle,
                                    marker_file = markerFilePath,
                                    db = org.Hs.eg.db,
                                    cds_gene_id_type = "ENSEMBL",
                                    num_unknown = 50,
                                    marker_file_gene_id_type = "ENSEMBL")

# classify cells
sce.monocle <- classify_cells(sce.monocle, classifier = classifier, 
                              db = org.Hs.eg.db, cluster_extend = TRUE, 
                              cds_gene_id_type = 'ENSEMBL')

# Visualize annotations
# Reported in Supplementary Figure S9B
# Graphical clusters
qplot(UMAP1, UMAP2, color = garnett_cluster, data = pData(sce.monocle)) + theme_bw() + theme(legend.direction = 'horizontal', legend.position = c(0.6,0.85)) 
# Cell type
qplot(UMAP1, UMAP2, color = cell_type, data = pData(sce.monocle)) + theme_bw() + theme(legend.direction = 'horizontal', legend.position = c(0.6,0.85)) 
# Cluster extended type
qplot(UMAP1, UMAP2, color = cluster_ext_type, data = pData(sce.monocle)) + theme_bw() + theme(legend.direction = 'horizontal', legend.position = c(0.6,0.85)) 

# Add annotations to the single cell object
sce$cell_type_garnett <- pData(sce.monocle)$cell_type
sce$cluster_ext_type_garnett <- pData(sce.monocle)$cluster_ext_type


# Save preprocessed files. ################################################################
saveRDS(sce, './results/20210316_preprocessed.rds')
