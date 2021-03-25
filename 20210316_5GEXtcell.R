# Installation.  ###########################################################################
# NA


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


idPosAndNegSubsets <- function(sceObject, genes) {
  # focus on genes of interest
  filt.sce <- logcounts(sceObject[genes,])
  
  # determine marker expression per cell
  genePos <- c()
  for (c in seq(1:ncol(filt.sce))) {
    if(c %% 1000 == 0) {print(c)}
    sub <- filt.sce[,c]
    filt <- sub[sub > 0]
    if (length(filt) > 0) {
      posMarkers <- paste(namesIDtoSymbol(sceObject, names(filt)), collapse = '_')
      genePos <- c(genePos, posMarkers)
    } else {
      genePos <- c(genePos, '-')
    }
  }
  # create colnames
  geneSymbol <- sprintf('sub_%s+', paste(namesIDtoSymbol(sceObject, genes), 
                                         collapse = '+'))
  # update sceObject
  colData(sceObject)[, geneSymbol] <- genePos
  sceObject
}


# Paths. ###########################################################################
path.to.files <- './results/preprocessing/'
path.out <- './results/tcells/'

dir.create(path.out)

# Load preprocessed data. ###########################################################################
sce <- readRDS(paste(path.to.files, '20210316_preprocessed.rds', sep = ''))


# Visualize immune landscape per patient. ###########################################################
# Reported in Supplementary Figure S9C
cell.to.pt <- as.matrix(table(sce$cluster_ext_type_garnett, sce$Sample))
cell.to.pt <- as.data.frame(t((t(cell.to.pt) / apply(cell.to.pt, 2, sum)) * 100))

p <- ggplot(cell.to.pt, aes(x = Var2, y = Freq, fill = Var1))
p <- p + geom_bar(stat = 'identity')
p <- p + xlab('Patient ID') + ylab('% of cells within patient') 
p <- p + scale_fill_brewer(name = 'Cell type:', palette = 'Set3')
p <- p + theme_bw()
p
rm(p, cell.to.pt)


# Select T cells in sce object. #################################################################
# Identify T-cell clusters
cell.to.clust <- as.matrix(table(sce$cluster_ext_type_garnett, sce$cluster.seurat))
cell.indices <- apply(cell.to.clust, 2, function(x) {which(x == max(x))})
cell.to.clust <- cbind.data.frame(cluster = names(cell.indices), 
                                  cell_type_garnett = rownames(cell.to.clust)[cell.indices])
rm(cell.indices)

cluster.t <- cell.to.clust[cell.to.clust$cell_type_garnett %in% c('CD8 T cells', 
                                                                  'CD4 T cells'), 
                           'cluster.seurat']

# Filter sce object
sce <- sce[,sce$cluster.seurat %in% cluster.t]
rm(cell.to.clust, cluster.t)


# Feature selection. #################################################################
# blocking on patient
dec.sce.sample <- modelGeneVarByPoisson(sce, block = sce$Sample)
hvgs.sample <- getTopHVGs(dec.sce.sample, prop = 0.1)
rowData(sce)$is.hvgs.t <- rownames(sce) %in% hvgs.sample


# Dimensionality reduction. #################################################################
# by PCA
reducedDim(sce, 'PCA.t') <- calculatePCA(sce, subset_row = hvgs.sample)
dPCs <- getDenoisedPCs(sce, technical = dec.sce.sample, subset.row = hvgs.sample)
dPCs <- ncol(dPCs$components)
reducedDim(sce, 'PCA') <- reducedDim(sce, 'PCA.t')[,1:dPCs]
rm(dec.sce.sample, dPCs, hvgs.sample)


# Batch correction. #################################################################
# harmony
sce <- RunHarmony(sce, group.by.vars = c('batch', 'Sample'), 
                  reduction.save = 'harmony.t')


# Graphical clustering. #################################################################
gs <- buildSNNGraph(sce, k = 10, use.dimred = 'harmony.t', type = c('jaccard'))
clust.t <- igraph::cluster_louvain(gs)$membership
sce$cluster.seurat.t <- factor(clust.t)
rm(gs, clust.t)


# Visualization. #################################################################
# UMAP
sce <- runUMAP(sce, ncomponents = 2, dimred = 'harmony.t', name = 'UMAP.t')


# Doublet removal. #################################################################
for (b in unique(sce$batch)) {
  sce.sub <- sce[,sce$batch == b]
  
  # doublet by simulation
  hvgs <- rownames(rowData(sce.sub)[rowData(sce.sub)$is.hvgs.t == T,])
  d.dens <- doubletCells(sce.sub, 
                         subset.row = hvgs,
                         d = ncol(reducedDim(sce, 'harmony.t')))

  sce.sub$doublet.score<- log10(d.dens + 1)
  p <- plotColData(sce.sub, x = "cluster.seurat.t", y = "doublet.score")
  p <- p + ggtitle(sprintf('Sequencing batch %i', b)) + theme_bw()
  p
  ggsave(paste(path.out, 
               sprintf('doubletDetection_simulation_batch%i.pdf', b), 
               sep = ''), useDingbats = FALSE)
  rm(hvgs, d.dens, p)
}
rm(b)

# Manually inspect graph to select the cluster with highest doublet score
# will be removed later on
d.cells <- colnames(sce[,sce$cluster.seurat.t == 5])
sce$doublets <- colnames(sce) %in% d.cells
rm(d.cells)


# Redo data processing WITHOUT doublets.  ###########################################################
sce <- sce[,!sce$cluster.seurat.t == 5]


# Feature selection - no doublets.  #################################################################
dec.sce.sample <- modelGeneVarByPoisson(sce, block = sce$Sample)
hvgs.sample <- getTopHVGs(dec.sce.sample, prop = 0.1)
rowData(sce)$is.hvgs.t <- rownames(sce) %in% hvgs.sample

# Dimensionality reduction - no doublets.  ##########################################################
# by PCA
reducedDim(sce, 'PCA.t.nd') <- calculatePCA(sce, subset_row = hvgs.sample)
dPCs <- getDenoisedPCs(sce, technical = dec.sce.sample, subset.row = hvgs.sample)
dPCs <- ncol(dPCs$components)
reducedDim(sce, 'PCA') <- reducedDim(sce, 'PCA.t.nd')[,1:dPCs]
rm(dec.sce.sample, dPCs, hvgs.sample)


# Batch correction - no doublets.  #################################################################
sce <- RunHarmony(sce, group.by.vars = c('batch', 'Sample'), 
                  reduction.save = 'harmony.t.nd')


# Graphical clustering - no doublets.  #################################################################
gs <- buildSNNGraph(sce, k = 10, use.dimred = 'harmony.t.nd', type = c('jaccard'))
clust.seurat <- igraph::cluster_louvain(gs)$membership
sce$cluster.seurat.t.nd <- factor(clust.seurat)
rm(gs, clust.seurat)


# Visualization - no doublets.  #################################################################
# Compute UMAP projection.
sce <- runUMAP(sce, ncomponents = 2, dimred = 'harmony.t.nd', name = 'UMAP.t.nd')

# Cluster renaming.
cluster.names <- data.frame(cluster.seurat.t.nd = seq(1:10), 
                            cluster.id = c('CD8.1', 'CD4.1', 'CD8.5', 'CD8.4', 
                                           'Treg.1', 'CD4.3', 'CD8.3', 'Treg.2', 'CD4.2', 'CD8.2'))
cluster.index <- match(sce$cluster.seurat.t.nd, cluster.names$cluster.seurat.t.nd)
sce$cluster.id <- cluster.names[cluster.index, 'cluster.id']
rm(cluster.index, cluster.names)

# Generate plot reported in Figure 2A
umap.plot <- plotReducedDim(sce, dimred = 'UMAP.t.nd', colour_by = 'cluster.seurat.t.nd',
                            text_by = 'cluster.id') 
umap.plot <- umap.plot + theme_bw() + theme(legend.position = 'none') 
umap.plot <- umap.plot + xlim(-9, 9) + ylim(-5, 5)
umap.plot <- umap.plot + xlab('UMAP 1') + ylab('UMAP 2')
umap.plot
rm(umap.plot)


# Annotation of T-cell states. - no doublets.  #######################################################
library(Seurat)
# Convert sce object to seurat object
sce.symb <- sce
rownames(sce.symb) <- rowData(sce.symb)$Symbol
sce.symb$ClusterID <- sce.symb$cluster.id
sce.seurat <- as.Seurat(sce.symb[!duplicated(rownames(sce.symb)), ], counts = 'counts', data = 'logcounts')
sce.seurat <- ScaleData(sce.seurat, features = rownames(sce.seurat))
rm(sce.symb)
# Assign cluster to cells
Idents(sce.seurat) <- sce.seurat$cluster.id
# Compute markers
markers <- FindAllMarkers(sce.seurat, only.pos = TRUE, min.pct = 0.25, 
                          logfc.threshold = 0.25)
# Generate heatmap presented in Supplementary Figure S10
library(dplyr)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
DoHeatmap(subset(sce.seurat, downsample = 50), features = top10$gene, raster = FALSE) + NoLegend()
rm(markers, top10)

# Evaluate expression of cytolytic molecules in CD8+ T cells
# Generate plot presented in Supplementary Figure S11A
VlnPlot(sce.seurat, features = c('GNLY','GZMK','GZMH','GZMB', 'GZMM', 'GZMA'), 
        split.by = 'ClusterID', cols = c('lightgrey'), 
        idents = c('CD8.1', 'CD8.2', 'CD8.3', 'CD8.4', 'CD8.5'))

# Evaluate expression of genes associated with tumor-infiltrating Tregs
# Generate plot presented in Supplementary Figure S11B
VlnPlot(sce.seurat, features = c('BATF', 'CCR8', 'TNFRSF8', 'IL1R2', 'IL21R', 'CD274', 'PDCD1LG2', 
                                 'FOXP3', 'IL2RA'), split.by = 'ClusterID', cols = c('lightgrey'), 
        idents = c('Treg.1', 'Treg.2'))


# Cell phenotyping (CD39, CD103, PD1) - no doublets.  ################################################
# Triple-positive, other, triple-negative.
sce <- idPosAndNegSubsets(sce, 
                          c('ENSG00000138185', 'ENSG00000083457', 'ENSG00000188389'))
sce$tp.pheno <- ifelse(sce$`sub_ENTPD1+ITGAE+PDCD1+` == 'ENTPD1_ITGAE_PDCD1', 'tp',
                       ifelse(sce$`sub_ENTPD1+ITGAE+PDCD1+` == '-', 'tn', 'other'))
# PD1 +  or -
sce <- idPosAndNegSubsets(sce, c('ENSG00000188389'))
sce$pd1.pheno <- ifelse(sce$`sub_PDCD1+` == '-', 'PD-1-', 'PD-1+')

# Generate plot reported in Figure 2B.
umap.plot <- plotReducedDim(sce, dimred = 'UMAP.t.nd', colour_by = 'tp.pheno',
                            text_by = 'cluster.id') 
umap.plot <- umap.plot + theme_bw()
umap.plot <- umap.plot + xlim(-9, 9) + ylim(-5, 5)
umap.plot <- umap.plot + xlab('UMAP 1') + ylab('UMAP 2')
umap.plot
rm(umap.plot)


# Evaluate the enrichment of triple-positive cells per cluster - no doublets.  ########################
# Triple-positive, other, triple-negative cells per cluster - counts
tp.to.cluster <- as.matrix(table(sce$tp.pheno, sce$cluster.id))

# Triple-positive, other, triple-negative cells per cluster - percentages
tp.to.cluster.perc <- as.data.frame(t((t(tp.to.cluster) / apply(tp.to.cluster,
                                                                2, sum)) * 100))

# Generate plot reported in Figure 2D.
p <- ggplot(tp.to.cluster.perc, aes(x = Var2, y = Freq, fill = Var1))
p <- p + geom_bar(stat = 'identity')
p <- p + xlab('cluster.id') + ylab('% of cells within cluster') 
p <- p + scale_fill_grey(name = 'Cell phenotype:')
p <- p + theme_bw()
p
rm(tp.to.cluster.perc, p)


# Assess statistical significance on count data by fisher's exact test.
# Generate plot reported in Figure 2E.
tp.to.cluster <- t(tp.to.cluster)
p.value <- c()
odds.ratio <- c()
for (cluster in rownames(tp.to.cluster)) {
  a <- tp.to.cluster[cluster, 'tp']
  b <- sum(tp.to.cluster[,'tp']) - a
  c <- sum(tp.to.cluster[cluster, c('other', 'tn')])
  d <- sum(tp.to.cluster) - (a+b+c)
  res <- fisher.test(matrix(c(a, b, c, d), ncol = 2))
  p.value <- c(p.value, res$p.value)
  odds.ratio <- c(odds.ratio, res$estimate)
}
tp.to.cluster <- cbind(tp.to.cluster, p.value, odds.ratio)
rm(p.value, odds.ratio, cluster, a, b, c, d, res)

is.significant <- ifelse(tp.to.cluster[,'p.value'] < 0.05, 
                         TRUE, FALSE)
log2.odds.ratio <- ifelse(tp.to.cluster[,'odds.ratio'] > 0, 
                          log2(tp.to.cluster[,'odds.ratio']), NA)

tp.to.cluster <- cbind(tp.to.cluster, is.significant, log2.odds.ratio)
tp.to.cluster <- as.data.frame(tp.to.cluster)

p <- ggplot(tp.to.cluster, aes(x = rownames(tp.to.cluster), 
                               y = log2.odds.ratio,
                               fill = as.factor(is.significant)))
p <- p + geom_bar(stat = 'identity')
p <- p + ylim(-6, 6) + scale_fill_grey(name = 'P-value:', labels = c('>= 0.05', 
                                                                     '< 0.05'))
p <- p + theme_bw() + xlab('Cluster ID') + ylab('log2(odds ratio)')
p
rm(p, tp.to.cluster, is.significant, log2.odds.ratio)


# Triple-positive, other, triple-negative cells captured by each cluster - percentages.
# Generate plot reported in Figure 2C.
cluster.to.tp <- as.matrix(table(sce$cluster.id, sce$tp.pheno))
cluster.to.tp <- as.data.frame(t((t(cluster.to.tp) / apply(cluster.to.tp,
                                                           2, sum)) * 100))

p <- ggplot(cluster.to.tp, aes(x = Var2, y = Freq, fill = Var1))
p <- p + geom_bar(stat = 'identity')
p <- p + xlab('Cell phenotype') + ylab('% of cells') 
p <- p + scale_fill_brewer(name = 'Cluster ID:', palette = 'Set3')
p <- p + theme_bw()
p
rm(cluster.to.tp, p)


# Differential gene expression analyses. - no doublets.  #############################################
# Compare gene expression within each T-cell subset (CD8, Treg)
# Triple-positive vs. other or triple-negative
# Generate volcano plots presented in Figure 3A,B and Figure 4A,B
fc.th = 1.5
pop <- c('CD8', 'Treg')

for (p in pop) {
  keep <- grepl(p, sce$cluster.id)
  sce.sub <- sce[, keep]
  
  markers <- findMarkers(sce.sub, groups = sce.sub$tp.pheno, 
                         test.type = 't', pval.type = 'any', block = sce.sub$Sample)
  
  markers.tp <- as.data.frame(markers[['tp']])
  markers.tp <- markers.tp[markers.tp$FDR < 0.05,]
  markers.tp$Gene <- namesIDtoSymbol(sce.sub, rownames(markers.tp))
  
  for (subpop in c('logFC.tn', 'logFC.other')) {
    other <- gsub('logFC\\.', '', subpop)
    x.axis.limit <- round(max(abs(min(markers.tp[,subpop])), 
                              abs(max(markers.tp[,subpop])))) + 1
    
    volcano <- ggplot(markers.tp, aes(x = markers.tp[,subpop], y = -log10(FDR)))
    volcano <- volcano + geom_point(shape = 16, size = 2, 
                                    colour = ifelse(abs(markers.tp[,subpop]) > 
                                                      log2(fc.th),
                                                    'steelblue','darkgrey'), 
                                    alpha = 0.8)
    volcano <- volcano + geom_text(aes(label = ifelse(abs(markers.tp[,subpop]) >
                                                        log2(fc.th), 
                                                      Gene, '')),
                                   hjust = 0, vjust = 0)
    volcano <- volcano + xlim(-x.axis.limit, x.axis.limit) + theme_bw()
    # volcano <- volcano + xlim(0.5, 1) + ylim(15, 35) # used to generate zoom
    volcano <- volcano + xlab(sprintf('log2(%s/%s)', 'tp', other))
    print(volcano)
    ggsave(paste(path.out, sprintf('%s_%svolcanoPlot_tpvs%s.pdf', 
                                   today, p, other), 
                 sep = ''), 
           useDingbats = FALSE)
    
    to.save <- markers.tp[,c('FDR', subpop, 'Gene')]
    to.save$selected <- ifelse(abs(to.save[,subpop]) > log2(fc.th), TRUE, FALSE)
    write.table(to.save, paste(path.out, sprintf('%s_%s_degsList_tpvs%s.txt',
                                                 today, p, other), sep = ''),
                quote = FALSE)
  }
}
rm(markers, markers.tp, sce.sub, to.save, volcano, fc.th, keep, other, p, 
   subpop, x.axis.limit)

# Compare gene expression within each T-cell subset (CD8, Treg)
# Triple-positive vs. all other subsets (3 single positives and 3 double positives)
# Generate volcano plots presented in Supplementary Figure S13A and S15A
fc.th <- 1.5
# S13A
p = 'CD8'
subs <- c('CD8.1', 'CD8.2', 'CD8.3', 'CD8.4', 'CD8.5')
## S15A
# p = 'Treg'
# subs <- c('Treg.1', 'Treg.2')
sce.sub.degs <- sce[,sce$cluster.id %in% subs]

markers.degs <- findMarkers(sce.sub.degs, groups = sce.sub.degs$`sub_ENTPD1+ITGAE+PDCD1+`, 
                            test.type = 't', pval.type = 'any', block = sce.sub.degs$Sample)
markers.degs.tp <- as.data.frame(markers.degs[['ENTPD1_ITGAE_PDCD1']])
markers.degs.tp <- markers.degs.tp[markers.degs.tp$FDR < 0.05,]
markers.degs.tp$Gene <- namesIDtoSymbol(sce.sub.degs, rownames(markers.degs.tp))

subpops <- colnames(markers.degs.tp)[6:11]
for (subpop in subpops) {
  other <- gsub('logFC\\.', '', subpop)
  x.axis.limit <- round(max(abs(min(markers.degs.tp[,subpop])), 
                            abs(max(markers.degs.tp[,subpop]))) +1)
  
  volcano <- ggplot(markers.degs.tp, aes(x = markers.degs.tp[,subpop], y = -log10(FDR)))
  volcano <- volcano + geom_point(shape = 16, size = 2,
                                  colour = ifelse(abs(markers.degs.tp[,subpop]) >
                                                    log2(fc.th),
                                                  'steelblue','darkgrey'),
                                  alpha = 0.8)
  volcano <- volcano + geom_text(aes(label = ifelse(abs(markers.degs.tp[,subpop]) >
                                                      log2(fc.th),
                                                    Gene, '')),
                                 hjust = 0, vjust = 0)
  volcano <- volcano + xlim(-x.axis.limit, x.axis.limit) + theme_bw()
  # volcano <- volcano + xlim(0.5, 1) + ylim(15, 35) # used to generate zoom
  volcano <- volcano + xlab(sprintf('log2(%s/%s)', 'tp', other))
  print(volcano)
  ggsave(paste(path.out, sprintf('/%svolcanoPlot_tpvs%s.pdf', p, other), 
               sep = ''), 
         useDingbats = FALSE)
}

# Compare gene expression within each T-cell subset (CD103-CD39-, CD103-CD39+, CD103+CD39-, CD103+CD39+)
# CD8 T cells vs. CD4 T cells (CD4 and Tregs)
# Generate volcano plots presented in Supplementary Figure S6B
fc.th <- 1.5
# p <- 'CD103n_CD39n'
# subs <- c('-', 'PDCD1')
# p <- 'CD103n_CD39p'
# subs <- c('ENTPD1_PDCD1', 'ENTPD1')
# p <- 'CD103p_CD39n'
# subs <- c('ITGAE_PDCD1', 'ITGAE')
p <- 'CD103p_CD39p'
subs <- c('ENTPD1_ITGAE', 'ENTPD1_ITGAE_PDCD1')
sce.sub.degs <- sce[,sce$`sub_ENTPD1+ITGAE+PDCD1+` %in% subs]

markers.degs <- findMarkers(sce.sub.degs, groups = sce.sub.degs$cluster_ext_type_garnett, 
                            test.type = 't', pval.type = 'any', block = sce.sub.degs$Sample)
markers.degs.tp <- as.data.frame(markers.degs[['CD8 T cells']])
markers.degs.tp <- markers.degs.tp[markers.degs.tp$FDR < 0.05,]
markers.degs.tp$Gene <- namesIDtoSymbol(sce.sub.degs, rownames(markers.degs.tp))

subpops <- colnames(markers.degs.tp)[5]
for (subpop in subpops) {
  other <- gsub('logFC\\.', '', subpop)
  x.axis.limit <- round(max(abs(min(markers.degs.tp[,subpop])), 
                            abs(max(markers.degs.tp[,subpop]))) +1)
  
  volcano <- ggplot(markers.degs.tp, aes(x = markers.degs.tp[,subpop], y = -log10(FDR)))
  volcano <- volcano + geom_point(shape = 16, size = 2,
                                  colour = ifelse(abs(markers.degs.tp[,subpop]) >
                                                    log2(fc.th),
                                                  'steelblue','darkgrey'),
                                  alpha = 0.8)
  volcano <- volcano + geom_text(aes(label = ifelse(abs(markers.degs.tp[,subpop]) >
                                                      log2(fc.th),
                                                    Gene, '')),
                                 hjust = 0, vjust = 0)
  volcano <- volcano + xlim(-x.axis.limit, x.axis.limit) + theme_bw()
  # volcano <- volcano + xlim(0.5, 1) + ylim(15, 35) # used to generate zoom
  volcano <- volcano + xlab(sprintf('log2(%s/%s)', 'CD8.T.cells', other))
  print(volcano)
  ggsave(paste(path.out, sprintf('/%svolcanoPlot_CD8.T.cellsvs%s.pdf', p, other), 
               sep = ''), 
         useDingbats = FALSE)
}

# Evaluate Trm gene signature in CD8 T cells
# Genes from: PMID: 30555481
# Triple-positive vs. other vs. triple-negative
# Generate violin plots presented in Supplementary Figure S12A and B
high <- c('ITGAE', 'CD69', 'ITGA1', 'CXCR6', 'PDCD1', 'HAVCR2', 'LAG3', 'TIGIT', 'CD44')
low <- c('KLRG1', 'CCR7', 'SELL','S1PR1')

sce.cd8 <- sce[,sce$cluster.id %in% c('CD8.1', 'CD8.2', 'CD8.3', 'CD8.4', 'CD8.5')]
rownames(sce.cd8) <- rowData(sce.cd8)$Symbol
plotExpression(sce.cd8, high, x = 'tp.pheno') + theme_bw()
ggsave(paste(path.out, sprintf('/%s_expressionTRMhigh.pdf', 'CD8'), 
             sep = ''), useDingbats = FALSE)
plotExpression(sce.cd8, low, x = 'tp.pheno') + theme_bw()
ggsave(paste(path.out, sprintf('/%s_expressionTRMlow.pdf', 'CD8'), 
             sep = ''), useDingbats = FALSE)
rm(sce.cd8)


# Save post-analysis T-cell file. ################################################################
saveRDS(sce, './results/tcells/20210118_postAnalysis_alltcells.rds')

