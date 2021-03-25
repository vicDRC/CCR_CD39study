# Installation.  ###########################################################################
BiocManager::install("scRepertoire")


# Libraries. ###########################################################################
library(scRepertoire)
library(SingleCellExperiment)
library(scater)


# Functions. ###########################################################################
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
path.to.tcr <- './data/'
path.to.gex <- './results/tcells/'
path.out <- './results/tcr/'

dir.create(path.out)


# Load data. ###########################################################################
# TCR
pt.list <- list('cd39-01-s1', 'cd39-01-s2', 'cd39-01-s3', 'cd39-01-s4', 'cd39-01-s5', 'cd39-01-s6')
sample.type <- rep('HGSC', 6)

df.list <- lapply(pt.list, function(x) {
  path.to.file <- paste(path.to.tcr, x, '/tcr/filtered_contig_annotations.csv', sep = '')
  tcr <- read.table(path.to.file, sep = ',', header = TRUE)
  return(tcr)
})
names(df.list) <- pt.list

df.combined <- combineTCR(df.list, samples = unlist(pt.list), ID = sample.type, 
                          cells = 'T-AB', 
                          removeNA = TRUE, filterMulti = TRUE)

# 5GEX
sce <- readRDS(paste(path.to.gex, '20210316_tcells.rds', sep = ''))

# Merge TCR and 5GEX
# Make barcode names identical to the one in the TCR object
colnames(sce) <- gsub('_', '_HGSC_', colnames(sce))
# Actual merge
sce.tcr <- combineExpression(df.combined, sce, cloneCall = 'gene+nt', groupBy = 'sample')
# Remove cells without TCR info
sce.tcr <- sce.tcr[,!is.na(sce.tcr$barcode)]


# Clone frequency analysis. ###########################################################################
# Generate dot plots presented in Figure 3C and 4C, as well as Supplementary Figure S13B and S15B
sce.tcr.df <- as.data.frame(colData(sce.tcr))
# Make sure each clone is here one time only (use CT strict because cloneCall='gene+nt')
sce.tcr.df <- sce.tcr.df[!duplicated(sce.tcr.df$CTstrict),]

for (p in c('CD8', 'Treg')) {
  # Keep T cells of interest
  sce.tcr.filt <- sce.tcr.df[grepl(p, sce.tcr.df$cluster.id),]
  
  # tp, other, tn
  sce.tcr.filt$tp.pheno <- factor(sce.tcr.filt$tp.pheno, levels = c('tp', 'other', 'tn'))
  g <- ggplot(sce.tcr.filt, aes(x = tp.pheno, y = Frequency))
  g <- g + geom_boxplot() + scale_y_log10()
  g <- g + theme_bw() + xlab('Cell phenotype') + ylab('Clonotype frequency')
  g
  ggsave(paste(path.out, '/', sprintf('%s_cloneFreq_3pop.pdf', p), sep = ''), 
         useDingbats = FALSE)
  # p-value
  sink(paste(path.out, '/', sprintf('%s_cloneFreq_3pop.txt', p), sep = ''))
  print(pairwise.wilcox.test(sce.tcr.filt$Frequency,
                             sce.tcr.filt$tp.pheno,
                             p.adjust.method = 'bonferroni'))
  sink()
  
  # number of cells per pheno
  # cannot use df, as it was filtered to get clonotype frequency
  tp.df <- as.data.frame(table(sce.tcr[,grepl(p, sce.tcr$cluster.id)]$tp.pheno))
  colnames(tp.df) <- c('phenotype', 'frequency')
  write.table(tp.df, paste(path.out, '/', sprintf('%s_nbCells_3pop.txt', p), sep = ''), 
              sep = '\t', quote = FALSE)
  
  # all phenotypes
  sce.tcr.filt$sub_ENTPD1.ITGAE.PDCD1. <- factor(sce.tcr.filt$sub_ENTPD1.ITGAE.PDCD1., 
                                                 levels = c('ENTPD1_ITGAE_PDCD1', 
                                                            'ENTPD1_ITGAE',
                                                            'ENTPD1_PDCD1', 'ENTPD1',
                                                            'ITGAE_PDCD1', 'ITGAE', 'PDCD1',
                                                            '-'))
  g <- ggplot(sce.tcr.filt, aes(x = sub_ENTPD1.ITGAE.PDCD1., y = Frequency))
  g <- g + geom_boxplot() + scale_y_log10()
  g <- g + theme_bw() + xlab('Cell phenotype') + ylab('Clonotype frequency')
  g
  ggsave(paste(path.out, '/', sprintf('%s_cloneFreq_all.pdf', p), sep = ''), 
         useDingbats = FALSE)
  # p-value
  sink(paste(path.out, '/', sprintf('%s_cloneFreq_all.txt', p), sep = ''))
  print(pairwise.wilcox.test(sce.tcr.filt$Frequency,
                             sce.tcr.filt$sub_ENTPD1.ITGAE.PDCD1.,
                             p.adjust.method = 'bonferroni'))
  sink()
  # number of cells per pheno
  # cannot use df, as it was filtered to get clonotype frequency
  all.df <- as.data.frame(table(sce.tcr[,grepl(p, 
                                               sce.tcr$cluster.id)]$`sub_ENTPD1+ITGAE+PDCD1+`))
  colnames(all.df) <- c('phenotype', 'frequency')
  write.table(all.df, paste(path.out, '/', sprintf('%s_nbCells_all.txt', p), sep = ''), 
              sep = '\t', quote = FALSE)
}


# TCR overlap analysis. ########################################################################### 
# Generate heatmap presented in Figure 3E and 4E, as well as Supplementary Figure S14A and B
for (p in c('CD8', 'Treg')) {
  # Keep T cells of interest
  sce.tcr.filt <- sce.tcr[,grepl(p, sce.tcr$cluster.id)]
  
  # 3 pop
  combined.tp <- expression2List(sce.tcr.filt, group = "tp.pheno")
  co2 <- clonalOverlap(combined.tp, cloneCall="gene+nt", method = "morisita") + theme_bw()
  co2
  ggsave(paste(path.out, '/', sprintf('%s_cloneOverlapM_tp.pdf', p), sep = ''), 
         useDingbats = FALSE)
  
  # all
  combined.all <- expression2List(sce.tcr.filt, group = "sub_ENTPD1.ITGAE.PDCD1.")
  co2 <- clonalOverlap(combined.all, cloneCall="gene+nt", method = "morisita") + theme_bw()
  co2
  ggsave(paste(path.out, '/', sprintf('%s_cloneOverlapM_all.pdf', p), sep = ''), 
         useDingbats = FALSE)
}
rm(combined.tp, combined.all, co2)

# Export TCR seq for Venn done online with Venny
# https://bioinfogp.cnb.csic.es/tools/venny/
# Generate Venn diagram presented in Figure 3F and 4F
for (p in c('CD8', 'Treg')) {
  # Keep T cells of interest
  sce.tcr.filt <- sce.tcr[,grepl(p, sce.tcr$cluster.id)]
  
  # 3 pop
  for (s in c('tp', 'tn', 'other')) {
    tcr <- unique(colData(sce.tcr.filt[, sce.tcr.filt$tp.pheno == s])$CTstrict)
    write.table(tcr, paste(path.out, '/', sprintf('%s_%s_CTstrict.txt', p, s), sep = ''),
                quote = FALSE, row.names = FALSE)
  }
}
rm(p, s, sce.tcr.filt, tcr)


# Shannon's diversity index analysis. 
# Generate dot plots presented in Figure 3D and 4D, as well as Supplementary Figure S13C,D and S15C,D
# evaluate phenotypic status
# cd39, cd103 and PD1
# c('ENSG00000138185', 'ENSG00000083457', 'ENSG00000188389')
# cd39 +/-
sce.tcr <- idPosAndNegSubsets(sce.tcr, c('ENSG00000138185'))
sce.tcr$`sub_ENTPD1+` <- ifelse(sce.tcr$`sub_ENTPD1+` == '-', 'neg', 'pos')
# cd103 +/-
sce.tcr <- idPosAndNegSubsets(sce.tcr, c('ENSG00000083457'))
sce.tcr$`sub_ITGAE+` <- ifelse(sce.tcr$`sub_ITGAE+` == '-', 'neg', 'pos')
# pd1 +/-
sce.tcr$`sub_PDCD1+` <- ifelse(sce.tcr$`sub_PDCD1+` == '-', 'neg', 'pos')

phenos <- c('sub_PDCD1.', 'sub_ITGAE.', 'sub_ENTPD1.', 'sub_ENTPD1.ITGAE.PDCD1.', 'tp.pheno')
pops <- c('CD8', 'Treg')

for (pheno in phenos) {
  for (p in pops) {
    # Keep T cells of interest
    sce.tcr.filt <- sce.tcr[,grepl(p, sce.tcr$cluster.id)]
    combined.div <- expression2List(sce.tcr.filt, group = c(pheno))
    
    # split by samples
    d.split <- lapply(combined.div, function(x){
      x.list <- split(x, x$Sample)
      return(x.list)
    })
    # flatten list of lists (l[[pheno]][[sample]] --> l[[pheno_sample]])
    res.list <- list()
    n.k <- c()
    k <- 1
    for (i in seq(1, length(d.split))) {
      for (j in seq(1, length(d.split[[i]]))) {
        n.i <- names(d.split)[i]
        n.j <- names(d.split[[i]])[j]
        n.k <- c(n.k, paste(n.i, n.j, sep = '_'))
        # print(c(k, n.k))
        # print(nrow(d.split[[i]][[j]]))
        
        d.split[[i]][[j]]$sample.pheno <- rep(n.k[length(n.k)], nrow(d.split[[i]][[j]]))
        res.list[[k]] <- d.split[[i]][[j]]
        
        k <- k +1
      }
    }
    names(res.list) <- n.k
    
    # compute clonal diversity indexes with scRepertoire function
    r <- clonalDiversity(res.list, group = pheno, exportTable = TRUE)
    r$Shannon <- as.numeric(r$Shannon)
    # plot and export Shannon results: graph + test
    sh.plot <- ggplot(r, aes(x = r[, pheno], y = Shannon))
    sh.plot <- sh.plot + geom_jitter(width = 0.1) + ylim(0, NA) + theme_bw()
    sh.plot <- sh.plot + xlab(sprintf('Cell phenotype: %s', pheno)) + ylab('Shannon\'s diversity index')
    sh.plot
    ggsave(paste(path.out, '/', sprintf('%s_%s_shannonDiversity.pdf', p, pheno), sep = ''), 
           useDingbats = FALSE)
    
    sink(paste(path.out, '/', sprintf('%s_%s_shannonDiversity.txt', p, pheno), sep = ''))
    print(pairwise.wilcox.test(r[, 'Shannon'], r[, pheno], p.adjust.method = 'bonferroni'))
    sink()
  }
}
rm(phenos, pops, sce.tcr.filt, combined.div, d.split, i, j, k, n.k, res.list, r, sh.plot)
rm(n.i, n.j, p, pheno)
