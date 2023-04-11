#libraries
library(ggplot2)
library(Matrix)
library(DropletUtils)
library(uwot)
library(igraph)
library(reshape2)
library(scran)

#functions
assembleCounts <- function(sce){
  counts <- counts(sce)
  rownames(counts) <- rowData(sce)$Symbol
  colnames(counts) <- sce$Barcode
  rownames(counts) <- make.names(rownames(counts), unique = TRUE)
  colnames(counts) <- make.names(colnames(counts), unique = TRUE)
  return(counts)
}

assembleMeta <- function(sce, labels){
  meta <- data.frame(cell = sce$Barcode)
  meta$cell <-  make.names(meta$cell, unique = TRUE)
  meta$sample <- sce$Sample
  return(meta)
}

lowQualityCells <- function(counts){
  br.out <- barcodeRanks(counts)
  o <- order(br.out$rank)
  m <- which(br.out$total > metadata(br.out)$knee)
  return(m)
}

normalize <- function(counts, pseudocount = 1e4){
  librarySize <- Matrix::colSums(counts)
  return(pseudocount * t(t(counts) / librarySize ))
}

compute_hvg <- function(counts.normalized, lowmean = 0.05, highmean = 0.8, fano = 0.65){
  mean.expression <- Matrix::rowMeans(counts.normalized)
  meansq.expression <- Matrix::rowMeans(counts.normalized^2)
  fano.factor <- (meansq.expression - (mean.expression^2))/mean.expression
  hvg <- counts.normalized[which((mean.expression >= quantile(mean.expression, lowmean)) & (fano.factor >= quantile(na.omit(fano.factor),fano)) & (mean.expression <= quantile(mean.expression,highmean))),]
  return(hvg)
}

umap_project <- function(meta, hvg, neighbors = 10, dist = 0.3, seed = 42, metric="cosine",pca=NULL,threads=16){
  X <- t(as.matrix(log2(hvg+1)))
  set.seed(seed)
  out <- umap(X, n_neighbors = neighbors, metric = metric, min_dist = dist, verbose = TRUE, n_threads = threads, n_sgd_threads = "auto", batch = TRUE, pca = pca)
  meta$x <- out[,1]
  meta$y <- out[,2]
  return(meta)
}




umap_cluster <- function(meta, hvg, neighbors = 10, steps = 4, dist = 0.3, seed1 = 42, seed2 = 42, metric="cosine",pca=NULL,threads=16){
  X <- t(as.matrix(log2(hvg+1)))
  set.seed(seed1)
  out <- umap(X, n_neighbors = neighbors, metric = metric, min_dist = dist, verbose = TRUE, n_threads = threads, n_sgd_threads = "auto", batch = TRUE, pca = pca, n_epochs = 0,ret_extra = "fgraph")
  g <- graph_from_adjacency_matrix(out$fgraph,weighted = TRUE,mode = "undirected")
  set.seed(seed2)
  meta$cluster <- cluster_walktrap(g,steps = steps)$membership
  return(meta)
}

#plotting functions
plotMeta <- function(meta, mode = "none"){
  if(mode == "none"){
    p <- ggplot( data = meta, mapping= aes(x = x, y = y)) +
      theme_bw() +
      geom_point(size=0.1, col = "black") +
      coord_equal() +
      theme(aspect.ratio = 1)+
      theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),aspect.ratio=1) 
  }
  
  else if(mode == "cluster"){
    centroid <- aggregate(cbind(x,y) ~ cluster, data=meta, FUN=median)
    p <- ggplot( data = meta, mapping= aes(x = x, y = y, col = factor(meta$cluster))) +
      theme_bw() +
      geom_point(size=0.1) +
      coord_equal() +
      theme(aspect.ratio = 1)+
      theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),aspect.ratio=1)  +
      theme(legend.position = "none") +
      geom_text(data = centroid, mapping=aes(x = x, y = y, label = cluster), col = "black")
    
  }
  else if(mode == "cycling"){
    p <- ggplot( data = meta, mapping= aes(x = x, y = y, col = factor(meta$CellCyclePhase))) +
      theme_bw() +
      geom_point(size=0.1) +
      coord_equal() +
      theme(aspect.ratio = 1)+
      theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),aspect.ratio=1)  
    
    
  }
  
  return(p)
}


plotGene <- function(meta, countn, geneName, lims = "none", size = 0.2, title = "none"){
  
  X <- meta$x
  Y <- meta$y
  atlas <- as.data.frame(cbind(X,Y))
  if(!(geneName %in% rownames(countn))){
    return()
  }
  gene <- countn[which(rownames(countn) == geneName),]
  if(max(gene) == 0){ p <- ggplot( data = atlas, mapping=aes(x=X,y=Y)) +
    geom_point(size=size,col=alpha("grey",1)) +
    theme_bw() +
    # scale_colour_gradient2(low="#ff6666", mid = "#e60000", high = "#990000") + 
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    
    theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", aspect.ratio = 1) +
    ggtitle(geneName)
  return(p)}
  log_counts = log10(gene+1)
  Xsub <- X[log_counts>0]
  Ysub <- Y[log_counts>0]
  log_counts = log_counts[log_counts>0]
  m <- order(log_counts)
  Xsub <- Xsub[m]
  Ysub <- Ysub[m]
  log_counts <- log_counts[m]
  subset <- as.data.frame(cbind(Xsub, Ysub))
  p <- ggplot( data = atlas, mapping=aes(x=X,y=Y)) +
    geom_point(size=size,col=alpha("grey",1)) +
    geom_point(data = subset,mapping = aes(x=Xsub, y=Ysub, col = log_counts), size = size) +
    theme_bw() +
    # scale_colour_gradient2(low="#ff6666", mid = "#e60000", high = "#990000") + 
    scale_colour_gradient(low="yellow", high = "red") + 
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    coord_equal() + 
    
    theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
    if (title == "none"){
      ggtitle(geneName)
    }
    else{
      ggtitle(title)
    }
  

  if(lims != "none"){
    p <- p + xlim(lims[1:2]) + ylim(lims[3:4])
  }
  
  return(p)
  
}

plotSample <- function(meta,m1, title = "none", mode = "none", size = 0.3, lims = "none", return = FALSE){
  umap <- as.data.frame(cbind(meta$x,meta$y))
  umap1 <- as.data.frame(cbind(meta$x[m1],meta$y[m1]))
  umap1$V3 <- meta$sample[m1]
  if(mode == "none"){ cols <- "red"}
  else if(mode == "batch"){ cols <- factor(umap1$V3)}
  p <- ggplot( data = umap, mapping= aes(x = V1, y = V2)) +
    theme_bw() +
    geom_point(size=size, col = "grey", stroke = 0) +
    scale_shape(solid = F) +
    coord_equal() +
    geom_point(data = umap1, mapping = aes(x = as.numeric(V1), y = as.numeric(V2), col = cols),size=size, fill = NA, stroke = 0) +
    theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
    guides(colour = guide_legend(override.aes = list(size=3,
                                                     alpha = 1))) + theme(legend.position = "none")
  if(title != "none"){
    p <- p + ggtitle(title)
  }
  if(lims != "none"){
    p <- p + xlim(lims[1:2]) + ylim(lims[3:4])
  }
  if(return == TRUE){
    return(p)
  }
  else{
    print(p)
  }
}

#Construct a heatmap for the average expression of different genes across different conditions/clusters
dotplot <- function(countn, genes = rownames(countn), cells = 1:length(colnames(countn)), condition, norm = "max", conditionList = "none", aspect.ratio = 1, xAngle = 45, yAngle = 0, xSize = 10, ySize = 10, xAdj = 1, plot = TRUE, 
                    title = "title", return = FALSE, normcells = 1:length(colnames(countn)), dot.scale = 2, collow = "lightgrey", colhigh = "blue"){
  
  
  genes <- rev(genes)
  
  
  subset <- log10(countn[na.omit(match(genes,rownames(countn))),cells]+1)
  nonzero <- countn[na.omit(match(genes,rownames(countn))),cells] > 0
  
  avgd <- (t(aggregate(t(as.matrix(subset)),list(condition[cells]),mean)))
  clusters <- avgd[1,]
  avgd <- avgd[-1,]
  avgd <- matrix(as.numeric(unlist(avgd)),nrow=nrow(avgd))
  rownames(avgd) <- rownames(subset)
  colnames(avgd) <- clusters
  
  pct <- (t(aggregate(t(as.matrix(nonzero)),list(condition[cells]),mean)))
  clusters <- pct[1,]
  pct <- pct[-1,]
  pct <- matrix(as.numeric(unlist(pct)),nrow=nrow(pct))
  rownames(pct) <- rownames(nonzero)
  colnames(pct) <- clusters
  
  
  
  if (norm == "max"){
    avgd_n <- (avgd / apply(avgd,1,function(x) max(x)))
  } else if (norm == "sum"){
    avgd_n <- avgd / rowSums(avgd)
  } 
  else if (norm == "totalMax"){
    subsetMax <- log10(countn[na.omit(match(genes,rownames(countn))),normcells]+1)
    avgdMax <- (t(aggregate(t(as.matrix(subsetMax)),list(condition[normcells]),mean)))
    clustersMax <- avgdMax[1,]
    avgdMax <- avgdMax[-1,]
    avgdMax <- matrix(as.numeric(unlist(avgdMax)),nrow=nrow(avgdMax))
    rownames(avgdMax) <- rownames(subsetMax)
    colnames(avgdMax) <- clustersMax
    avgd_n <- (avgd / apply(avgdMax,1,function(x) max(x)))
    
  } else {    avgd_n <- avgd
  }
  
  if(conditionList != "none"){
    avgd_n <- avgd_n[,match(conditionList,colnames(avgd_n))]
    pct <- pct[,match(conditionList,colnames(pct))]
    
  }
  
  if(length(which(is.na(rowSums(avgd_n)))) > 0){
    pct <- pct[-which(is.na(rowSums(avgd_n))),]
    avgd_n <- avgd_n[-which(is.na(rowSums(avgd_n))),]
  }
  
  melted_ <- melt(t(avgd_n))
  pct_ <- melt(t(pct))
  melted_$pct <- 100*pct_$value
  
  p <- ggplot(data = melted_) + 
    theme_bw() +
    geom_point(mapping = aes(x = Var1, y = Var2,size = pct,col = value))+
    scale_color_gradient(low = collow,high = colhigh)+
    scale_size(range = c(0,dot.scale), limits = c(0, 100)) +
    theme(axis.text.x = element_text(angle = xAngle, hjust = xAdj, size = xSize), axis.text.y = element_text(angle = yAngle, size = ySize), aspect.ratio = aspect.ratio, axis.title.x = element_blank(), axis.title.y = element_blank()) +
    theme(panel.background = element_blank()) + theme(panel.grid = element_blank())
  
  p <- p + ggtitle(title)
  if(plot == TRUE){
    print(p)
  }
  if (return) {return(avgd_n)}
}



savePDF <- function(filename,scale = 1, width = 5, height = 5, dpi = 300){
  ggsave(paste0(filename,".pdf"), plot = last_plot(), device = "pdf", scale = scale, width = width, height = height, units = "in", dpi = dpi,  useDingbats = FALSE)
}

saveTIFF <- function(filename,scale = 1, width = 5, height = 5, dpi = 300){
  ggsave(paste0(filename,".tiff"), plot = last_plot(), device = "tiff", scale = scale, width = width, height = height, units = "in", dpi = dpi)
}

saveFig <- function(name = "name",mode = "pdf",scale = 1, width = 5, height = 5, dpi = 300){ 
  name = paste0("panels/", name)
  if(plot){
    if(mode == "pdf"){savePDF(name,scale = 1, width = 5, height = 5, dpi = 300)}
    else if (mode == "tiff"){saveTIFF(name,scale = 1, width = 5, height = 5, dpi = 300)}
  }
}


