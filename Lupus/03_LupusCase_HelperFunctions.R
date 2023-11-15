prepData <- function(data,
                    splitVar,
                    filter = c("none", "edgeR_joint", "edgeR_sep"),
                    edgeR_filter_spec = list(group = NULL,
                                             min.count = 1,
                                             min.total.count = 0,
                                             large.n = 10,
                                             min.prop = 0)){

  split_ids <- data[[splitVar]]

  ###### no feature-level filtering
  if(filter == "none"){
    sce_per_cl <- lapply(levels(split_ids), function(id){
      data[,split_ids == id]
    })
    names(sce_per_cl) <- levels(split_ids)
    return(sce_per_cl)

    ###### edgeR feature-level filtering on all clusters jointly
  } else if(filter == "edgeR_joint"){
    if(!is.null(edgeR_filter_spec$groupVar)){
      grouping <- data[[groupVar]]
    } else {
      grouping <- NULL
    }
    keep <- filterByExpr(data,
                         group = grouping,
                         min.count = edgeR_filter_spec$min.count,
                         min.total.count = edgeR_filter_spec$min.total.count,
                         large.n = edgeR_filter_spec$large.n,
                         min.prop = edgeR_filter_spec$min.prop)
    data <- data[keep,]
    sce_per_cl <- lapply(levels(split_ids), function(id){
      data[,split_ids == id]
    })
    names(sce_per_cl) <- levels(split_ids)
    return(sce_per_cl)

    ###### edgeR feature-level filtering on all clusters separately
  } else if(filter == "edgeR_sep"){
    sce_per_cl <- lapply(levels(split_ids), function(id){
      data[,split_ids == id]
    })

    sce_per_cl <- lapply(sce_per_cl, function(sce_per_cl_i){

      if(!is.null(edgeR_filter_spec$groupVar)){
        grouping <- data[[groupVar]]
      } else {
        grouping <- NULL
      }

      #TODO large n still computed based on full data, could make it adaptive
      keep <- filterByExpr(sce_per_cl_i,
                           group = grouping,
                           min.count = edgeR_filter_spec$min.count,
                           min.total.count = edgeR_filter_spec$min.total.count,
                           large.n = edgeR_filter_spec$large.n,
                           min.prop = edgeR_filter_spec$min.prop)

      sce_per_cl_i[keep,]
    })

    names(sce_per_cl) <- levels(split_ids)
    return(sce_per_cl)
  }
}

combine_sce <- function(data){

  # if all rownames are the same (none or joint filtering)
  if(length(unique(lapply(data, rownames))) == 1){

    counts <- do.call(cbind, lapply(data, function(element){
      assay(element)}))
    cd <- do.call(rbind, lapply(data, function(element){
      colData(element)}))
    rd <- rowData(data[[1]])
    if(length(data) > 1){
      for(i in 2:length(data)){
        cn <- tail(colnames(rowData(data[[i]])),1)
        rd[[cn]]  <- rowData(data[[i]])[[cn]]
      }
    }

    sce_combined <- SingleCellExperiment(list(counts = counts),
                                         colData = cd,
                                         rowData = rd)

    narm <- is.na(rowSums(counts))
    if(sum(narm > 0)){
      warning("Removing ", sum(narm), " features with NA rowsums")
      sce_combined <- sce_combined[!narm,]
    }

    return(sce_combined)
  } else { # with separate filtering per cell type

    # get assay
    allFeatures <- unique(unlist(lapply(data, rownames)))
    allCols <- unique(unlist(lapply(data, colnames)))
    counts <- matrix(data = 0, nrow=length(allFeatures), ncol = length(allCols),
                     dimnames = list(allFeatures, allCols))
    for(i in seq_along(data)){
      adata <- assay(data[[i]])
      counts[match(rownames(adata),allFeatures),
             match(colnames(adata),allCols)] <- as.matrix(adata)
    }

    # get coldata
    cd <- do.call(rbind, lapply(data, function(element){
      colData(element)}))

    # TODO make less verbose
    # get rowdata
    rd <- data.frame(features = allFeatures)
    rownames(rd) <- rd$features
    cn_common <- colnames(rowData(data[[1]]))[1:ncol(rowData(data[[1]]))-1]
    rd[cn_common] <- NA
    for(i in seq_along(data)){
      rd[match(rownames(data[[i]]), rownames(rd)),cn_common] <- as.data.frame(rowData(data[[i]])[,cn_common])

      altered <- tail(colnames(rowData(data[[i]])),1)
      rd[[altered]] <- rowData(data[[i]])[match(rownames(rd),rownames(data[[i]])),altered]

      rd[[altered]][is.na(rd[[altered]])] <- FALSE
    }
    rd$features <- NULL

    sce_combined <- SingleCellExperiment(list(counts = counts),
                                         colData = cd,
                                         rowData = rd)

    narm <- is.na(rowSums(counts))
    if(sum(narm > 0)){
      warning("Removing ", sum(narm), " features with NA rowsums")
      sce_combined <- sce_combined[!narm,]
    }

    return(sce_combined)
  }
}

violin_DGE <- function(pb, genes){
    cd <- colData(pb)
    libSize <- pb %>%
        counts %>%
        colSums %>%
        unname

    gg_list <- lapply(genes, function(gene){

        data <- data.frame(cpm = log2(counts(pb)[gene,]+0.5) - log2(libSize+1)+log2(1e6),
                           status = cd$SLE_status,
                           batch = cd$batch_cov)
        data$cpm[is.infinite(data$cpm)] <- NA

        data$status_batch <- factor(paste0(data$status, "_", data$batch))
        data$dotposition <- data$status_batch
        levels(data$dotposition) <- list("0.65" = levels(data$status_batch)[1],
                                         "0.79" = levels(data$status_batch)[2],
                                         "0.93" = levels(data$status_batch)[3],
                                         "1.07" = levels(data$status_batch)[4],
                                         "1.21" = levels(data$status_batch)[5],
                                         "1.35" = levels(data$status_batch)[6],
                                         "1.65" = levels(data$status_batch)[7],
                                         "1.79" = levels(data$status_batch)[8],
                                         "1.93" = levels(data$status_batch)[9],
                                         "2.07" = levels(data$status_batch)[10],
                                         "2.21" = levels(data$status_batch)[11],
                                         "2.35" = levels(data$status_batch)[12])
        data$dotposition <- as.numeric(as.character(data$dotposition))

        gg <- ggplot(data = data, aes(x=status,y=cpm)) +
            geom_violin() +
            geom_jitter(aes(x=dotposition,
                            y=cpm,
                            col=batch),
                        size=2,
                        width = 0.03) +
            theme_bw() +
            theme(plot.title = element_text(size=11))
        return(gg)
    })
    return(gg_list)
}

violin_DD <- function(pb, genes){
    bin_counts <- assay(pb)
    cd <- colData(pb)
    of <- colMeans(sweep(bin_counts, 2, cd$ncells, "/"))
    logprop <- sweep(log(bin_counts, base = 2),
                     2,
                     (log(cd$ncells, base = 2) + log(of, base = 2)),
                     "-")

    gg_list <- lapply(genes, function(gene){

        data <- data.frame(logprop = logprop[gene,],
                           status = cd$SLE_status,
                           batch = cd$batch_cov)
        data$logprop[is.infinite(data$logprop)] <- NA

        data$status_batch <- factor(paste0(data$status, "_", data$batch))
        data$dotposition <- data$status_batch
        levels(data$dotposition) <- list("0.65" = levels(data$status_batch)[1],
                                         "0.79" = levels(data$status_batch)[2],
                                         "0.93" = levels(data$status_batch)[3],
                                         "1.07" = levels(data$status_batch)[4],
                                         "1.21" = levels(data$status_batch)[5],
                                         "1.35" = levels(data$status_batch)[6],
                                         "1.65" = levels(data$status_batch)[7],
                                         "1.79" = levels(data$status_batch)[8],
                                         "1.93" = levels(data$status_batch)[9],
                                         "2.07" = levels(data$status_batch)[10],
                                         "2.21" = levels(data$status_batch)[11],
                                         "2.35" = levels(data$status_batch)[12])
        data$dotposition <- as.numeric(as.character(data$dotposition))

        gg <- ggplot(data = data,aes(x=status, y= logprop)) +
            geom_violin() +
            geom_jitter(aes(x=dotposition,
                            y=logprop,
                            col=batch),
                        size=2,
                        width = 0.05) +
            theme_bw() +
            theme(plot.title = element_text(size=11))
        return(gg)
    })
    return(gg_list)
}


## Tidy visualizations (violins and scatter) for manuscript
scatter_violin_paper <- function(scatterplot,
                                 scatterplotData,
                                 violinplots_DD,
                                 violinplots_DE,
                                 genes,
                                 legend_position,
                                 nudge_x,
                                 nudge_y){

    # scatter plot of -log10 DD p-values against -log10 DE p-values
    # topleft panel of full figure
    topleft <- scatterplot +
        theme(strip.text.x = element_blank(),
              legend.position = legend_position,
              legend.title = element_text(size = 10),
              legend.text = element_text(size = 9),
              legend.spacing.y = unit(2, "mm"),
              legend.key.size = unit(0.5, 'lines')
        ) +
        ggrepel::geom_label_repel(
            data = scatterplotData %>% filter(gene %in% genes),
            aes(label=gene),
            size = 3,
            fill = "white",
            color = "black",
            nudge_x = nudge_x,
            nudge_y = nudge_y,
            box.padding   = 1.5,
            point.padding = 0.5,
            segment.color = 'grey50',
            min.segment.length = 0)  +
        ggtitle("Panel A: pval-pval plot") +
        theme(plot.title = element_text(size = 11, face="bold"),
              plot.margin = margin(10, 50, 10, 10, "pt"))

    # Example of a gene that is both DE and DD (topright panel)
    topright <- ((violinplots_DE[[1]] +
                      ggtitle(paste0("Panel B: Gene ", genes[1])) +
                      theme(plot.title = element_text(size = 11, face="bold"))) +
                     violinplots_DD[[1]]) +
        theme(plot.margin = margin(10, 10, 10, 10, "pt")) +
        plot_layout(guides = "collect")

    # Example of a gene that is DE but not DD (bottomleft panel)
    bottomleft <- ((violinplots_DE[[2]] +
                        ggtitle(paste0("Panel C: Gene ", genes[2])) +
                        theme(plot.title = element_text(size = 11, face="bold"))) +
                       violinplots_DD[[2]]) +
        theme(plot.margin = margin(30, 50, 10, 10, "pt"))

    # Example of a gene that is DD but not DE (bottomright panel)
    bottomright <- ((violinplots_DE[[3]] +
                         ggtitle(paste0("Panel D: Gene ", genes[3])) +
                         theme(plot.title = element_text(size = 11, face="bold"))) +
                        violinplots_DD[[3]]) +
        theme(plot.margin = margin(30, 10, 10, 10, "pt"))

    return(list(topleft = topleft,
                topright = topright,
                bottomleft = bottomleft,
                bottomright = bottomright))
}
