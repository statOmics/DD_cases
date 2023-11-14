# GSEA helper function
run_gsea <- function(genes,
                     geneInfo,
                     universe){
    genes <- genes[!is.na(genes)]
    names(genes) <- geneInfo$entrezgene_id[match(names(genes),
                                                 geneInfo$hgnc_symbol)]
    genes <- names(genes[genes<0.05])
    gsea <- goana(de = genes,
                  universe = universe,
                  species = "Hs",
                  prior.prob = NULL,
                  covariate = NULL,
                  plot = FALSE)
    gsea <- topGO(gsea, ontology="BP", number = Inf)
    gsea$Rank <- 1:nrow(gsea)
    return(gsea)
}

# Visualization helper functions

## Wrangle data for creating scatterplots
scatterplot_data <- function(data_bin,
                             data_count,
                             DD_res,
                             DE_res,
                             stageR_res,
                             group1,
                             group2){

    # Differential detection data
    select <- data_bin$Status_on_day_collection_summary %in% c(group1, group2)
    bin_counts <- assay(data_bin[,select])
    cd <- colData(data_bin[,select])
    cd <- droplevels(cd)
    detection <- rowMeans(sweep(bin_counts, 2, cd$ncells, "/"))

    # Differential expression data
    select <- data_count$Status_on_day_collection_summary %in% c(group1, group2)
    counts <- assay(data_count[,select])
    libSize <- counts %>%
        colSums %>%
        unname
    expression <- rowMeans(sweep(log2(counts+0.5), 2, (-log2(libSize+1) + log2(1e6)), "+"))

    # Combine data summary statistics with results
    plotData <- data.frame(gene = rownames(DD_res),
                           DD_pval = DD_res$PValue,
                           DD_FDR = p.adjust(DD_res$PValue,
                                             method = "BH"),
                           DE_pval = DE_res$PValue,
                           DE_FDR = p.adjust(DE_res$PValue,
                                             method = "BH"),
                           stage_FDR = stageR_res[,1],
                           detection = detection,
                           logCPM = expression)
    plotData$DD_FDR[is.na(plotData$DD_FDR)] <- 1
    plotData$Category <- case_when(plotData$DD_FDR >= 0.05 & plotData$DE_FDR >= 0.05 ~ "Neither",
                                   plotData$DD_FDR < 0.05 & plotData$DE_FDR >= 0.05 ~ "DD only",
                                   plotData$DD_FDR >= 0.05 & plotData$DE_FDR < 0.05 ~ "DE only",
                                   plotData$DD_FDR < 0.05 & plotData$DE_FDR < 0.05 ~ "Both")
    plotData$Stagewise <- ifelse(plotData$stage_FDR <= 0.05, "yes", "no")
    return(plotData)
}

## Create scatterplots
scatterplots_visualize <- function(plotData){



    # Scatterplot 1: -log10 of DGE p-values against logCPM (normalized
    # expression, averaged across pseudobulk samples)
    plot1 <- ggplot(data = plotData,
                    aes(x = logCPM,
                        y = -log10(DE_pval),
                        color = Category,
                        fill = Stagewise,
                        shape = Stagewise)) +
        geom_point(size = 0.9,stroke=0.8) +
        scale_shape_manual(values=c(19,21)) +
        scale_color_manual(values=c("Both" = "#009E73",
                                    "DD only" = "#CC79A7",
                                    "DE only" = "#D55E00",
                                    "Neither" = "black")) +
        scale_fill_manual(values=c(NA,"#F0E442"),na.translate = FALSE) +
        theme_bw()

    # Scatterplot 2: -log10 DGE p-values against -log10 of DD p-values
    plot2 <- ggplot(data = plotData,
                    aes(x = -log10(DD_pval),
                        y = -log10(DE_pval),
                        color = Category,
                        fill = Stagewise,
                        shape = Stagewise)) +
        geom_point(size = 0.9,stroke=0.8) +
        scale_shape_manual(values=c(19,21)) +
        scale_color_manual(values=c("Both" = "#009E73",
                                    "DD only" = "#CC79A7",
                                    "DE only" = "#D55E00",
                                    "Neither" = "black")) +
        scale_fill_manual(values=c(NA,"#F0E442"),na.translate = FALSE) +
        theme_bw()

    # Scatterplot 3: logCPM (normalized expression, averaged across pseudobulk
    # samples) against average fraction of gene detection across pseudobulk
    # samples
    plot3 <- ggplot(data = plotData,
                    aes(x = logCPM,
                        y = detection,
                        color = Category,
                        fill = Stagewise,
                        shape = Stagewise)) +
        geom_point(size = 0.9,stroke=0.8) +
        scale_shape_manual(values=c(19,21)) +
        scale_color_manual(values=c("Both" = "#009E73",
                                    "DD only" = "#CC79A7",
                                    "DE only" = "#D55E00",
                                    "Neither" = "black")) +
        scale_fill_manual(values=c(NA,"#F0E442"),na.translate = FALSE) +
        theme_bw()

    # Scatterplot 1: -log10 of DD p-values against average fraction of gene
    # detection across pseudobulk samples
    plot4 <- ggplot(data = plotData,
                    aes(x = detection,
                        y = -log10(DD_pval),
                        color = Category,
                        fill = Stagewise,
                        shape = Stagewise)) +
        geom_point(size = 0.9,stroke=0.8) +
        scale_shape_manual(values=c(19,21)) +
        scale_color_manual(values=c("Both" = "#009E73",
                                    "DD only" = "#CC79A7",
                                    "DE only" = "#D55E00",
                                    "Neither" = "black")) +
        scale_fill_manual(values=c(NA,"#F0E442"),na.translate = FALSE) +
        theme_bw()

    return(list(plot1, plot2, plot3, plot4))
}

## Violin plots for DE (logCPM)
violin_DGE <- function(data_count,
                       group1,
                       group2,
                       genes){

    # filter binary data to requested groups
    select <- which(data_count$Status_on_day_collection_summary %in% c(group1, group2))
    pb_count <- data_count[,select]
    colData(pb_count) <- droplevels(colData(pb_count))

    # compute library size
    libSize <- pb_count %>%
        counts %>%
        colSums %>%
        unname

    gg_list <- lapply(genes, function(gene){

        data <- data.frame(logCPM = log2(counts(pb_count)[gene,]+0.5) - log2(libSize+1)+log2(1e6),
                           status = pb_count$Status_on_day_collection_summary,
                           batch = paste0(pb_count$Site,"_",pb_count$sex))
        data$logCPM[is.infinite(data$logCPM)] <- NA

        levels(data$batch) <- c("Cambridge_female", "Cambridge_male",
                                "Ncl_female", "Ncl_male")

        data$status_batch <- factor(paste0(data$status, "_", data$batch),
                                    levels = paste0(rep(levels(data$status),each=4),
                                                    "_",
                                                    rep(levels(data$batch), times=2)))
        data$dotposition <- data$status_batch
        levels(data$dotposition) <- list("0.7" = levels(data$status_batch)[1],
                                         "0.9" = levels(data$status_batch)[2],
                                         "1.1" = levels(data$status_batch)[3],
                                         "1.3" = levels(data$status_batch)[4],
                                         "1.7" = levels(data$status_batch)[5],
                                         "1.9" = levels(data$status_batch)[6],
                                         "2.1" = levels(data$status_batch)[7],
                                         "2.3" = levels(data$status_batch)[8])
        data$dotposition <- as.numeric(as.character(data$dotposition))

        gg <- ggplot(data = data, aes(x=status,y=logCPM)) +
            geom_violin() +
            geom_jitter(aes(x=dotposition,
                            y=logCPM,
                            col=batch),
                        size=2,
                        width = 0.05) +
            scale_color_manual(values=c("#E69F00", "#56B4E9", "#009E73", "#CC79A7")) +
            theme_bw() +
            theme(plot.title = element_text(size=11))
        return(gg)
    })
    return(gg_list)
}

## Violin plots for DD (normalized, log-transformed proportion)
violin_DD <- function(data_bin,
                      group1,
                      group2,
                      genes,
                      normalize = TRUE){

    # filter binary data to requested groups
    select <- which(data_bin$Status_on_day_collection_summary %in% c(group1, group2))
    pb_bin <- data_bin[,select]

    bin_counts <- assay(pb_bin)
    colData(pb_bin) <- droplevels(colData(pb_bin))

    # if normalize = TRUE, compute normalized log-proportions
    if(normalize){
        of <- colMeans(sweep(bin_counts, 2, pb_bin$ncells, "/")) # CDR normalization offset
        logprop <- sweep(log(bin_counts, base = 2),
                         2,
                         (log(pb_bin$ncells, base = 2) + log(of, base = 2)),
                         "-")
    } else {
        logprop <- log(sweep(bin_counts, 2, pb_bin$ncells, "/"), base = 2)
    }


    gg_list <- lapply(genes, function(gene){

        data <- data.frame(logprop = logprop[gene,],
                           status = pb_bin$Status_on_day_collection_summary,
                           batch = paste0(pb_bin$Site,"_",pb_bin$sex))
        data$logprop[is.infinite(data$logprop)] <- NA

        levels(data$batch) <- c("Cambridge_female", "Cambridge_male",
                                "Ncl_female", "Ncl_male")

        data$status_batch <- factor(paste0(data$status, "_", data$batch),
                                    levels = paste0(rep(levels(data$status),each=4),
                                                    "_",
                                                    rep(levels(data$batch), times=2)))
        data$dotposition <- data$status_batch
        levels(data$dotposition) <- list("0.7" = levels(data$status_batch)[1],
                                         "0.9" = levels(data$status_batch)[2],
                                         "1.1" = levels(data$status_batch)[3],
                                         "1.3" = levels(data$status_batch)[4],
                                         "1.7" = levels(data$status_batch)[5],
                                         "1.9" = levels(data$status_batch)[6],
                                         "2.1" = levels(data$status_batch)[7],
                                         "2.3" = levels(data$status_batch)[8])
        data$dotposition <- as.numeric(as.character(data$dotposition))

        gg <- ggplot(data = data,aes(x=status, y= logprop)) +
            geom_violin() +
            geom_jitter(aes(x=dotposition,
                            y=logprop,
                            col=batch),
                        size=2,
                        width = 0.05) +
            scale_color_manual(values=c("#E69F00", "#56B4E9", "#009E73", "#CC79A7")) +
            theme_bw() +
            theme(plot.title = element_text(size=11))
        return(gg)
    })
    return(gg_list)
}

## Tidy visualizations (violins and scatter) for manuscript
scatter_violin_paper <- function(scatterplots,
                                 scatterplotData,
                                 violinplots_DD,
                                 violinplots_DE,
                                 genes,
                                 legend_position,
                                 nudge_x,
                                 nudge_y){

    # scatter plot of -log10 DD p-values against -log10 DE p-values
    # topleft panel of full figure
    scatter <- scatterplots[[2]]
    topleft <- scatter +
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

