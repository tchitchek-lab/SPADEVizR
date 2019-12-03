#' @title Computing of the fraction of clusters with low number of associated cells
#' 
#' @description 
#' Computes the fraction of clusters having a number of associated cells less than a specific threshold (th.size parameter).
#'
#' @details 
#' This function can also generate a PDF file summarizing the results.
#'
#' @param Results a Results object
#' @param clusters a character vector containing the clusters names to be computed (by default all clusters will be computed)
#' @param th.size a numeric value specifying the minimum number of cells needed for a cluster to be considered a a small cluster
#' @param PDFfile a character specifying the location of the output PDF path
#' @param width a numeric specifying the plot width in the output PDF file
#' @param height a numeric specifying the plot height in the output PDF file
#' @param tile.color a character specifying the border color of the tiles (NA to remove tile borders)
#'  
#' @return a list containing a numeric (perc numeric element) containing the fraction of clusters having a number of associated cells less than the threshold and a dataframe (small.clusters element) specifying the clusters having a number of associated cells less than the threshold
#'
#' @export
qcSmallClusters <- function(Results,
                            clusters   = NULL,
                            th.size    = 50,
                            PDFfile    = "qcSmallClusters.pdf",
                            width      = ncol(Results@cluster.abundances),
                            height     = nrow(Results@cluster.abundances)/4,
							tile.color = "black"){
    message("[BEGIN] - generating qc small clusters")
    
    data <- Results@cluster.abundances
    
    if (is.null(clusters)) {
        clusters <- Results@cluster.names
    } else if (all(clusters %in% Results@cluster.names)) {
        if (typeof(clusters) != "character") {
            stop("Error in qcSmallClusters: 'clusters' parameter must be a character vector")
        }
        clusters <- unique(clusters)
        data <- data[clusters, , drop = FALSE]
    } else {
        stop("Error in qcSmallClusters:\nUnknown clusters : ", paste(setdiff(unique(clusters), Results@cluster.names), collapse = " "))
    }
    
    total_cells <- apply(data, 1, FUN = function(x){ifelse((sum(x) < th.size), TRUE, FALSE)})
    
    data[data < th.size]  <- TRUE
    data[data >= th.size] <- FALSE
    
    data           <- cbind(data, total_cells = total_cells)
    data           <- apply(data, 2, as.logical)
    rownames(data) <- clusters
    perc           <- sum(data[, "total_cells"] == TRUE)/nrow(data) * 100
    
    data.melted           <- reshape2::melt(data)
    colnames(data.melted) <- c("cluster", "sample", "small")
    data.melted$cluster   <- factor(data.melted$cluster, levels = rev(unique(data.melted$cluster)))
    
    on.exit(dev.off())
    pdf(PDFfile, width = width, height = height)
    
    title    <- paste0("QC: SmallClusters")
    subtitle <- paste0("percentage of clusters having a small number of cells (<", th.size, " cells) = ", 
                       format(round(perc, 2), nsmall = 2), "%")
    
    plot <- ggplot2::ggplot(data = data.melted) +
        ggplot2::ggtitle(bquote(atop(.(title), atop(italic(.(subtitle)), "")))) +
        ggplot2::geom_tile(ggplot2::aes_string(x = "sample", y = "cluster", fill = "small"), colour = tile.color) +            
        ggplot2::scale_x_discrete(expand = c(0, 0)) + 
        ggplot2::scale_y_discrete(expand = c(0, 0)) +
        ggplot2::scale_fill_manual(values = c("green", "red")) +
        ggplot2::geom_vline(xintercept = (ncol(data) + 0.5), colour = "black", size = 2) +
        ggplot2::geom_vline(xintercept = (ncol(data) - 0.5), colour = "black", size = 1) +
        ggplot2::geom_vline(xintercept = 0.5, colour = "black", size = 1) +
        ggplot2::geom_hline(yintercept = 0.5, colour = "black", size = 1) +
        ggplot2::geom_hline(yintercept = length(Results@cluster.names) + 0.5, colour = "black", size = 1) +
        ggplot2::xlab("samples") +
        ggplot2::ylab("clusters") +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x      = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
                       legend.key       = ggplot2::element_blank(),
                       axis.line        = ggplot2::element_blank(),
                       panel.background = ggplot2::element_blank(),
                       panel.border     = ggplot2::element_blank(),
                       panel.grid.major = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank(),
                       plot.background  = ggplot2::element_blank())
    plot(plot)
    
    message("[END] - generating qc small clusters")
    invisible(list(perc = perc, small.clusters = data))
}

#' @title Computing of the fraction of cluster with uniform phenotypes
#' 
#' @description 
#' The `qcUniformClusters()` function allows to report as PDF files 2 kinds of information. 
#' Firstly, the identification of clusters having non-unimodal marker expression densities is performs using a Hartigan's dip test. 
#' Secondly, markers of cluster having a large spread of expressions (above a threshold) are identified.
#' 
#' If all the clustering marker expressions of a cluster checks this assessments (depending of the 'uniform.test' parameter), the cluster is considered as uniform.
#'
#' @details
#' The test of multimodal distribution is performs using a Hartigan's dip test. 
#' The null hypothesis of this test assume the distribution as unimodal 
#' Supplemental parameters passed to the dip.test function are ('simulate.p.value' and 'B').
#' 'simulate.p.value' is a logical indicating whether to compute p-values by Monte Carlo simulation (FALSE by default).
#' 'B' is an integer specifying the number of replicates used in the Monte Carlo test (2000 by default).
#' 
#' The spread of a marker expression is considered as high if the IQR (interquartile range) is higher than a provided threshold th.IQR.
#' 
#' The 'uniform.test' parameter can take 3 values: 
#' "unimodality": check if the distribution of the marker expression is unimodal or not using a Hartigans dip test
#' "spread": check if the spread of the marker expression is below an IQR threshold (interquantile range)
#' "both": check "unimodality" and "spread"
#'
#' @param Results a Results object (with 'flowset' slot not null)
#' @param clusters a character vector containing the clusters names to be computed (by default all clusters will be computed)
#' @param only.clustering_markers a logical specifying if only clustering marker densities must be displayed.
#' @param uniform.test a character specifying the qc assessment to perform
#' @param th.pvalue a numeric value specifying the pvalue threshold of the Hartigan's dip test (multimodal if p.value < th.pvalue)
#' @param th.IQR a numeric value specifying the IQR (interquartile range) threshold to assume a distribution as uniform
#' @param density.PDFfile a character specifying the output path of the marker expression density plots (or NULL to avoid this step)
#' @param density.PDFfile.dim a numeric vector specifying the width and the height of the density PDF file
#' @param heatmap.PDFfile a character specifying the output path of the marker expression accuracy heatmap (or NULL to avoid this step)
#' @param heatmap.PDFfile.dim a numeric vector specifying the width and the height of the heatmap PDF file
#' @param tile.color a character specifying the border color of the tiles (NA to remove tile borders)
#' @param verbose a boolean specifying if some messages must be displayed during the generation of the qc report
#' @param ... supplemental parameters passed to the dip.test function
#' 
#' @return a list containing a numeric (perc numeric element) specifying the percentage of clusters having only uniform clustering marker expressions and a dataframe (accuracy.matrix element) specifying if each marker of each cluster is uniform or not.
#' 
#' @export
#'
#' @import diptest
#' @importFrom gtools mixedsort
#' @importFrom flowCore fsApply
qcUniformClusters <- function(Results,
                              clusters                = NULL,
                              only.clustering_markers = FALSE,
                              uniform.test            = "both",
                              th.pvalue               = 0.05,
                              th.IQR                  = 2,
                              density.PDFfile         = "qcUniformClusters_density.pdf",
                              density.PDFfile.dim     = c(17, 10),
                              heatmap.PDFfile         = "qcUniformClusters_heatmap.pdf",
                              heatmap.PDFfile.dim     = c(length(Results@marker.names), length(Results@cluster.names)/4),
							  tile.color              = "black",
                              verbose                 = TRUE,
                              ...){
    
    if (is.null(Results)) {
        stop("Error in qcUniformClusters: The 'Results' parameter can not be NULL")
    } else if (class(Results)[1] != "Results") {
        stop("Error in qcUniformClusters: The 'Results' parameter required a 'Results' object")
    }
    
    if (th.pvalue <= 0) {
        stop("Error in qcUniformClusters: The 'th.pvalue' parameter can not be negative")
    }
    
    if (th.IQR <= 0) {
        stop("Error in qcUniformClusters: The 'th.sd' parameter can not be negative")
    }
    
    flowset <- Results@flowset
    
    if (is.null(flowset)) {
        stop("Error in qcUniformClusters: The 'flowset' slot of the 'Results' object is not loaded, if the slot 'fcs.files' is not null, use the function load.flowSet before using the 'biplotViewer' function")
    }
    
    if (!is.null(heatmap.PDFfile) && class(heatmap.PDFfile)[1] != "character") {
        stop("Error in qcUniformClusters: The 'heatmap.PDFfile' parameter must be null or a character")
    }
    
    if (!is.null(density.PDFfile) && class(density.PDFfile)[1] != "character") {
        stop("Error in qcUniformClusters: The 'density.PDFfile' parameter must be null or a character")
    }

    if (is.null(uniform.test) || !is.element(uniform.test, c("unimodality", "spread", "both"))) {
        stop("Error in qcUniformClusters: The 'uniform.test' parameter must be a character among 'unimodality', 'spread' or 'both'")
    }
    
    if (!is.null(density.PDFfile)) {
        on.exit(dev.off())
        pdf(density.PDFfile, width = density.PDFfile.dim[1], height = density.PDFfile.dim[2])
    }

    if (is.null(clusters)) {
        clusters <- gtools::mixedsort(Results@cluster.names)
    } else if (all(clusters %in% Results@cluster.names)) {
        if (typeof(clusters) != "character") {
            stop("Error in qcUniformClusters: 'clusters' parameter must be a character vector")
        }
    } else {
        stop("Error in qcUniformClusters:\nUnknown clusters : ", paste(setdiff(unique(clusters), Results@cluster.names), collapse = " "))
    }
        
    message("[BEGIN] - generating Uniform Phenotypes QC")

    clustering.markers <- Results@clustering.markers
    
    if (only.clustering_markers) {
        markers <- clustering.markers
    } else {
        markers <- Results@marker.names
    }

    flowset <- Results@flowset
        
    accuracy.matrix <- matrix(nrow = length(clusters), ncol = length(markers), dimnames = list(clusters, markers))
    
    min <- floor(min(flowCore::fsApply(flowset[, flowCore::colnames(flowset) != "cluster"], flowCore::each_col, base::min, na.rm = TRUE), na.rm = TRUE))
    max <- ceiling(max(flowCore::fsApply(flowset[, flowCore::colnames(flowset) != "cluster"], flowCore::each_col, base::max, na.rm = TRUE), na.rm = TRUE))

    ordered.markers    <- c(gtools::mixedsort(clustering.markers),gtools::mixedsort(setdiff(markers, clustering.markers)))
    
    bold.markers       <- ifelse(is.element(ordered.markers, clustering.markers), "bold", "plain")
    colored.markers    <- ifelse(is.element(ordered.markers, clustering.markers), "blue", "black")

	count <- 0 
    for (cluster in clusters) {
        
        data.facet <- data.frame(ind = c(), subtitle = c(), upper = c(), lower = c(), median = c(), pinnacle = c())

        if (verbose) {
			count <- count + 1
            message(paste0("Cluster: ", count, " on ", length(clusters)))
        }
        expressions.list <- vector("list", length(flowset))
        for (sample in 1:length(flowset)) {
            frame <- flowset[[sample]]@exprs
            expressions.list[[sample]] <- frame[frame[, "cluster"] == cluster, colnames(frame) %in% ordered.markers]
        }
        expressions <- do.call(rbind, expressions.list)
        for (marker in ordered.markers) {
            if (nrow(expressions) > 1) {

                marker.expression <- expressions[, marker]
                subtitle <- ""
                if (uniform.test == "unimodality" || uniform.test == "both") {
                    # bug on some platforms: diptest message cannot be captured by utils::capture.output ...
					p.value <- diptest::dip.test(marker.expression, ...)$p.value
                    if (p.value < th.pvalue) {
                        uniform <- FALSE
                        subtitle <- paste0(subtitle," Non-unimodale distribution detected (p.value = ", round(p.value, 2), ")")
                    } else {
                        uniform <- TRUE
                        subtitle <- paste0(subtitle," Unimodale distribution detected (p.value = ", round(p.value, 2), ")")
                    }
                    facet <- data.frame(ind      = marker,
                                        subtitle = subtitle,
                                        upper    = NA,
                                        lower    = NA,
                                        median   = NA,
                                        pinnacle = NA)
                    
                }else {
					uniform <- TRUE
				} 
				
                if (uniform.test == "spread" || uniform.test == "both") {

                    quantile  <- quantile(marker.expression)
                    IQR       <- quantile[4] - quantile[2]
                    pinnacle  <- computemode(marker.expression)$y
                    
                    if (IQR < th.IQR) {
                        uniform <- ifelse(uniform, TRUE, FALSE)
                        subtitle <- paste0(subtitle,"\n"," Low spread detected (IQR < ", th.IQR, ")")
                    } else {
                        uniform <- FALSE
                        subtitle <- paste0(subtitle,"\n"," High spread detected (IQR > ", th.IQR, ")")
                    }

                    facet <- data.frame(ind      = marker,
                                        subtitle = subtitle,
                                        upper    = quantile[4],
                                        lower    = quantile[2],
                                        median   = quantile[3],
                                        pinnacle = pinnacle)
                }
                data.facet <- rbind(data.facet, facet)
                accuracy.matrix[cluster, marker] <- uniform
                
            } else {
                accuracy.matrix[cluster, marker] <- NA
            }
            
        }
        
        if (!is.null(density.PDFfile)) {

            data.expressions     <- as.data.frame(expressions)
            data.expressions     <- utils::stack(data.expressions)
            data.expressions$ind <- factor(data.expressions$ind, levels = ordered.markers)

            subtitle             <- paste(levels(data.expressions$ind), "\n", data.facet$subtitle)
            names(subtitle)      <- levels(data.expressions$ind)
            title                <- paste0("Densities of marker expressions for cluster : ", cluster, " (", format(nrow(expressions), big.mark = " "), " cells)")
            fill                 <- ifelse(accuracy.matrix[cluster, ], "green", "red")
            names(fill)          <- colnames(accuracy.matrix)

            data.facet$cm  <- data.facet$ind %in% clustering.markers
            data.facet$ind <- factor(data.facet$ind, levels = ordered.markers)
            
            plot <- ggplot2::ggplot() +
                    ggplot2::ggtitle(title) +
                    ggplot2::geom_density(data = data.expressions, ggplot2::aes_string(x = "values", fill = "ind"), alpha = 0.1, show.legend = FALSE) +
                    ggplot2::scale_fill_manual(name = "uniform", values = fill) +
                    ggplot2::scale_x_continuous(limits = c(min, max), expand = c(0.01, 0)) +
                    ggplot2::scale_y_continuous(expand = c(0.02, 0))
            
            if (uniform.test == "spread" || uniform.test == "both") {
                data.facet$IQR <- round(data.facet$upper - data.facet$lower, 2)

                plot <- plot + ggplot2::geom_vline(data = data.facet, ggplot2::aes_string(xintercept = "median"), color = "blue", size = 0.2) +
                               ggplot2::geom_vline(data = data.facet, ggplot2::aes_string(xintercept = "upper"), color = "blue", linetype = "dashed", size = 0.1) +
                               ggplot2::geom_vline(data = data.facet, ggplot2::aes_string(xintercept = "lower"), color = "blue", linetype = "dashed", size = 0.1) +
                               ggplot2::geom_segment(data = data.facet, ggplot2::aes_string(x = "median", y = "pinnacle*0.75", xend = "upper", yend = "pinnacle*0.75"), color = "blue", size = 0.1, arrow = ggplot2::arrow(length = ggplot2::unit(0.1, "cm"))) +
                               ggplot2::geom_segment(data = data.facet, ggplot2::aes_string(x = "median", y = "pinnacle*0.75", xend = "lower", yend = "pinnacle*0.75"), color = "blue", size = 0.1, arrow = ggplot2::arrow(length = ggplot2::unit(0.1, "cm"))) +
                               ggplot2::geom_label(data = data.facet, ggplot2::aes_string(x = "median", y = "pinnacle*0.75", label = "IQR"), color = "blue", size = 1, label.padding = ggplot2::unit(0.1, "lines"))
            }

            plot <- plot + ggplot2::geom_rect(data = data.facet, ggplot2::aes_string(color = "cm"), fill = NA, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, size = 1.5, show.legend = FALSE) +
                           ggplot2::scale_color_manual(values = c("grey80", "blue")) +
            ggplot2::facet_wrap(~ind, scales = "free", labeller = ggplot2::as_labeller(subtitle)) +
                           ggplot2::xlab("marker expressions") +
                           ggplot2::theme_bw() +
                           ggplot2::theme(strip.text.x = ggplot2::element_text(size = 5))
            plot(plot)

        }
    }
    if (!is.null(density.PDFfile)) {
        dev.off()
    }
    
    quality.matrix.cm       <- accuracy.matrix[, colnames(accuracy.matrix) %in% clustering.markers]
    uniform.cluster         <- apply(quality.matrix.cm, 1, all)
    perc                    <- length(uniform.cluster[uniform.cluster == TRUE])/length(uniform.cluster) * 100
    accuracy.matrix         <- data.frame(accuracy.matrix, uniform.cluster = uniform.cluster, check.names = FALSE)
    accuracy.matrix$cluster <- rownames(accuracy.matrix)
    
    if (!is.null(heatmap.PDFfile)) {
        pdf(heatmap.PDFfile, width = heatmap.PDFfile.dim[1], height = heatmap.PDFfile.dim[2])
        on.exit(dev.off())

        accuracy.matrix.melted           <- reshape2::melt(accuracy.matrix, id = "cluster")
        colnames(accuracy.matrix.melted) <- c("cluster", "marker", "uniform")
        
        ordered.markers    <- c(ordered.markers, "uniform.cluster")
        bold.markers       <- ifelse(is.element(ordered.markers, clustering.markers), "bold", "plain")
        colored.markers    <- ifelse(is.element(ordered.markers, clustering.markers), "blue", "black")
        
        accuracy.matrix.melted$cluster <- factor(accuracy.matrix.melted$cluster, levels = rev(clusters))
        accuracy.matrix.melted$marker  <- factor(accuracy.matrix.melted$marker, levels = ordered.markers, ordered = TRUE)
        
        details <- ""
        
        if (uniform.test == "unimodality" || uniform.test == "both") {
            details <- paste0(details, "having unimodal distribution (p.value < ", th.pvalue, ")")
        }
        if (uniform.test == "both") {
            details <- paste0(details, " and ")
        }
        if (uniform.test == "spread" || uniform.test == "both") {
            details <- paste0(details, "having low spread (th.IQR < ", th.IQR, ")")
        }
        
        title    <- paste0("QC: Uniform Phenotypes")
        subtitle <- paste0("percentage of clusters having only uniform clustering marker distributions = ", 
                           format(round(perc, 2), nsmall = 2), "% (marker ", details,")")
        
        plot <- ggplot2::ggplot(data = accuracy.matrix.melted) +
                ggplot2::ggtitle(bquote(atop(.(title), atop(italic(.(subtitle)), "")))) +
                ggplot2::geom_tile(ggplot2::aes_string(x = "marker", y = "cluster", fill = "uniform"), colour = tile.color) +
                ggplot2::scale_fill_manual(values = c("FALSE" = "red", "TRUE" = "green"), na.value = "grey50") +
                ggplot2::scale_x_discrete(expand = c(0, 0)) + 
                ggplot2::scale_y_discrete(expand = c(0, 0)) +
                ggplot2::xlab("markers") +
                ggplot2::ylab("clusters") +
                ggplot2::geom_vline(xintercept = (ncol(accuracy.matrix) - 1.5), colour = "black", size = 2) +
                ggplot2::geom_vline(xintercept = (length(clustering.markers) + 0.5), colour = "grey40", size = 1.5) +
                ggplot2::geom_vline(xintercept = (ncol(accuracy.matrix) - 0.5), colour = "black", size = 1) +
                ggplot2::geom_vline(xintercept = 0.5, colour = "black", size = 1) +
                ggplot2::geom_hline(yintercept = 0.5, colour = "black", size = 1) +
                ggplot2::geom_hline(yintercept = length(clusters) + 0.5, colour = "black", size = 1) +
                ggplot2::theme_bw() +
                ggplot2::theme(axis.text.x      = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, face = bold.markers, color = colored.markers),
                               legend.key       = ggplot2::element_blank(),
                               axis.line        = ggplot2::element_blank(),
                               panel.background = ggplot2::element_blank(),
                               panel.border     = ggplot2::element_blank(),
                               panel.grid.major = ggplot2::element_blank(),
                               panel.grid.minor = ggplot2::element_blank(),
                               plot.background  = ggplot2::element_blank())
        plot(plot)
    }
    message("[END] - generating Uniform Phenotypes QC")
    invisible(list(perc = perc, accuracy.matrix = accuracy.matrix))
    
}