#' @title Visualization of abundant clusters
#'
#' @description 
#' Generates a scatter plot representation showing for each cluster: its mean abundance and associated p-value.
#' 
#' This representation displays the p-value (shown as -log10(p-value)) in the X-axis and the mean of cells abundance in the Y-axis, in a two-dimensional chart. 
#' Each dot in the representation corresponds to a cell cluster, and both p-value and mean thresholds are shown using red dashed lines. 
#' Abundant Clusters are highlighted in red and labeled. 
#' The size of the dots is proportional to the total number of associated cells in the considered samples.
#' 
#' @details
#' By default, only significant abundant clusters are labeled. 
#' Labels for all clusters can be displayed by setting the 'show.all_labels' parameter to TRUE.
#' 
#' @param AC an object of class AC (object returned by the 'computeAC()' function)
#' @param show.cluster_sizes a logical specifying if dot sizes are proportional to number of associated cells
#' @param show.all_labels a logical specifying if all cluster labels must be shown. Only labels of significant clusters are displayed otherwise
#' @param show.on_device a logical specifying if the ggplot representation must be displayed on the device 
#' @param max.dots_size a numeric specifying the number of associated cells in the largest dot
#'
#' @return a 'ggplot' object
#' 
#' @export
#' 
#' @import ggplot2 ggrepel
abundantClustersViewer <- function(AC,
                                   show.cluster_sizes = TRUE,
                                   show.all_labels    = FALSE,
                                   show.on_device     = TRUE,
                                   max.dots_size      = max(AC@cluster.size)){

    if (!is.logical(show.cluster_sizes)) {
        stop("Error in abundantClustersViewer: The 'show.cluster_sizes' parameter must be a logical")
    }

    if (!is.logical(show.all_labels)) {
        stop("Error in abundantClustersViewer: The 'show.all_labels' parameter must be a logical")
    }

    if (!is.logical(show.on_device)) {
        stop("Error in abundantClustersViewer: The 'show.on_device' parameter must be a logical")
    }


    AC@results <- cbind(AC@results, cluster.size = AC@cluster.size)
    data.text <- AC@results

    if (!show.all_labels) {
        data.text <- subset(AC@results, AC@results$significant)
    }

    title      <- paste("Abundant Clusters")
    subtitle   <- ifelse(AC@use.percentages, "based on relative abundances", "based on absolutes abundances")
    if (AC@method == "t.test") {
        subtitle   <- paste0(subtitle, " using a Student t-test")
    } else if (AC@method == "wilcox.test") {
        subtitle   <- paste0(subtitle, " using a Mann-Whitney u-test")
    }
    
    if (AC@method.adjust != "none") {
        subtitle   <- paste0(subtitle, "(corrected using a ", AC@method.adjust, " method)")
    }
    
    AC@results <- AC@results[order(AC@results$cluster.size, decreasing = TRUE),]

    plot <-  ggplot2::ggplot(data = AC@results) +
             ggplot2::ggtitle(bquote(atop(.(title), atop(italic(.(subtitle)), "")))) +
             ggplot2::geom_hline(yintercept = AC@mu,
                                 linetype   = "dashed",
                                 alpha      = 0.3,
                                 color      = "red",
                                 size       = 1) +  
             ggplot2::geom_vline(xintercept = -log10(AC@th.pvalue),
                                 linetype   = "dashed",
                                 alpha      = 0.3,
                                 color      = "red",
                                 size       = 1)

    if (show.cluster_sizes) {

        plot <- plot + ggplot2::geom_point(ggplot2::aes_string(x = "-log10(pvalue)", y = "mean", fill = "significant", size = "cluster.size"),
                                           shape  = 21,
                                           colour = "grey40",
                                           stroke = 1) +
                       ggplot2::scale_size_area(name     = "number of cells",
                                                   limits   = c(0, max.dots_size),
                                                guide    = ggplot2::guide_legend(order = 2))
    }else{
        plot <- plot + ggplot2::geom_point(ggplot2::aes_string(x = "-log10(pvalue)", y = "mean", fill = "significant"),
                                           shape = 21,
                                           colour = "grey40",
                                           stroke = 1)
    }                  
    
    pvalues  <- c( -log10(AC@th.pvalue), -log10(AC@results$pvalue))
    x.max    <- ceiling(max(pvalues[!is.infinite(pvalues)]))
    x.breaks <- c(seq(0, x.max, by = 1), round(-log10(AC@th.pvalue), 2))
    
    y.max    <- ceiling(max(AC@mu, AC@results$mean))
    y.breaks <- seq(0, y.max, by = 1)

    plot <- plot +  ggrepel::geom_text_repel(data          = data.text, 
                                             ggplot2::aes_string(x = "-log10(pvalue)", y = "mean", label = "cluster"),
                                             size          = 3,
                                             box.padding   = grid::unit(0.35, "lines"),
                                             point.padding = grid::unit(0.3, "lines")) +
                    ggplot2::scale_fill_manual(values = c("grey", "red"), guide = ggplot2::guide_legend(order = 1)) +
                    ggplot2::scale_x_continuous(limits = c(0, x.max), minor_breaks = NULL, breaks = x.breaks) + 
                    ggplot2::scale_y_continuous(limits = c(0, y.max), minor_breaks = NULL, breaks = y.breaks) +
                    ggplot2::xlab("-log10(p-value)") +
                    ggplot2::ylab(ifelse(AC@use.percentages, "mean (% of cells)", "mean (# of cells)")) +                  
                    ggplot2::theme_bw() +
                    ggplot2::theme(legend.key = ggplot2::element_blank())
    
    if (show.on_device) {
        plot(plot)
    }

    invisible(plot)

}


#' @title Visualization of differentially abundant clusters
#'
#' @description
#' Generates a volcano plot representation showing for each cluster: its mean abundance fold-change and associated p-value. 
#' 
#' This representation displays the p-value (shown as -log10(p-value)) in the Y-axis and the fold-change of cell abundances, in the X-axis in a two-dimensional chart. 
#'Each dot in the representation corresponds to a cell cluster, and both p-value and fold-change thresholds are shown using red dashed lines.
#'Differentially Abundant Clusters are highlighted in red and labeled.
#'The size of dots is proportional to the total number of associated cells in the 2 conditions merged.
#'
#' @details 
#' By default, only significant differentially abundant clusters are labeled. 
#' Labels for all clusters can be displayed by setting the 'all.label' parameter to TRUE. 
#'
#' @param DAC an object of class 'DAC' (object returned by the 'computeDAC()' function)
#' @param fc.log2 a logical specifying if the fold-change or the log2(fold-change) must be used 
#' @param show.cluster_sizes a logical specifying if dot sizes are proportional to number of associated cells
#' @param show.all_labels a logical specifying if all cluster labels must be shown. Only labels of significant clusters are displayed otherwise
#' @param show.on_device a logical specifying if the ggplot representation must be displayed on the device 
#' @param max.dots_size a numeric specifying the number of associated cells in the largest dot
#'
#' @return a 'ggplot' object
#' 
#' @export
#'  
#' @import ggplot2
volcanoViewer <- function(DAC                = NULL,
                          fc.log2            = TRUE,
                          show.cluster_sizes = TRUE,
                          show.all_labels    = FALSE,
                          show.on_device     = TRUE,
                          max.dots_size      = max(DAC@cluster.size)) {

    if (!is.logical(fc.log2)) {
        stop("Error in volcanoViewer: The 'fc.log2' parameter must be a logical")
    }

    if (!is.logical(show.cluster_sizes)) {
        stop("Error in volcanoViewer: The 'show.cluster_sizes' parameter must be a logical")
    }

    if (!is.logical(show.all_labels)) {
        stop("Error in volcanoViewer: The 'show.all_labels' parameter must be a logical")
    }

    if (!is.logical(show.on_device)) {
        stop("Error in volcanoViewer: The 'show.on_device' parameter must be a logical")
    }
    
    th.fc <- DAC@th.fc
    
    if (fc.log2) {
        th.fc <- log2(th.fc)
        temp <- DAC@results$fold.change
        for (i in 1:nrow(DAC@results)) {
            DAC@results$fold.change[i] <- ifelse(temp[i] > 0, log2(DAC@results$fold.change[i]), -log2(abs(DAC@results$fold.change[i])))
        }
    }
    

    DAC@results <- cbind(DAC@results, cluster.size = DAC@cluster.size)

    DAC@results <- DAC@results[order(DAC@results$cluster.size, decreasing = TRUE), ]

    data.text <- DAC@results
    if (!show.all_labels) {
        data.text <- subset(DAC@results, DAC@results$significant)
    }

    if (all(is.na(DAC@results$fold.change))) {
        stop("Error, all cluster fold-changes are NA")
    }

    fc    <- DAC@results$fold.change
	fc    <- fc[!is.infinite(fc)]
	fc    <- fc[!is.na(fc)]
	x.min <- floor(min(fc))
	x.max <- ceiling(max(fc))
    
	x.max    <- max(x.max, abs(x.min), na.rm = TRUE)
    x.breaks <- c(round(c( -th.fc, th.fc), 2), seq( -x.max, x.max, by = 1))

    pvalues  <- c( -log10(DAC@th.pvalue), -log10(DAC@results$pvalue))
    y.max    <- ceiling(max(pvalues[!is.infinite(pvalues)]))
    y.breaks <- c(seq(0, y.max, by = 1), round(-log10(DAC@th.pvalue), 2))

    title.details <- ifelse(DAC@use.percentages, "based on % of cells", "based on # of cells")
    title.details <- paste()

    title <- paste("Differentially Abundant Clusters")
    subtitle <- ifelse(DAC@use.percentages, "based on relative abundances", "based on absolutes abundances")
    
    if (DAC@method == "t.test") {
        subtitle   <- paste0(subtitle, " using a Student ", ifelse(DAC@method.paired, "paired", "unpaired"), " t-test")
    } else if (DAC@method == "wilcox.test") {
        subtitle   <- paste0(subtitle, " using a Mann-Whitney", ifelse(DAC@method.paired, "paired", "unpaired"), " u-test")
    }
    
    if (DAC@method.adjust != "none") {
        subtitle   <- paste0(subtitle, "(corrected using a ", DAC@method.adjust, " method)")
    }

    plot <- ggplot2::ggplot(data = DAC@results, ggplot2::aes_string(x = "fold.change", y = "-log10(pvalue)")) +
            ggplot2::ggtitle(bquote(atop(.(title), atop(italic(.(subtitle)), "")))) +
            ggplot2::geom_vline(xintercept = c(th.fc, -th.fc),
                                linetype   = "dashed",
                                alpha      = 0.3,
                                color      = "red",
                                size       = 1) +  
            ggplot2::geom_hline(yintercept = -log10(DAC@th.pvalue),
                                linetype   = "dashed",
                                alpha      = 0.3,
                                color      = "red",
                                size       = 1)
    if (show.cluster_sizes) {

        plot <- plot + ggplot2::geom_point(ggplot2::aes_string(fill = "significant", size = "cluster.size"),
                                           shape  = 21,
                                           colour = "grey40",
                                           stroke = 1) +
                       ggplot2::scale_size_area(name     = "number of cells",
                                                limits   = c(0, max.dots_size),
                                                guide    = ggplot2::guide_legend(order = 2))
    } else {
        plot <- plot + ggplot2::geom_point(ggplot2::aes_string(fill = "significant"),
                                           shape = 21,
                                           colour = "grey40", 
                                           stroke = 1)
    }

    cond1 <- gsub("^mean.", "", colnames(DAC@results)[2])
    cond2 <- gsub("^mean.", "", colnames(DAC@results)[4])

    orientation <- paste0("\n", cond2, " <- enriched -> ", cond1)


    plot <- plot + ggrepel::geom_text_repel(data          = data.text,
                                            ggplot2::aes_string(label = "cluster"),
                                            size          = 3,
                                            box.padding   = grid::unit(0.35, "lines"),
                                            point.padding = grid::unit(0.3, "lines")) +
                   ggplot2::scale_fill_manual(values = c("grey", "red")) +
                   ggplot2::scale_x_continuous(limits = c(-x.max, x.max), minor_breaks = NULL, breaks = x.breaks) +
                   ggplot2::scale_y_continuous(limits = c(0, y.max), minor_breaks = NULL, breaks = y.breaks) +
                   ggplot2::xlab(paste0(ifelse(fc.log2, "log2(fold.change)", "fold.change"), orientation)) +
                   ggplot2::ylab("-log10(p-value)") +
                   ggplot2::theme_bw() +
                   ggplot2::theme(legend.key = ggplot2::element_blank())

    if (show.on_device) {
        plot(plot)
    }

    invisible(plot)
}


#' @title Visualization of correlated clusters
#'
#' @description
#' Generates a volcano plot representation showing for each cluster: its correlation coefficient and associated p-value.
#'
#'This representation displays the p-value (shown as -log10(p-value)) in the Y-axis and the correlation coefficient in the X-axis, in a two-dimensional chart.
#'Each dot in the representation corresponds to a cell cluster, and both correlation coefficient (positive and negative) and p-value thresholds are shown using red dashed lines. 
#'Correlated Clusters are highlighted in red and labeled. 
#'The size of dots is proportional to the total number of associated cells in the considered samples.  
#' 
#' @details 
#' By default, only significant correlated clusters are labeled.
#' Labels for all clusters can be displayed by setting the 'all.label' parameter to TRUE. 
#'
#' @param CC an object of class 'CC' (object returned by the 'computeCC()' function)
#' @param show.cluster_sizes a logical specifying if dot sizes are proportional to number of associated cells
#' @param show.all_labels a logical specifying if all cluster labels must be shown. Only labels of significant clusters are displayed otherwise
#' @param show.on_device a logical specifying if the ggplot representation must be displayed on the device 
#' @param max.dots_size a numeric specifying the number of associated cells in the largest dot
#' @param varname a character indicating the name of the variable to display in the title
#'
#' @return a 'ggplot' object
#' 
#' @export
#' 
#' @import ggplot2 ggrepel
correlatedClustersViewer <- function(CC,
                                     show.cluster_sizes = TRUE,
                                     show.all_labels    = FALSE,
                                     show.on_device     = TRUE,
                                     max.dots_size      = max(CC@cluster.size),
									 varname            = NULL) {

    if (!is.logical(show.cluster_sizes)) {
        stop("Error in correlatedClustersViewer: The 'show.cluster_sizes' parameter must be a logical")
    }

    if (!is.logical(show.all_labels)) {
        stop("Error in correlatedClustersViewer: The 'show.all_labels' parameter must be a logical")
    }

    if (!is.logical(show.on_device)) {
        stop("Error in correlatedClustersViewer: The 'show.on_device' parameter must be a logical")
    }

    CC@results <- cbind(CC@results, cluster.size = CC@cluster.size)

    CC@results <- CC@results[order(CC@results$cluster.size, decreasing = TRUE),]

    data.text <- CC@results
    if (!show.all_labels) {
        data.text <- subset(CC@results, CC@results$significant)
    }
    
    title    <- paste("Correlated Clusters")
    subtitle <- ifelse(CC@use.percentages, "based on relative abundances", "based on absolutes abundances")
    method <- CC@method
    substring(method, 1) <- toupper(substring(method, 1, 1))
       subtitle   <- paste0(subtitle, " using a ", method, " correlation method")

    if (CC@method.adjust != "none") {
        subtitle   <- paste0(subtitle, "(corrected using a ", CC@method.adjust, " method)")
    }
		
    plot <- ggplot2::ggplot(data = CC@results) +
            ggplot2::geom_hline(yintercept = -log10(CC@th.pvalue),
                                linetype   = "dashed",
                                alpha      = 0.3,
                                color      = "red",
                                size       = 1) +  
            ggplot2::geom_vline(xintercept = c(-CC@th.correlation, CC@th.correlation),
                                linetype   = "dashed",
                                alpha      = 0.3,
                                color      = "red",
                                size       = 1)
	
	if(!is.null(varname)){
		subtitle_var <- paste0("with variable: ",varname)
		plot <- plot+ggplot2::ggtitle(bquote(atop(.(title), atop(italic(.(subtitle)),italic(.(subtitle_var)), ""))))
	}else{
		plot <- plot+ggplot2::ggtitle(bquote(atop(.(title), atop(italic(.(subtitle)), ""))))
	}

    if (show.cluster_sizes) {

        plot <- plot + ggplot2::geom_point(ggplot2::aes_string(x = "correlation", y = "-log10(pvalue)", fill = "significant", size = "cluster.size"),
                                           shape = 21,
                                           colour = "grey40",
                                           stroke = 1) +
                        ggplot2::scale_size_area(name     = "number of cells",
                                                 limits   = c(0, max.dots_size),
                                                 guide    = ggplot2::guide_legend(order = 2))
    }else{
        plot <- plot + ggplot2::geom_point(ggplot2::aes_string(x = "correlation", y = "-log10(pvalue)", fill = "significant"),
                                           shape = 21,
                                           colour = "grey40",
                                           stroke = 1)
    }

    x.breaks <- c( -CC@th.correlation, CC@th.correlation, seq(-1, 1, by = 0.1))
    pvalues  <- c(-log10(CC@th.pvalue), -log10(CC@results$pvalue))
    y.max    <- ceiling(max(pvalues[!is.infinite(pvalues)]))
    
    y.breaks <- c(seq(0, y.max, by = 1), round(-log10(CC@th.pvalue), 2))
    
    plot <- plot +  ggrepel::geom_text_repel(data          = data.text, 
                                             ggplot2::aes_string(x = "correlation", y = "-log10(pvalue)", label = "cluster"),
                                             size          = 3,
                                             box.padding   = grid::unit(0.35, "lines"),
                                             point.padding = grid::unit(0.3, "lines")) +
                    ggplot2::scale_fill_manual(values = c("grey", "red")) +
                    ggplot2::scale_x_continuous(minor_breaks = NULL, limits = c(-1, 1), breaks = x.breaks) +
                    ggplot2::scale_y_continuous(limits = c(0, y.max), minor_breaks = NULL, breaks = y.breaks) +
                    ggplot2::xlab(paste(CC@method, "coefficient of correlation")) +
                    ggplot2::ylab("-log10(p-value)") +
                    ggplot2::theme_bw() +
                    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5), legend.key = ggplot2::element_blank())

    if (show.on_device) {
        plot(plot)
    }

    invisible(plot)

}


#' @title Visualization of abundance profiles classification
#'
#' @description
#' Generates a colored circle packing representation showing for each cluster: its associated class based on its abundance profile.  
#' 
#' @details 
#' The sizes of the dots are proportionals to the total number of associated cells to each cluster.
#' Circle packing classes are sorted by the number of clusters in each class. 
#'
#' @param AP an object of class 'AP' (object returned by the 'classifyAbundanceProfiles()' function)
#' @param show.cluster_sizes a logical specifying if dot sizes are proportional to number of associated cells
#' @param show.on_device a logical specifying if the representation will be displayed on device 
#'
#' @return a 'ggplot' object
#' 
#' @export
#'
#' @import ggplot2 grDevices gridExtra
circlesPackingViewer <- function(AP,
                                 show.cluster_sizes = TRUE,
                                 show.on_device     = TRUE) {

    if (!is.logical(show.on_device)) {
        stop("Error in circlesPackingViewer: The 'show.on_device' parameter must be a logical")
    }

    classes <- AP@classes
    classes <- stats::na.omit(classes)
    sorted.classes <- names(sort(table(classes$class), decreasing = TRUE))

    colours <- grDevices::rainbow(n = length(sorted.classes))

    plots <- list()
    
    for (i in sorted.classes) {

        same.class <- classes[classes$class == i,]
        
        if(show.cluster_sizes){
            cluster.sizes <- AP@cluster.size[as.character(same.class$cluster)]
        } else {
            cluster.sizes <- rep(1000, nrow(same.class))
        }
        
        same.class <- data.frame(cluster = same.class$cluster, size = cluster.sizes, stringsAsFactors = FALSE)
        plots[[i]] <- buildCircles(circles = same.class, colours[as.numeric(i)], class = i)

    }

    title    <- paste("Circles Packing Viewer")
    subtitle <- paste0("based on abundance profiles using ", AP@method, " method")

    grob.title <- grid::textGrob(bquote(atop(.(title), atop(italic(.(subtitle)), ""))))

    if (length(plots) == 0) {
        grobs <- gridExtra::arrangeGrob(grobs = list(grid::textGrob("Empty AP object")),
                                        top   = grob.title)
    } else {
        if(show.cluster_sizes) {
            extra_space <- (length(plots) %% 3)
    
            if (extra_space) {
                empty_space <- 3 - extra_space
                for (i in 1:empty_space) {
                    plots[[length(plots) + 1]] <- grid::rectGrob(gp = grid::gpar(col = 0))
                }
            }
    
            plots[[length(plots) + 1]] <- grid::rectGrob(gp = grid::gpar(col = 0))
            plots[[length(plots) + 1]] <- buildCirclesLegend()
            plots[[length(plots) + 1]] <- grid::rectGrob(gp = grid::gpar(col = 0))
            
        }
        
        grobs <- gridExtra::arrangeGrob(grobs = plots, ncol = 3,
                                        top   = grob.title)

    }
        
    if (show.on_device) {
        grid::grid.newpage()
        grid::grid.draw(grobs)
    }

    invisible(grobs)

}


#' @title Visualization of cluster sizes
#'
#' @description 
#' The Count Viewer aims to visualize the number of cells associated to each cluster. 
#' This representation displays the clusters in the X-axis and the total number of associated cells in the Y-axis.
#' Additionally, the numbers of cells associated to each cluster for each sample are also displayed using a jitter representation.
#' The size of the dots is proportional to the total number of associated cells. 
#' 
#' @details 
#' By default, all clusters will be displayed but the representation can be restricted to a set of selected samples (using the `samples` parameter) or to a set of selected clusters (using the `clusters` parameter).
#'
#' @param Results a 'Results' object
#' @param samples a character vector providing the sample names to used (all samples by default)
#' @param clusters a character vector containing the clusters names to be visualized (by default all clusters will be displayed)
#' @param min.cells a numeric specifying the minimum number of cell (sum of all selected samples) to display a cluster
#' @param sort a logical specifying if the clusters will be to be sorted (descending) based on the sum of all selected samples for each cluster
#' @param max.dots_size a numeric specifying the number of cells in the largest dot
#' @param show.samples a logical specifying if the number of cells for all selected samples will be displayed
#' @param show.on_device a logical specifying if the representation will be displayed on device 
#'
#' @return a 'ggplot' object
#' 
#' @export
#' 
#' @import ggplot2 reshape2
countViewer <- function(Results,
                        samples        = NULL,
                        clusters       = NULL,
                        min.cells      = 0,
                        sort           = TRUE,
                        max.dots_size  = NULL,
                        show.samples   = TRUE,
                        show.on_device = TRUE) {

    if (is.null(Results)) {
        stop("Error in countViewer: 'Results' parameter can not be NULL")
    } else if (class(Results)[1] != "Results") {
        stop("Error in countViewer: 'Results' parameter must be a 'Results' object")
    }
    
    data <- Results@cluster.abundances

    if (is.null(samples)) {
        data <- data
    } else if (!all(samples %in% Results@sample.names)) {
        stop("Error in countViewer: 'samples' parameter must contain only samples names\n Unknown sample names: ",
             paste(setdiff(unique(samples), Results@sample.names), collapse = " "))
    } else {
        data <- data[, samples, drop = FALSE]
    }

    if (is.null(clusters)) {
        clusters <- Results@cluster.names
    } else if (all(clusters %in% Results@cluster.names)) {
        if (typeof(clusters) != "character") {
            stop("Error in countViewer: 'clusters' parameter must be a character vector")
        }
        clusters <- unique(clusters)
        data <- data[clusters, , drop = FALSE]
    } else {
        stop("Error in countViewer:\nUnknown clusters : ", paste(setdiff(unique(clusters), Results@cluster.names), collapse = " "))
    }

    if (min.cells < 0) {
        stop("Error in countViewer: 'min.cells' parameter must be a positive integer")
    }
    if (!is.logical(sort)) {
        stop("Error in countViewer: 'sort' parameter must be a logical")
    }
    if (!is.logical(show.samples)) {
        stop("Error in countViewer: 'show.samples' parameter must be a logical")
    }
    if (!is.logical(show.on_device)) {
        stop("Error in countViewer: 'show.on_device' parameter must be a logical")
    }

    cells.number <- sum(colSums(data))

    data <- cbind(data, "sum.of.samples" = apply(data, 1, sum))
    data <- data[data$sum.of.samples > min.cells,]

    data <- cbind("clusters" = rownames(data), data)

    if (sort) {
		order <- order(data$sum.of.samples,decreasing=TRUE)
        data$clusters <- factor(data$clusters,levels=data$clusters[order])
    } else {
        data$clusters <- as.factor(data$clusters)
    }

    
    
    data.melted <- reshape2::melt(data, id = c("clusters"))
    colnames(data.melted) <- c("clusters", "samples", "values")
    data.melted$total <- ifelse(data.melted[, "samples"] == "sum.of.samples", "sum for selected samples", "")

    if (is.null(max.dots_size)) {
        max.dots_size <- max(data.melted$values, na.rm = TRUE)
    }
    
    plot <- ggplot2::ggplot(data = data.melted) +
            ggplot2::ggtitle(paste("Count Viewer (", format(cells.number, big.mark = " "), " cells)", sep = ""))
    if (show.samples) {
        plot <- plot + ggplot2::geom_jitter(data = subset(data.melted, samples != "sum.of.samples"),
                                            ggplot2::aes_string(x = "clusters", y = "values", size = "values", fill = "samples"),
                                            height = 0,
                                            width  = 0.5,
                                            shape  = 21,
                                            alpha  = 0.4)
    }

    plot <- plot + ggplot2::geom_point(data = subset(data.melted, samples == "sum.of.samples"),
                                       ggplot2::aes_string(x = "clusters", y = "values", size = "values", shape = "total"),
                                       fill = "grey40") +
                   ggplot2::scale_size_area(name = "# of cells", limits = c(0, max.dots_size)) +
                   ggplot2::scale_shape_manual(values = 21) +
                   ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, 1.1 * max(data.melted$values))) +
                   ggplot2::xlab("clusters") +
                   ggplot2::ylab("# of cells") +
                   ggplot2::theme_bw() +
                   ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5), legend.key = ggplot2::element_blank())

    if (show.on_device) {
        plot(plot)
    }

    invisible(plot)
}


#' @title Visualization of combined SPADE trees
#' 
#' @description 
#' The Tree Viewer aims to visualize the SPADE tree representations. 
#' This representation displays the SPADE cell clusters using this minimal spanning tree layout computed by SPADE. 
#' In such tree, each node represents a cell cluster and nodes are linked based on theirs phenotype similarities.
#' This viewer improves the original SPADE tree representations by allowing to combine SPADE trees from multiple samples. 
#' 
#' Significant clusters can be highlighted (node borders are then colored in blue) by providing a `AC`, `DAC`, or `CC` object (using the `highlight` parameter).
#' As with the original SPADE tree representations, nodes can be colored based on the marker median expression of a specific marker (using the `marker` parameter). 
#' 
#' @details 
#' The size of tree nodes is proportional to the number of cells in each cluster. 
#' If the 'stat.object' parameter is provided node outlines are colored according to clusters significance.
#' 
#' @param Results a 'Results' object (with 'graph' and 'graph.layout' slots not null)
#' @param samples a character vector providing the sample names to used (all samples by default)
#' @param marker a character specifying the marker to be displayed 
#' @param highlight an AC, DAC or CC object to highlight identified significant clusters in the SPADE tree
#' @param show.on_device a logical specifying if the representation will be displayed on device 
#'
#' @return a 'ggplot' object
#'
#' @export
#'
#' @import reshape2 ggplot2 grid igraph
treeViewer <- function(Results,
                       samples        = NULL,
                       highlight      = NULL,
                       marker         = NULL,
                       show.on_device = TRUE) {

    if (is.null(Results)) {
        stop("Error in treeViewer: 'Results' parameter can not be NULL")
    } else if (class(Results)[1] != "Results") {
        stop("Error in treeViewer: 'Results' parameter required a 'Results' object")
    } else if (is.null(Results@graph)) {
        stop("Error in treeViewer: 'Results' parameter required a 'Results' object with a not null 'graph' slot")
    } else if (is.null(Results@graph.layout)) {
        stop("Error in treeViewer: 'Results' parameter required a 'Results' object with a not null 'graph.layout' slot")
    }  
    
    if(length(Results@marker.names) == 0){
        stop("Error in treeViewer: 'Results' object must contain phenotypes")
    }
    
    data   <- Results@cluster.abundances
    
    if (is.null(samples)) {
        data <- data
    } else if (!all(samples %in% Results@sample.names)) {
        stop("Error in treeViewer: 'samples' parameter must contains only samples names\n Unknown sample names: ",
             paste(setdiff(unique(samples), Results@sample.names), collapse = " "))
    } else {
        data <- data[, samples, drop = FALSE]
    }
    highlight.name <- class(highlight)[1]
    if (!is.null(highlight) && !is.element(highlight.name, c("AC", "DAC", "CC"))) {
        stop("Error in treeViewer: 'highlight' parameter must be NULL or either a AC, DAC or CC object")
    }

    if (!is.null(marker) && !is.element(marker, Results@marker.names)) {
        stop(paste0("Error in treeViewer: 'marker' parameter: ", marker, ", is unknown"))
    }

    if (!is.logical(show.on_device)) { stop("Error in treeViewer: 'show.on_device' parameter must be a logical") }

    vertex.size <- apply(data, 1, sum)
    
    pos.vertex  <- data.frame(id   = as.character(1:nrow(Results@graph.layout)),
                              x    = Results@graph.layout[, 1],
                              y    = Results@graph.layout[, 2],
                              size = vertex.size)
    
    if (!is.null(marker)) {
        if (!is.null(samples)) {
            expr   <- subset(Results@cluster.phenotypes, sample %in% samples, drop = FALSE)
        }else{
            expr   <- Results@cluster.phenotypes
        }

        expr           <- expr[, c("cluster", "sample", marker)]
        expr           <- reshape2::dcast(expr, cluster ~ sample, value.var = marker)
        
        rownames(expr) <- expr$cluster
        expr           <- expr[, -1]
        mean.expr      <- apply(expr, 1, mean, na.rm = TRUE)

        pos.vertex     <- cbind(pos.vertex, marker = mean.expr)
        colnames(pos.vertex)[ncol(pos.vertex)] <- marker

        max.mean.expr <- ceiling(max(mean.expr, na.rm = TRUE))
        seq.mean.expr <- seq(from = -1, to = max.mean.expr, by = 1)

    }   
    
    edges    <- igraph::get.edgelist(Results@graph,names = FALSE)
    pos.edge <- data.frame(x    = Results@graph.layout[edges[, 1], 1],
                           xend = Results@graph.layout[edges[, 2], 1],
                           y    = Results@graph.layout[edges[, 1], 2],
                           yend = Results@graph.layout[edges[, 2], 2])
    
    cells.number <- sum(colSums(data))
    
    plot <- ggplot2::ggplot() +
            ggplot2::ggtitle(paste("Tree Viewer (", format( cells.number, big.mark = " "), " cells)", sep = "")) +
            ggplot2::geom_segment(data = pos.edge, ggplot2::aes_string(x = "x", xend = "xend", y = "y", yend = "yend"))
    
    if (!is.null(highlight)) {
        pos.vertex[, highlight.name] <- highlight@results$significant
        if (!is.null(marker)) {
            plot <- plot + ggplot2::geom_point(data   = pos.vertex,
                                               ggplot2::aes_string(x = "x", y = "y", size = "size", fill = marker, colour = highlight.name),
                                               stroke = 2.5,
                                               shape = 21) +
                           ggplot2::scale_fill_gradient(low = "#ECE822", high = "#EE302D", limits = c(-1, max.mean.expr), breaks = seq.mean.expr)
        }else{
            plot <- plot + ggplot2::geom_point(data   = pos.vertex,
                                               ggplot2::aes_string(x = "x", y = "y", size = "size", colour = highlight.name),
                                               fill   = "grey80",
                                               stroke = 2.5,
                                               shape  = 21)
        }
    }else{
        if (!is.null(marker)) {
            plot <- plot + ggplot2::geom_point(data   = pos.vertex,
                                               ggplot2::aes_string(x = "x", y = "y", size = "size", fill = marker),
                                               stroke = 2.5,
                                               shape  = 21) +
                           ggplot2::scale_fill_gradient(low = "#ECE822", high = "#EE302D", limits = c(-1, max.mean.expr), breaks = seq.mean.expr) 
        } else {
            plot <- plot + ggplot2::geom_point(data   = pos.vertex,
                                               ggplot2::aes_string(x = "x", y = "y", size = "size"),
                                               fill   = "grey80",
                                               stroke = 2.5,
                                               shape  = 21)
        }
    }
    plot <- plot + ggplot2::scale_color_manual(values = c("black", "blue")) +
                   ggplot2::scale_size_area(name = "# of cells", max_size = 15) +
                   ggrepel::geom_label_repel(data          = pos.vertex, 
                                             ggplot2::aes_string(x = "x", y = "y", label = "id"),
                                             size          = 4,
                                             color         = "black",
                                             box.padding   = grid::unit(0.1, "lines"),
                                             point.padding = grid::unit(0.1, "lines")) +
                   ggplot2::coord_fixed() +
                   ggplot2::theme(panel.background = ggplot2::element_blank(),
                                  panel.border     = ggplot2::element_blank(),
                                  axis.text.x      = ggplot2::element_blank(),
                                  axis.text.y      = ggplot2::element_blank(),
                                  axis.title.x     = ggplot2::element_blank(),
                                  axis.title.y     = ggplot2::element_blank(),
                                  panel.grid.minor = ggplot2::element_blank(),
                                  panel.grid.major = ggplot2::element_blank(),
                                  axis.ticks       = ggplot2::element_blank(),
                                  legend.position  = "right",
                                  legend.key       = ggplot2::element_blank())

    if (show.on_device) {
        plot(plot)
    }

    invisible(plot)
}


#' @title Visualization of all clusters phenotypes using categorical heatmap
#' 
#' @description
#' The Heatmap Viewer aims to visualize the phenotypes of all cell clusters.
#' This representation displays marker expressions of all clusters using a categorical heatmap (5 categories are computed by default). 
#' The range expression of each cell marker is discretized in several categories between bounds of marker expressions.
#' Each marker of each cluster is then categorized into one category based on the mean of median marker expressions.
#' Markers used as clustering markers are shown in blue.
#' 
#' @details
#' Both cell clusters and cell markers can be clustered using a hierarchical clustering (using the `dendrograms`, generating marker and cluster dendrograms by default).
#'
#' It is to note that the Heatmap Viewer is the default viewer for `Results` objects.
#'
#' The structures of the marker and cluster hierarchical clustering dendrograms are also returned by this fucntion.
#'
#' @param Results a 'Results' object
#' @param clusters a character vector providing the clusters to be used (all clusters by default)
#' @param markers a character vector providing the markers to be used (all markers by default)
#' @param num a numeric value specifying the number of markers expression categories to be used
#' @param dendrograms a character specifying which dendrograms must be build ("markers", "clusters", "both" or "none")
#' @param tile.color a character specifying the border color of the tiles (NA to remove tile borders)
#' @param show.on_device a logical specifying if the representation will be displayed on device
#'
#' @return a 'ggplot' object. This object also contains the categorical phenotipical heatmap (pheno.table attribute) and the marker and cluster hierarchical clustering dendrograms (markers.hc and clusters.hc attributes)
#'
#' @export
#'
#' @import grid reshape2
#' @importFrom gtools mixedsort
heatmapViewer <- function(Results,
                          clusters       = NULL,
                          markers        = NULL,
                          num            = 5,
                          dendrograms    = "both",
						  tile.color     = "black",
                          show.on_device = TRUE) {
    
    if (is.null(Results)) {
        stop("Error in heatmapViewer: 'Results' parameter can not be NULL")
    } else if (class(Results)[1] != "Results") {
        stop("Error in heatmapViewer: 'Results' parameter must be a 'Results' object")
    }
    
    if(length(Results@marker.names) == 0){
        stop("Error in heatmapViewer: 'Results' object must contain phenotypes")
    }
    
    data <- Results@cluster.phenotypes
	
    if (is.null(clusters)) {
        clusters <- gtools::mixedsort(Results@cluster.names)
    } else if (all(clusters %in% Results@cluster.names)) {
        if (typeof(clusters) != "character") {
            stop("Error in heatmapViewer: 'clusters' parameter must be a character vector")
        }
        clusters.select <- data[, "cluster"] %in% clusters
        data            <- data[clusters.select,]
    } else {
        stop("Error in heatmapViewer:\nUnknown clusters : ", paste(setdiff(unique(clusters), Results@cluster.names), collapse = " "))
    }
    
	if (is.null(markers)) {
        markers <- gtools::mixedsort(Results@marker.names)
    } else if (all(markers %in% Results@marker.names)) {
        data <- data[, c("sample", "cluster", markers)]
    } else {
        stop("Error in heatmapViewer: Unknown markers :", paste(setdiff(unique(markers), Results@marker.names), collapse = " "))
    }
    
    if (num <= 0) {
        stop("Error in heatmapViewer: 'num' parameter must be a strictly positive integer")
    }
    
    if (is.null(dendrograms) || !is.element(dendrograms, c("markers", "clusters", "both", "none"))) {
        stop("Error in heatmapViewer: 'dendrograms' parameter must be a character among 'markers', 'clusters', 'both' or 'none'")
    }
    
    if (!is.logical(show.on_device)) {
        stop("Error in heatmapViewer: 'show.on_device' parameter must be a logical")
    }
    
	pheno.table <- computePhenoTable(data, bounds = Results@bounds, num = num)
    
	if(all(is.na((pheno.table)))){
		stop("Error in heatmapViewer: filtered phenotypical heatmap is full of small clusters")
	}
	
    unrepresented <- setdiff(clusters, unique(pheno.table$cluster))
    
    for (i in seq_len(length(unrepresented))) {
        na.represented <- data.frame(cluster = rep(unrepresented[i], length(markers)), marker = markers, value = rep(-1, length(markers)))
        pheno.table    <- rbind(pheno.table, na.represented)
    }

    pheno.table <- reshape2::dcast(pheno.table, cluster ~ marker)
    
    cluster.colnames <- pheno.table$cluster

    pheno.table <- pheno.table[, 2:ncol(pheno.table)]
    
    pheno.table <- t(pheno.table)

    pheno.table           <- as.matrix(pheno.table)
    colnames(pheno.table) <- cluster.colnames

    pheno.table <- pheno.table[, clusters]
    pheno.table <- pheno.table[markers, ]
	
    dendrograms <- ifelse(dendrograms == "markers", "row", dendrograms)
    dendrograms <- ifelse(dendrograms == "clusters", "col", dendrograms)
    
    plot.elements <- ggheatmap(pheno.table, num = num, clustering.markers = Results@clustering.markers, dendrograms = dendrograms)
    row.hc        <- plot.elements$row.hc
	col.hc        <- plot.elements$col.hc
	plot.elements$row.hc <- NULL
	plot.elements$col.hc <- NULL
	
    title    <- "Heatmap Viewer"
    bounds   <- as.numeric(row.names(Results@bounds))
    subtitle <- paste0("catogories are computed between ", (bounds[1] * 100), "% and ", (bounds[2] * 100), "% percentiles of expression.")

    top  <- grid::textGrob(bquote(atop(.(title), atop(italic(.(subtitle)), ""))))
    plot <- ggheatmap.combine(plot.elements, top = top)

    if (show.on_device) {
        grid::grid.newpage()
        grid::grid.draw(plot)
    }else{
		grDevices::dev.off()
	}

	plot$markers.hc  <- row.hc
	plot$clusters.hc <- col.hc
	plot$pheno.table <- pheno.table
	
    invisible(plot)

}


#' @title Visualization of cluster enrichment profiles conditions
#'
#' @description
#' The Boxplot Viewer aims to visualize and compare the abundances between several biological conditions for one single cluster or for a set of combined clusters.
#'
#' This representation displays cell cluster abundances of each biological condition using boxplots.
#' Additionally, the abundance density of each biological condition can be visualize using violin representation.
#' Cluster abundance of each sample can also be visualized via colored dots.
#'
#' Abundance from several clusters can be a by providing one or several cluster names to the `clusters` parameter.
#' Biological conditions information must be assigned to the samples in order to use this viewer (provided by the slot `assignments` in the `Results` object).
#' The cell cluster abundances could be displayed as percentages or absolute numbers using the `use.percentages` parameter (TRUE by default).
#' 
#' @details
#' Biological conditions are specified in the Results object at the importation or using the 'assignContext()' function.
#' Cells clusters are colored based on theirs associated biological samples.
#' 
#' @param Results a 'Results' object
#' @param samples a character vector providing the sample names to used (all samples by default)
#' @param clusters a character vector containing the clusters names to be visualized (required)
#' @param use.percentages a logical specifying if the visualization must be performed on percentage
#' @param show.legend a logical specifying if the legend must be displayed
#' @param show.violin a logical specifying if the count distribution must be displayed
#' @param show.on_device a logical specifying if the representation will be displayed on device 
#'
#' @return a 'ggplot' object
#'
#' @export
#' 
#' @import grid reshape2 ggplot2
boxplotViewer <- function(Results,
                          samples         = NULL,
                          clusters        = NULL,
                          use.percentages = TRUE,
                          show.legend     = TRUE,
                          show.violin     = TRUE,
                          show.on_device  = TRUE) {

    if (is.null(Results)) {
        stop("Error in boxplotViewer: 'Results' parameter can not be NULL")
    } else if (class(Results)[1] != "Results") {
        stop("Error in boxplotViewer: 'Results' parameter must be a 'Results' object")
    }
    
    if (is.null(samples)) {
        samples     <- Results@sample.names
        cluster.abundances <- Results@cluster.abundances
    } else if (!all(samples %in% Results@sample.names)) {
        stop("Error in boxplotViewer: 'samples' parameter must contains only samples names\n Unknown sample names: ",
             paste(setdiff(unique(samples), Results@sample.names), collapse = " "))
    } 
    
    data        <- Results@cluster.abundances[, samples, drop = FALSE]
    cluster.abundances <- data

    if (!is.logical(use.percentages)) {
        stop("Error in boxplotViewer: 'use.percentages' parameter must be a logical")
    }

    if (use.percentages) {
        data.percent <- prop.table(as.matrix(data), 2) * 100
        data         <- data.frame(data.percent, check.names = FALSE)
        legendy      <- "% of cells relative to parent"
    } else {
        legendy      <- "# of cells"
    }

    if (is.null(clusters)) {
        stop("Error in boxplotViewer: 'clusters' parameter is required")
    } else if (all(clusters %in% Results@cluster.names)) {
        if (typeof(clusters) != "character") {
            stop("Error in boxplotViewer: 'clusters' parameter must be a character vector")
        }
        clusters    <- unique(clusters)
        cluster.abundances <- cluster.abundances[clusters,]
        data        <- data[clusters,]
        data        <- apply(data, 2, sum)
    } else {
        stop("Error in boxplotViewer:\nUnknown clusters : ", paste(setdiff(unique(clusters), Results@cluster.names), collapse = " "))
    }

    if (!is.logical(show.on_device)) {
        stop("Error in boxplotViewer: 'show.on_device' parameter must be a logical")
    }

    names.palette  <- unique(colnames(Results@cluster.abundances))
    palette        <- ggcolors(length(names.palette))
    names(palette) <- names.palette

    assignments <- Results@assignments
    if (!is.null(assignments)) {
        order <- unique(assignments$bc)
        assignments <- assignments[samples, , drop = FALSE]
        order       <- intersect(order, unique(assignments$bc))
    } else if (is.element("bc", colnames(assignments))) {
        stop("Error in boxplotViewer: 'assignments' slot must contain the column 'bc' in the provided 'Results' object")
    } else {
        stop("Error in boxplotViewer: 'assignments' slot in the provided Results object is required")
    }

    data      <- data.frame(samples = names(data), value = data, cond = assignments[names(data), "bc"])
    data$cond <- factor(data$cond, levels = order)

    max.value <- max(data$value, na.rm = TRUE)
    max.value <- max.value + 0.1 * max.value + 1

    cells.number <- sum(colSums(cluster.abundances))

    plot <- ggplot2::ggplot(data = data, ggplot2::aes_string(x = "cond", y = "value")) +
            ggplot2::ggtitle(paste("Boxplot Viewer - cluster: ", paste0(clusters, collapse = ", "), " (", format(cells.number, big.mark = " "), " cells) ", sep = "")) +
            ggplot2::geom_boxplot() +
            ggplot2::geom_jitter(ggplot2::aes_string(color = "samples"), width = 0.2, show.legend = show.legend) +
            ggplot2::scale_color_manual(values = palette)
            
    if (show.violin) {
        plot <- plot + ggplot2::geom_violin(alpha = 0.05, fill = "red", colour = "red")
    }

    if (use.percentages) {
        plot <- plot + ggplot2::scale_y_continuous(limits = c(0, max.value), breaks = round(seq(0, max.value)), minor_breaks = NULL)
    } else {
        plot <- plot + ggplot2::scale_y_continuous(limits = c(0, max.value))
    }

    plot <- plot + ggplot2::ylab(legendy) +
                   ggplot2::xlab("biological conditions") +
                   ggplot2::theme_bw() +
                   ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 0, vjust = 0.5),
                                  legend.text = ggplot2::element_text(size = 6),
                                  legend.key = ggplot2::element_blank())

    if (show.on_device) {
        grid::grid.draw(plot)
    }

    invisible(plot)

}


#' @title Visualization of cluster enrichment profiles kinetics
#' 
#' @description 
#' The Kinetics Viewer aims to visualize the cell cluster abundances in a kinetics manner. 
#' This representation displays the cell abundances over the time for each individual. 
#' The cell cluster abundances could be displayed as percentages or absolute numbers using the `use.percentages` parameter (TRUE by default).
#' 
#' @details
#' Timepoints and indivuals information must be assigned to the samples in order to use this viewer (provided by the slot `assignments` in the `Results` object).
#' 
#' @param Results a 'Results' object
#' @param samples a character vector providing the sample names to used (all samples by default)
#' @param clusters a character or a character vector containing the cluster to be shown or if multiple clusters which will be merged
#' @param use.percentages a logical specifying if the visualization should be performed on percentage
#' @param show.on_device a logical specifying if the representation will be displayed on device 
#'
#' @return a 'ggplot' object
#'
#' @export
#' 
#' @import grid reshape2 ggplot2
kineticsViewer <- function(Results,
                           samples         = NULL,
                           clusters        = NULL,
                           use.percentages = TRUE,
                           show.on_device  = TRUE) {

    if (is.null(Results)) {
        stop("Error in kineticsViewer: 'Results' parameter can not be NULL")
    } else if (class(Results)[1] != "Results") {
        stop("Error in kineticsViewer: 'Results' parameter must be a 'Results' object")
    }
    
    if (is.null(samples)) {
        samples     <- Results@sample.names
        cluster.abundances <- Results@cluster.abundances
    } else if (!all(samples %in% Results@sample.names)) {
        stop("Error in kineticsViewer: 'samples' parameter must contains only samples names\n Unknown sample names: ",
             paste(setdiff(unique(samples), Results@sample.names), collapse = " "))
    }
    
    data        <- Results@cluster.abundances[, samples, drop = FALSE]
    cluster.abundances <- data

    if (!is.logical(use.percentages)) { stop("Error in kineticsViewer: The 'use.percentages' parameter must be a logical") }

    if (use.percentages) {
        data.percent <- prop.table(as.matrix(data), 2) * 100
        data         <- data.frame(data.percent, check.names = FALSE)
        legendy      <- "% of cells relative to parent"
    } else {
        legendy      <- "# of cells"
    }
    
    if (is.null(clusters)) {
        stop("Error in kineticsViewer: 'clusters' parameter is required")
    } else if (all(clusters %in% Results@cluster.names)) {
        if (typeof(clusters) != "character") {
            stop("Error in kineticsViewer: 'clusters' parameter must be a character vector")
        }
        clusters    <- unique(clusters)
        cluster.abundances <- cluster.abundances[clusters, ]
        data        <- data[clusters, ]
        data        <- apply(data, 2, sum)
    } else {
        stop("Error in kineticsViewer:\nUnknown clusters : ", paste(setdiff(unique(clusters), Results@cluster.names), collapse = " "))
    }

    if (!is.logical(show.on_device)) { stop("Error in kineticsViewer: 'show.on_device' parameter must be a logical") }

    assignments <- Results@assignments
    if (!is.null(assignments)) {
        order       <- unique(assignments$tp)
        assignments <- assignments[samples, , drop = FALSE]
        order       <- intersect(order, unique(assignments$tp))
        
        names.palette  <- unique(assignments$ind)
        palette        <- ggcolors(length(names.palette))
        names(palette) <- names.palette
        
    } else if (all(c("tp", "ind") %in% colnames(assignments))) {
        stop("Error in kineticsViewer: 'assignments' slot must contain the column \"tp\" and \"ind\" in the provided 'Results' object")
    } else {
        stop("Error in kineticsViewer: 'assignments' slot in the provided 'Results' object is required")
    }

    data            <- data.frame(samples = names(data), value = data, ind = assignments[names(data), "ind"])
    data$tp <- assignments[data$samples, "tp"]
    data$tp <- factor(data$tp, levels = order)

    max.value <- max(data$value, na.rm = TRUE)
    max.value <- max.value + max.value * 0.1 + 1

    cells.number <- sum(colSums(cluster.abundances))

    plot <- ggplot2::ggplot(data = data, ggplot2::aes_string(x = "as.factor(tp)", y = "value", group = "ind", color = "ind")) +
            ggplot2::ggtitle(paste("Kinetics Viewer - cluster: ", paste0(clusters, collapse = ", "), " (", format(cells.number, big.mark = " "), " cells) ", sep = "")) +
            ggplot2::geom_line() +
            ggplot2::scale_color_manual(name = "individuals", values = palette) +
            ggplot2::geom_point(na.rm = TRUE) +
            ggplot2::scale_x_discrete(expand = c(0, 0.05))
   if (use.percentages) {
       plot <- plot + ggplot2::scale_y_continuous(limits = c(0, max.value), breaks = round(seq(0, max.value)), minor_breaks = NULL)
   } else {
       plot <- plot + ggplot2::scale_y_continuous(limits = c(0, max.value))
   }

    plot <- plot + ggplot2::ylab(legendy) +
                   ggplot2::xlab("timepoints") +
                   ggplot2::theme_bw() +
                   ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 0, vjust = 0.5),
                                  legend.text = ggplot2::element_text(size = 6),
                                  legend.key  = ggplot2::element_blank())

    if (show.on_device) {
        grid::grid.draw(plot)
    }

    invisible(plot)
}


#' @title Visualization of cluster abundance dynamics
#'
#' @description
#' The Streamgraph Viewer aims to visualize both absolute and relative abundance of clusters across the samples.
#' This representation displays cell abundances using a stacked area graph which is placed around a central axis.
#' 
#' The cell clusters to visualize must be specified using the `clusters` parameter. 
#' Moreover, samples to be represented and their orders can be specified using the `samples` parameter.
#'
#' @details
#' The order of samples in the 'samples' vector correspond to the order where the sample will be displayed.
#'
#' @param Results a 'Results' object
#' @param samples a character vector providing the sample names to used (all samples by default)
#' @param clusters a character vector containing the clusters names to be visualized
#' @param use.relative a logical specifying if the visualization should be performed on relative abundance
#' @param show.on_device a logical specifying if the representation will be displayed on device
#'
#' @return a 'ggplot' object
#' 
#' @export
#'
#' @import ggplot2
#' @importFrom data.table .SD setDT :=
#' @importFrom reshape2 dcast melt
streamgraphViewer <- function(Results,
                              samples        = NULL,
                              clusters       = NULL,
                              use.relative   = FALSE,
                              show.on_device = TRUE) {
    
    if (is.null(Results)) {
        stop("Error in streamgraphViewer: 'Results' parameter can not be NULL")
    } else if (class(Results)[1] != "Results") {
        stop("Error in streamgraphViewer: 'Results' parameter must be a 'Results' object")
    }
    
    data <- Results@cluster.abundances
    
    if (is.null(samples)) {
        data <- data
    } else if (!all(samples %in% Results@sample.names)) {
        stop("Error in streamgraphViewer: 'samples' parameter must contains only samples names\n Unknown sample names: ",
             paste(setdiff(unique(samples), Results@sample.names), collapse = " "))
    } else {
        data <- data[, samples, drop = FALSE]
    }
    
    if (is.null(clusters)) {
        stop("Error in streamgraphViewer: 'clusters' parameter is required")
    } else if (all(clusters %in% Results@cluster.names)) {
        if (typeof(clusters) != "character") {
            stop("Error in streamgraphViewer: 'clusters' parameter must be a character vector")
        }
        clusters <- unique(clusters)
        data     <- data[clusters, ]
    } else {
        stop("Error in streamgraphViewer:\nUnknown clusters : ", paste(setdiff(unique(clusters), Results@cluster.names), collapse = " "))
    }

    if (!is.logical(use.relative)) { stop("Error in streamgraphViewer: 'use.relative' parameter must be a logical") }
    if (!is.logical(show.on_device)) { stop("Error in streamgraphViewer: 'show.on_device' parameter must be a logical") }

    cells.number <- sum(colSums(data))
    
    if (use.relative) {
        data         <- prop.table(as.matrix(data), 2) * 100
        data         <- data.frame(data, check.names = FALSE)
    }
    
    data                  <- cbind(clusters = rownames(data), data)
    melted.data           <- reshape2::melt(data, id = "clusters")
    colnames(melted.data) <- c("clusters", "samples", "values")

    melted.data$clusters <- factor(melted.data$clusters, levels = clusters)
    melted.data          <- melted.data[order(melted.data$clusters, decreasing = TRUE),]
    melted.data          <- melted.data[order(melted.data$samples),]

    dt.melted.data      <- melted.data 
    res <- plyr::ddply(dt.melted.data, "samples", function(df){cumsum(df$values) - (sum(df$values) / 2)})
    res$samples <- NULL
    dt.melted.data$ymax <- as.numeric(t(as.matrix(res)))
    
    res <- plyr::ddply(dt.melted.data, "samples", function(df){df$ymax - df$values})
    res$samples <- NULL
    dt.melted.data$ymin <- as.numeric(t(as.matrix(res)))
    
    res <- plyr::ddply(dt.melted.data, "samples", function(df){format(round(df$ymax - min(df$ymin), 2), big.mark = " ")})
    res$samples <- NULL
    dt.melted.data$label <- as.character(t(as.matrix(res)))
    
    title <- paste("Streamgraph using ", ifelse(use.relative, "relative", "absolute"), " abundance (", format(cells.number, big.mark = " "), " cells)", sep = "")
    plot  <- ggplot2::ggplot(data = dt.melted.data) +
             ggplot2::ggtitle(title) +
             ggplot2::geom_ribbon(ggplot2::aes_string(x = "samples", ymin = "ymin", ymax = "ymax", group = "clusters", fill = "clusters"), color = "grey40", size = 0.1) +
             ggplot2::geom_point(ggplot2::aes_string(x = "samples", y = "ymax", group = "clusters"), shape = 45) +
             ggplot2::geom_text(ggplot2::aes_string(x = "samples", y = "ymax", label = "label"), check_overlap = TRUE, angle = 360, hjust = 1.1, size = 3) +
             ggplot2::xlab("samples") +
             ggplot2::ylab(ifelse(use.relative, "% of cells", "# of cells")) +
             ggplot2::theme_bw() +
             ggplot2::theme(legend.text      = ggplot2::element_text(size = 6),
                            axis.text.x      = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
                            axis.line        = ggplot2::element_blank(),
                            axis.text.y      = ggplot2::element_blank(),
                            axis.ticks       = ggplot2::element_blank(),
                            panel.background = ggplot2::element_blank(),
                            panel.border     = ggplot2::element_blank(),
                            panel.grid.major = ggplot2::element_blank(),
                            panel.grid.minor = ggplot2::element_blank())

    if (show.on_device) {
        plot(plot)
    }

    invisible(plot)
    
}


#' @title Visualization of cluster phenotypes
#' 
#' @description 
#' The Pheno Viewer aims to visualize the phenotype of single cluster or a set of combined clusters.
#' This representation displays median marker expressions of each sample using parallel coordinates. 
#' In such viewer, each line corresponds to a biological sample and lines are positioned on a space where the X-axis represents the cell markers and the Y-axis represents the marker expressions. 
#' If biological conditions have been assigned to the biological samples, then samples belonging to the same condition will be shown using the same color.
#' 
#' Marker expressions from several clusters can be combined by providing one or several cluster names to the `clusters` parameter.
#' In this case, the displayed marker expressions are the means of the median expressions for each cluster of each marker.
#' 
#' Markers used as clustering markers are shown in blue.
#' A dashed line indicates the mean marker expressions of all samples.
#' Importantly, a grey ribbon indicates the bounds of marker expressions in the overall dataset.
#' The visualization can be restricted to specific markers using the `markers` parameter.
#' 
#' @details 
#' The ranges of value between marker bounds (using the 'bounds' slot) will be displayed using a grey ribbon.
#' 
#' The 'show.mean' parameter allows to visualize three kinds of information:
#' \itemize{
#' \item "none" value will show marker median expressions for each selected samples;
#' \item "only" value will show only the mean of median maker expressions for all selected samples (displayed as black dashed line);
#' \item "both" value will show marker median expressions for each selected samples together with the mean of median maker expressions for all selected samples.
#' }
#' 
#' @param Results a Results object
#' @param samples a character vector providing the sample names to used (all samples by default)
#' @param clusters a character vector containing the clusters names to be visualized (by default all clusters will be displayed)
#' @param markers a character vector specifying the markers to be displayed assignments
#' @param show.mean a character specifying if marker means expression should be displayed, possible value are among: "none", "only" or "both"
#' @param show.on_device a logical specifying if the representation will be displayed on device 
#'
#' @return a 'ggplot' object
#'
#' @export
#'
#' @import grid reshape2 ggplot2 
#' @importFrom gtools mixedsort
phenoViewer <- function(Results,
                        samples        = NULL,
                        clusters       = NULL,
                        markers        = NULL,
                        show.mean      = "both",
                        show.on_device = TRUE) {

    if (is.null(Results)) {
        stop("Error in phenoViewer: 'Results' parameter can not be NULL")
    } else if (class(Results)[1] != "Results") {
        stop("Error in phenoViewer: 'Results' parameter must be a 'Results' object")
    }

    if(length(Results@marker.names) == 0){
        stop("Error in phenoViewer: 'Results' object must contain phenotypes")
    }
    
    if (is.null(samples)) {
        samples     <- Results@sample.names
        data        <- Results@cluster.phenotypes
        cluster.abundances <- Results@cluster.abundances
    } else if (!all(samples %in% Results@sample.names)) {
        stop("Error in phenoViewer: 'samples' parameter must contains only samples names\n Unknown sample names: ",
             paste(setdiff(unique(samples), Results@sample.names), collapse = " "))
    } else {
        data        <- subset(Results@cluster.phenotypes, sample %in% samples, drop = FALSE)
        cluster.abundances <- Results@cluster.abundances[, samples, drop = FALSE]
    }
    
    data <- stats::na.omit(data)

    if (is.null(clusters)) {
        stop("Error in phenoViewer: 'clusters' parameter is required")
    } else if (all(clusters %in% Results@cluster.names)) {
        if (typeof(clusters) != "character") {
            stop("Error in phenoViewer: 'clusters' parameter must be a character vector")
        }
        clusters        <- unique(clusters)
        clusters.select <- data[, "cluster"] %in% clusters
        data            <- data[clusters.select,]
        cluster.abundances     <- cluster.abundances[clusters,]
    } else {
        stop("Error in phenoViewer:\nUnknown clusters : ", paste(setdiff(unique(clusters), Results@cluster.names), collapse = " "))
    }
    
    data <- plyr::ddply(data, c("sample"), function(df) {
        apply(df[, 3:ncol(df)], 2, mean, na.rm = TRUE)
    }) 

    if (is.null(markers)) {
        markers <- Results@marker.names
    } else if (all(markers %in% Results@marker.names)) {
        markers <- unique(markers)
        data <- data[, c("sample", markers)]
    } else {
        stop("Error in phenoViewer: Unknown markers :", paste(setdiff(unique(markers), Results@marker.names), collapse = " "))
    }

    if (show.mean != "none" && show.mean != "both" && show.mean != "only") {
        stop("Error in phenoViewer: 'show.mean' parameter must contain only one of these : 'none', 'both' or 'only'")
    }

    if (!is.logical(show.on_device)) { stop("Error in phenoViewer: 'show.on_device' parameter must be a logical") }

    data           <- reshape2::melt(data, id = c("sample"), stringsAsFactors = FALSE)
    colnames(data) <- c("samples", "marker", "value")
    
    names.palette  <- unique(Results@cluster.phenotypes$sample)
    palette        <- ggcolors(length(names.palette))
    names(palette) <- names.palette
    
    assignments <- Results@assignments
    
    if (!is.null(assignments)) {
        
        order       <- unique(assignments$bc)
        assignments <- assignments[samples, , drop = FALSE]
        data$bc <- assignments[data$samples, "bc"]
        order       <- intersect(order, unique(assignments$bc))
        data$bc <- factor(data$bc, levels = order)
        
        names.palette  <- unique(assignments$bc)
        palette        <- ggcolors(length(names.palette))
        names(palette) <- names.palette

    } else if (is.element("bc", colnames(assignments))) {
        warning("Warning in phenoViewer: 'assignments' slot do not contain the column 'bc' in the provided 'Results' object. Consequently, the samples names will be used in remplacement")
    } else {
        warning("Warning in phenoViewer: 'assignments' slot in the provided 'Results' object is absent. Consequently, the samples names will be used in remplacement")
    }
    
    clustering.markers  <- Results@clustering.markers
    ordered.markers    <- c(gtools::mixedsort(clustering.markers),gtools::mixedsort(setdiff(Results@marker.names, clustering.markers)))
    bold.markers       <- ifelse(is.element(ordered.markers, clustering.markers), "bold", "plain")
    colored.markers    <- ifelse(is.element(ordered.markers, clustering.markers), "blue", "black")
    data$marker        <- factor(data$marker, levels = ordered.markers, ordered = TRUE)

    for (i in seq_len(nrow(data))) {
        data[i, "lower.bound"] <- Results@bounds[1, as.character(data[i, "marker"])]
        data[i, "upper.bound"] <- Results@bounds[2, as.character(data[i, "marker"])]
    }
    
    cells.number <- sum(colSums(cluster.abundances))
    
    title    <- paste("Pheno Viewer - cluster: ", paste0(clusters, collapse = ", "), " (", format(cells.number, big.mark = " "), " cells)", sep = "")
    bounds   <- as.numeric(row.names(Results@bounds))
    subtitle <- paste0("Grey ribbon displays from ", (bounds[1] * 100), "% to ", (bounds[2] * 100), "% percentiles of the range expression")
    
    plot <- ggplot2::ggplot(data = data) +
            ggplot2::ggtitle(bquote(atop(.(title), atop(italic(.(subtitle)), ""))))

    if (nrow(data) == 0) {
        data <- data.frame(marker = Results@marker.names)
        plot <- plot + ggplot2::geom_blank(data = data, ggplot2::aes_string(x = "marker"))
    } else {

        max.value <- -1
        min.value <- -1
        
        max.value <- max(c(data$value, data$upper.bound), na.rm = TRUE)
        min.value <- min(c(data$value, data$lower.bound), na.rm = TRUE)

        max.value <- max.value + 0.1 * max.value
        min.value <- min.value + 0.1 * min.value
        
        if (show.mean == "both" || show.mean == "none") {
            color <- ifelse(is.null(assignments), "samples", "bc")
            plot  <- plot + ggplot2::geom_line(ggplot2::aes_string(x = "marker", y = "value", group = "samples", color = color),
                                               size = 0.8,
                                               alpha = 1) +
                     ggplot2::scale_colour_manual(name = ifelse(is.null(assignments), "samples", "biological conditions"), 
                                                   values = palette)
        }
        if (show.mean == "only" || show.mean == "both") {
            means <- plyr::ddply(data,
                                c("marker"),
                                function(df){mean(df$value, na.rm = TRUE)})
            colnames(means) <- c("marker", "means")
            plot <- plot + ggplot2::geom_line(data  = means,
                                              ggplot2::aes_string(x = "marker", y = "means"),
                                              group = 1,
                                              linetype = "dashed",
                                              size  = 1)
        }
        if (show.mean == "only") {
            plot <- plot + ggplot2::theme(legend.position = "none")
        }

        plot <- plot + ggplot2::geom_ribbon(ggplot2::aes_string(x = "as.numeric(marker)", ymin = "lower.bound", ymax = "upper.bound"),
                                            alpha = 0.1,
                                            fill  = "grey20") +
                       ggplot2::scale_y_continuous(limits = c(min.value, max.value), breaks = round(seq(0, max.value, by = 1), 0)) +
                       ggplot2::theme_bw()
    }
    plot <- plot + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, face = bold.markers, color = colored.markers)) +
                   ggplot2::theme(legend.text = ggplot2::element_text(size = 6),
                                  legend.key  = ggplot2::element_blank()) +
                   ggplot2::xlab("markers") +
                   ggplot2::ylab("marker expressions") +
                   ggplot2::guides(col = ggplot2::guide_legend(ncol = 1))

    if (show.on_device) {
        grid::grid.draw(plot)
    }

    invisible(plot)

}


#' @title Visualization of SPADE cluster or sample similarities using MDS
#'
#' @description 
#' Multidimensional Scaling (MDS) methods aim to represent the similarities and differences among high-dimensional objects into a space of a lower number of dimensions, generally in two or three dimensions for visualization purposes.
#' In MDS representations, the Kruskal Stress (KS) indicates the percentage of information lost during the dimensionality reduction process.
#' 
#' The MDS Viewer aims to visualize the similarities between samples or clusters based on their abundances. 
#' In such representation, each dot represents a sample or a cluster and the distances between the dots are proportional to the Euclidean distance between these objects.
#' 
#' The representation space can be specified using the `space` parameter  ("samples" or "clusters").
#' 
#' @details 
#' In the case of "samples" space, biological conditions information can be assigned to the samples (provided by the slot `assignments` in the `Results` object).
#' This parameter must be a dataframe with sample names in row names and 2 other columns specifying the biological conditions and individuals. 
#' 
#' @param Results a 'Results' object
#' @param space a character specifying the space ("clusters" or "samples", "cluster" by default)
#' @param samples a character vector providing the sample names to used (all samples by default)
#' @param clusters a character vector containing the clusters names to be visualized (by default all clusters will be displayed)
#' @param use.percentages a logical specifying if the visualization should be performed on percentage
#' @param dist.method a character string containing the name of the distance measure to use
#' @param show.on_device a logical specifying if the representation will be displayed on device 
#'
#' @return a 'ggplot' object
#'
#' @export
#' 
#' @import MASS ggplot2 ggrepel grDevices
#' @importFrom data.table .SD := data.table
MDSViewer <- function(Results,
                      space           = "clusters",
                      samples         = NULL,
                      clusters        = NULL,
                      use.percentages = TRUE,
                      dist.method     = "euclidean",
                      show.on_device  = TRUE) {

    if (is.null(Results)) {
        stop("Error in MDSViewer: 'Results' parameter can not be NULL")
    } else if (class(Results)[1] != "Results") {
        stop("Error in MDSViewer: 'Results' parameter must be a 'Results' object")
    }
    
    if (is.null(samples)) {
        samples     <- Results@sample.names
    } else if (!all(samples %in% Results@sample.names)) {
        stop("Error in MDSViewer: 'samples' parameter must contain only samples names\n Unknown sample names: ",paste(setdiff(unique(samples), Results@sample.names), collapse = " "))
    }
    
    data <- Results@cluster.abundances[, samples, drop = FALSE]

    if (use.percentages) {
        data.percent <- prop.table(as.matrix(data), 2) * 100
        data         <- data.frame(data.percent, check.names = FALSE)
        legendy      <- "% of cells relative to parent"
    } else {
        legendy      <- "# of cells"
    }
    
    if (is.null(clusters)) {
        clusters <- Results@cluster.names
    } else if (all(clusters %in% Results@cluster.names)) {
        if (typeof(clusters) != "character") {
            stop("Error in MDSViewer: 'clusters' parameter must be a character vector")
        }
        clusters    <- unique(clusters)
        data        <- data[clusters, ]
    } else {
        stop("Error in MDSViewer:\nUnknown clusters: ", paste(setdiff(unique(clusters), Results@cluster.names), collapse = " "))
    }

    if (!is.logical(use.percentages)) { stop("Error in MDSViewer: 'use.percentages' parameter must be a logical") }

    data <- cbind(cluster = rownames(data), data)

    if (space == "samples") {
        
        assignments <- Results@assignments
        if (!is.null(assignments)) {

        } else if (all(c("tp", "ind") %in% colnames(assignments))) {
            warning("warning in MDSViewer: 'assignments' slot not contain the column \"tp\" or/and \"ind\" in the provided 'Results' object. Consequently the samples names will be used in remplacement")
        } else {
            warning("warning in MDSViewer: 'assignments' slot in the provided 'Results' object is absent. Consequently the samples names will be used in remplacement")
        }

        data <- t(data[, colnames(data) != "cluster"])
    } else if (space != "clusters") {
        stop("Error in MDSViewer: 'space' parameter must be equal to 'clusters' or 'samples'")
    } else if (nrow(data) <= 2) {
        stop("Error in MDSViewer: if 'space' parameter is equal to 'clusters', at least 2 clusters must be selected using the 'clusters' parameter")
    }

    if (!is.logical(show.on_device)) { stop("Error in MDSViewer: 'show.on_device' parameter must be a logical") }

    dist   <- dist(data, method = dist.method)
    
    fit    <- MASS::isoMDS(dist, k = 2, trace = FALSE)
    stress <- fit$stress
    fit    <- fit$point

    x       <- fit[, 1]
    y       <- fit[, 2]
    
    data_i = data.frame(x = x,y = y)
    
    min.lim <- min(min(x), min(y)) * 1.1
    max.lim <- max(max(x), max(y)) * 1.1
    lim     <- max(abs(min.lim), abs(max.lim))
    min.lim <- -lim
    max.lim <- lim
    
    if (space == "samples") {

        title    <- "MDS at the sample level"
        subtitle <- paste0("Kruskal Stress : ", format(round(stress, 2), nsmall = 2))
        
        plot <- ggplot2::ggplot() +
                ggplot2::ggtitle(bquote(atop(.(title), atop(italic(.(subtitle)), ""))))

        if (is.null(assignments)) {
            
            data_i$samples <- rownames(data_i)
            plot <- plot + ggrepel::geom_label_repel(data = data_i, ggplot2::aes_string(x = "x", y = "y", label = "samples"), size = 5) +
                           ggplot2::geom_point(data = data_i, ggplot2::aes_string(x = "x", y = "y"), size = 4)
            
            
        } else {
            assignments <- assignments[samples, , drop = FALSE]
            
            data_i$ind  <- assignments[rownames(data_i), "ind"]
            data_i$ind  <- as.factor(data_i$ind)
            
            data_i$bc <- assignments[rownames(data_i), "bc"]
            data_i$bc <- factor(data_i$bc, levels = unique(assignments$bc))
            data.table_i                 <- data.table::data.table(data_i, key = "bc")
            hulls                        <- data.table_i[, .SD[grDevices::chull(x, y)], by = "bc"]
            
            plot <- plot + ggplot2::geom_polygon(data   = hulls,
                                                 ggplot2::aes_string(x = "x", y = "y", group = "bc", fill = "bc"),
                                                 colour = "black",
                                                 alpha  = 0.3) +
                ggplot2::geom_point(data = data_i, ggplot2::aes_string(x = "x", y = "y", colour = "bc", shape = "ind"), size = 4) +
                ggplot2::scale_shape(name = "individuals") +
                ggplot2::scale_colour_hue(name = "biological conditions") +
                ggplot2::scale_fill_hue(name = "biological conditions")
                           
       }

        plot <- plot + ggplot2::geom_hline(yintercept = (min.lim + max.lim) / 2, linetype = "dashed") +
                       ggplot2::geom_vline(xintercept = (min.lim + max.lim) / 2, linetype = "dashed") +
                       ggplot2::xlim(min.lim, max.lim) +
                       ggplot2::ylim(min.lim, max.lim) +
                       ggplot2::coord_fixed() +
                       ggplot2::theme_bw() +
                       ggplot2::theme(panel.background = ggplot2::element_blank(),
                                      panel.border     = ggplot2::element_rect(fill = NA),
                                      axis.text.x      = ggplot2::element_blank(),
                                      axis.text.y      = ggplot2::element_blank(),
                                      axis.title.x     = ggplot2::element_blank(),
                                      axis.title.y     = ggplot2::element_blank(),
                                      panel.grid.minor = ggplot2::element_blank(),
                                      panel.grid.major = ggplot2::element_blank(),
                                      axis.ticks       = ggplot2::element_blank(),
                                      legend.position  = "right",
                                      legend.key       = ggplot2::element_blank())

    } else {
        data_i <- cbind(data_i, clusters = data[, "cluster"])
        data_i$clusters <- as.factor(data_i$clusters)

        title    <- "MDS at the cluster level"
        subtitle <- paste0("Kruskal Stress : ", signif(stress, 2))

        plot <- ggplot2::ggplot(data = data_i) +
                ggplot2::ggtitle(bquote(atop(.(title), atop(italic(.(subtitle)), "")))) +
                ggplot2::geom_hline(yintercept = (min.lim + max.lim) / 2, linetype = "dashed") +
                ggplot2::geom_vline(xintercept = (min.lim + max.lim) / 2, linetype = "dashed") +
                ggplot2::geom_point(ggplot2::aes_string(x = "x", y = "y"),
                                    size = 2) +
                ggrepel::geom_text_repel(ggplot2::aes_string(x = "x", y = "y", label = "clusters"),
                                         size = 5) +
                ggplot2::xlim(min.lim, max.lim) +
                ggplot2::ylim(min.lim, max.lim) +
                ggplot2::coord_fixed() +
                ggplot2::theme_bw() +
                ggplot2::theme(panel.background = ggplot2::element_blank(),
                               panel.border     = ggplot2::element_rect(fill = NA),
                               axis.text.x      = ggplot2::element_blank(),
                               axis.text.y      = ggplot2::element_blank(),
                               axis.title.x     = ggplot2::element_blank(),
                               axis.title.y     = ggplot2::element_blank(),
                               panel.grid.minor = ggplot2::element_blank(),
                               panel.grid.major = ggplot2::element_blank(),
                               axis.ticks       = ggplot2::element_blank(),
                               legend.position  = "none",
                               legend.key       = ggplot2::element_blank())
    }

    if (show.on_device) {
        plot(plot)
    }

    invisible(plot)
}


#' @title Visualization of marker co-expression
#'
#' @description 
#' The Biplot Viewer aims to visualize co-expressions between 2 markers using a biplot representation. 
#' In such representation, each cell is represented by a dot which is positioned in a two-dimensional space corresponding to the marker expressions.
#' 
#' @details 
#' Cells coming from specific clusters or samples can be selected using the `clusters` and `samples` parameters. 
#' Moreover, samples can be displayed independently (default behaviour) or merged.
#' In order to seed-up the computations, the number of cells to display can be down-sampled using the `resample.ratio` parameter. 
#' 
#' @param Results a 'Results' object (with 'flowset' slot not null)
#' @param x.marker a character indicating the marker name of the first dimension
#' @param y.marker a character indicating the marker name of the second dimension
#' @param samples a character vector providing the sample names to used (all samples by default)
#' @param clusters a character vector containing the clusters names to be visualized (by default all clusters will be used)
#' @param resample.ratio a numeric ratio (between 0 and 1) specifying the down-sampling ratio to show less dots (or NULL)
#' @param sample.merge a logical specifying if the selected samples must be merged in a single biplot
#' @param show.on_device a logical specifying if the representation will be displayed on device 
#'
#' @return a 'ggplot' object
#'
#' @export
#'
#' @import grDevices ggplot2 
#' @importFrom flowCore sampleNames
biplotViewer <- function(Results,
                         x.marker,
                         y.marker,
                         samples        = NULL, 
                         clusters       = NULL,
                         resample.ratio = NULL,
                         sample.merge   = FALSE,
                         show.on_device = TRUE) {

    if (is.null(Results)) {
        stop("Error in biplotViewer: 'Results' parameter can not be NULL")
    } else if (class(Results)[1] != "Results") {
        stop("Error in biplotViewer: 'Results' parameter required a 'Results' object")
    }

    if(length(Results@marker.names) == 0){
        stop("Error in biplotViewer: 'Results' object must contain phenotypes")
    }
        
    flowset <- Results@flowset

    if (is.null(flowset)) {
        stop("Error in biplotViewer: 'flowset' slot of the 'Results' object is not loaded, if the slot 'fcs.files' is not null, use the function load.flowSet before using the 'biplotViewer' function")
    }

    if (!is.null(x.marker) && !is.element(x.marker, Results@marker.names)) {
        stop(paste0("Error in biplotViewer: 'x.marker' parameter: ", x.marker, ", is unknown"))
    }
    if (!is.null(y.marker) && !is.element(y.marker, Results@marker.names)) {
        stop(paste0("Error in biplotViewer: 'y.marker' parameter: ", y.marker, ", is unknown"))
    }

    if (is.null(samples)) {
        samples <- Results@sample.names
    } else if (!all(samples %in% Results@sample.names)) {
        stop("Error in biplotViewer: 'samples' parameter must contains only samples names\n Unknown sample names: ",
            paste(setdiff(unique(samples), Results@sample.names), collapse = " "))
    } else {
        samples <- unique(samples)
    }

    if (is.null(clusters)) {
        clusters <- Results@cluster.names
    } else if (all(clusters %in% Results@cluster.names)) {
        if (typeof(clusters) != "character") {
            stop("Error in biplotViewer: 'clusters' parameter must be a character vector")
        }
        clusters <- unique(clusters)
    } else {
        stop("Error in biplotViewer:\nUnknown clusters : ", paste(setdiff(unique(clusters), Results@cluster.names), collapse = " "))
    }

    if (!is.logical(sample.merge)) { stop("Error in biplotViewer: 'sample.merge' parameter must be a logical") }
    if (!is.logical(show.on_device)) { stop("Error in biplotViewer: 'show.on_device' parameter must be a logical") }

    x.data <- c()
    y.data <- c()
    facet  <- c()
    flowset.samples <- flowCore::sampleNames(flowset)

    
    plots <- list()
    for (sample in samples) {

        flowframe <- flowset[[which(sample == flowset.samples), ]]
        exprs     <- flowframe@exprs
        if (!is.null(clusters)) {
            exprs <- subset(exprs, exprs[, "cluster"] %in% clusters)
        }
        x.data <- c(x.data, exprs[, x.marker])
        y.data <- c(y.data, exprs[, y.marker])

        if (!sample.merge) {
            facet <- c(facet, rep(sample, length(exprs[, x.marker])))
        }

    }
    
    if (is.null(samples)) {
        cluster.abundances <- Results@cluster.abundances
    }else{
        cluster.abundances <- Results@cluster.abundances[, samples, drop = FALSE]
    }
    cells.number.by.sample <- colSums(cluster.abundances)

    if (sample.merge) {
        data <- data.frame(x = x.data, y = y.data)
    } else {
        data <- data.frame(x = x.data, y = y.data, facet = paste0(facet, " (" , format(cells.number.by.sample[facet], big.mark = " "), " cells)"))
    }
    
    if (!is.null(resample.ratio)) {
        if (resample.ratio > 0 && resample.ratio < 1) {
            data <- data[sample(nrow(data), round((nrow(data) * resample.ratio))),]
        } else {
            stop("resample.ratio must be > 0 and < 1 or null")
        }
    }
    
    x.max              <- max(data["x"]) * 1.1
    y.max              <- max(data["y"]) * 1.1
    
    colramp          <- grDevices::colorRampPalette(c("yellow", "red"))
    data$cols        <- grDevices::densCols(data$x, data$y, colramp = colramp)
    
    cells.number <- sum(cells.number.by.sample)
    plot <- ggplot2::ggplot(data = data) +
            ggplot2::ggtitle(paste0(" Biplot Viewer (", format(cells.number, big.mark = " "), " cells)", sep = "")) +
            ggplot2::geom_point(ggplot2::aes_string(x = "x", y = "y", colour = "cols"), size = 0.25) +
            ggplot2::stat_density2d(ggplot2::aes_string(x = "x", y = "y"), size = 0.2, colour = "blue", linetype = "dashed") +
            ggplot2::scale_color_identity() +
            ggplot2::xlab(x.marker) +
            ggplot2::ylab(y.marker) +
            ggplot2::coord_cartesian(xlim = c(-1, x.max), ylim = c(-1, y.max)) +
            ggplot2::theme_bw() +
            ggplot2::theme(panel.grid.major = ggplot2::element_line(color = "black", linetype = "dotted"), legend.key = ggplot2::element_blank())
    if (!sample.merge) {
        plot <- plot + ggplot2::facet_wrap(~facet, scales = "free")
    }

    if (show.on_device) {
        plot(plot)
    }

    invisible(plot)
    
}


#' @title Visualization of marker co-expressions
#'
#' @description 
#' The Distogram Viewer aims to visualize the pairwise co-expressions between all markers using a distogram representation. 
#' In such representation, each tile corresponds to the co-expression between 2 markers and is gradient-colored based on the Pearson correlation between the expressions of those 2 markers (as stored in the phenotype matrix).
#' Markers used as clustering markers are shown in blue.
#' 
#' @details 
#' The Pearson correlation is computed based on the marker expressions.
#' The visualization can be restricted to specific clusters, samples and markers by using respectively the `clusters`, `samples` and `markers` parameters.
#'
#' @param Results a 'Results' object
#' @param clusters a character vector containing the clusters names to be use (by default all clusters will be used)
#' @param samples a character vector providing the sample names to used (all samples by default)
#' @param markers a character vector specifying the markers to be displayed
#' @param show.on_device a logical specifying if the representation will be displayed on device 
#'
#' @return a 'ggplot' object.  This object also contains the correlations for each pair of marker (cor attribute).
#' 
#' @export
#' 
#' @import reshape2 ggplot2
distogramViewer <- function(Results,
                            samples        = NULL,
                            clusters       = NULL,
                            markers        = NULL,
                            show.on_device = TRUE) {

    if (is.null(Results)) {
        stop("Error in distogramViewer: 'Results' parameter can not be NULL")
    } else if (class(Results)[1] != "Results") {
        stop("Error in distogramViewer: 'Results' parameter must be a 'Results' object")
    }
    
    if(length(Results@marker.names) == 0){
        stop("Error in distogramViewer: 'Results' object must contain phenotypes")
    }
    
    if (is.null(samples)) {
        data        <- Results@cluster.phenotypes
    } else if (!all(samples %in% Results@sample.names)) {
        stop("Error in distogramViewer: 'samples' parameter must contains only samples names\n Unknown sample names: ",
             paste(setdiff(unique(samples), Results@sample.names), collapse = " "))
    } else {
        data        <- subset(Results@cluster.phenotypes, sample %in% samples, drop = FALSE)
    }
    
    if (is.null(clusters)) {
        
    } else if (all(clusters %in% Results@cluster.names)) {
        if (typeof(clusters) != "character") {
            stop("Error in distogramViewer: 'clusters' parameter must be a character vector")
        }
        clusters        <- unique(clusters)
        clusters.select <- data[, "cluster"] %in% clusters
        data            <- data[clusters.select,]

    } else {
        stop("Error in distogramViewer:\nUnknown clusters : ", paste(setdiff(unique(clusters), Results@cluster.names), collapse = " "))
    }

    data <- data[, -c(1, 2)]
    data <- stats::na.omit(data)

    if (is.null(markers)) {
        markers <- Results@marker.names
    } else if (all(markers %in% Results@marker.names)) {
        markers <- unique(markers)
        data    <- data[, markers]
    } else {
        stop("Error in distogramViewer: unknown markers :", paste(setdiff(unique(markers), Results@marker.names), collapse = " "))
    }
    
    if (!is.logical(show.on_device)) { stop("Error in distogramViewer: 'show.on_device' parameter must be a logical") }

    if (nrow(data) < 2) {
        stop("Error in distogramViewer: can not compute correlations using only one cluster of one sample")
    }
    cormat <- round(stats::cor(data, method = "pearson"), 2)
    dist   <- stats::as.dist(1 - cormat)
    hc     <- stats::hclust(dist)
    cormat <- cormat[hc$order, hc$order]
    cormat[upper.tri(cormat, diag = TRUE)] <- NA
    
    markers          <- colnames(cormat)
    dimnames(cormat) <- NULL
    melted.cormat    <- reshape2::melt(cormat)
    
    clustering.markers <- is.element(markers, Results@clustering.markers)
    bold.markers       <- ifelse(clustering.markers, "bold", "plain")
    colored.markers       <- ifelse(clustering.markers, "blue", "black")

    plot <- ggplot2::ggplot(data = melted.cormat, ggplot2::aes_string(x = "Var1", y = "Var2", fill = "value")) + 
            ggplot2::ggtitle("Distogram Viewer for marker phenotypes correlations") +
            ggplot2::geom_tile(color = "white") +
            ggplot2::scale_fill_gradient2(low = "green", high = "red", mid = "black", 
                                          midpoint = 0, limit = c(-1, 1), na.value = 'white',
                                          name = "Pearson correlation") +
            ggplot2::annotate(geom     = "text",
                              x        = 1:length(markers),
                              y        = 1:length(markers),
                              angle    = -45,
                              size     = 4,
                              label    = markers,
                              hjust    = 1,
                              fontface = bold.markers,
                              color    = colored.markers) +
            ggplot2::coord_fixed() +
            ggplot2::theme(axis.line            = ggplot2::element_blank(),
                           axis.text.x          = ggplot2::element_blank(),
                           axis.text.y          = ggplot2::element_blank(),
                           axis.ticks           = ggplot2::element_blank(),
                           axis.title.x         = ggplot2::element_blank(),
                           axis.title.y         = ggplot2::element_blank(),
                           panel.background     = ggplot2::element_blank(),
                           panel.border         = ggplot2::element_blank(),
                           panel.grid.major     = ggplot2::element_blank(),
                           panel.grid.minor     = ggplot2::element_blank(),
                           plot.background      = ggplot2::element_blank(),
                           legend.justification = c(1, 0),
                           legend.position      = c(0.4, 0.7),
                           legend.direction     = "horizontal",
                           legend.key           = ggplot2::element_blank()) +
            ggplot2::guides(fill = ggplot2::guide_colorbar(barwidth = 7,
                                                           barheight      = 1,
                                                           title.position = "top",
                                                           title.hjust    = 0.5)) 
    
    if (show.on_device) {
        plot(plot)
    }

	rownames(cormat) <- markers
	colnames(cormat) <- markers
	
	plot$cor = cormat
	
    invisible(plot)

}