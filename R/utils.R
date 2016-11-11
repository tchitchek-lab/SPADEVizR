# @title Internal - Generates marker expression scores describing phenotypes
# 
# @description 
# This function is used internally to generate a melted numeric dataframe of discrete expression categories for each marker of each cluster.
# 
# @details 
# The function calculates the mean of median expressions between samples (NA values are removed).
# The function assign calculated mean of median expressions to a category (the 'num' parameter provides the number of categories) between bounds of marker expressions.
# The bound parameter must contain the marker name in colnames and two rows. 
# For each marker, first row is the lower bound value et second row is the upper bound value.
# The resulting matrix of this function contains 3 columns: "cluster", "marker" and "value".
# 
# @param cluster.phenotypes a dataframe containing the marker median expressions for each cluster of each sample
# @param bounds a dataframe containing the bounds for each marker
# @param num a numeric value specifying the number of markers expression categories
#  
# @return a numeric matrix of expression scores
# 
#' @import plyr
computePhenoTable <- function(cluster.phenotypes, bounds, num = 5){

    cluster.phenotypes        <- stats::na.omit(cluster.phenotypes)
    cluster.phenotypes.melted <- reshape2::melt(cluster.phenotypes, id.vars = c("sample", "cluster"))
    
    colnames(cluster.phenotypes.melted) <- c("sample", "cluster", "marker", "value")
    cluster.phenotypes.melted$marker    <- as.vector(cluster.phenotypes.melted$marker)
    means                               <- plyr::ddply(cluster.phenotypes.melted,c("cluster", "marker"),function(df){mean(df$value, na.rm = TRUE)})
    
    colnames(means)       <- c("cluster", "marker", "value")

    for (i in seq_len(nrow(means))) {
        
        cluster <- means[i, "cluster"]
        value   <- means[i, "value"]

        min     <- bounds[1, means[i, "marker"]]
        max     <- bounds[2, means[i, "marker"]]

        seq               <- seq(from = min, to = max, length.out = num)
        means[i, "value"] <- which.min(abs(value - seq))
        
    }

    return(means)
}


# @title Internal - Create a list of elements allowing to build a heatmap
#
# @description 
# This function is used internally to build the elements needed for an heatmap.
# Clustering markers are displyed in blue.
# 
# Dendrograms are computed based on the Euclidean metrix and the Ward linkage method (by default).
# Others distances could be specified using the distance parameter among : "maximum", "manhattan", "canberra", "binary" or "minkowski".
# 
# @param matrix a numeric matrix containing the markers expression categories
# @param dendrogram.type a character specifying the look of dendrograms ("rectangle" or "triangle", "rectangle" by default)
# @param num a numeric value specifying the number of markers expression categories
# @param xlab a character specifying the X-axis label
# @param ylab a character specifying the Y-axis label
# @param legend.title a character specifying the legend title
# @param clustering.markers a character vector of clustering markers
# @param dendrograms a character specifying if "row", "col", "both" or "none" dendrograms must be build
# @param method a character specifying the agglomeration method used to compute the hierarchical dendrograms
# @param distance a character specifying the measure of distances to be used
# @param tile.color a character specifying the border color of the tiles (NA to remove tile borders)
# @param ... further parameters passed to the R dist method
#
# @return a list of 3 plots (top dendrogram, right dendrogram, heatmap) and the structures of row and column dendrograms (row.hc and col.hc)
#
#' @import ggplot2 grid reshape2 grDevices
ggheatmap <- function(matrix,
                      dendrogram.type    = "rectangle",
                      num                = 5,
                      xlab               = "clusters",
                      ylab               = "markers",
                      legend.title       = "relative expression",
                      clustering.markers = NULL,
                      dendrograms        = "both",
                      method             = "ward.D",
                      distance           = "euclidean",
                      tile.color         = "black",
					  ...) {

    if (dendrograms == "both" || dendrograms == "row") {
        row.hc     <- stats::hclust(stats::dist(matrix, method = distance, ...), method = method)
        row.dendro <- ggdendro::dendro_data(stats::as.dendrogram(row.hc), type = dendrogram.type)
        row.plot   <- ggdendrogram(row.dendro, row = TRUE)
        row.ord    <- match(row.dendro$labels$label, rownames(matrix))
    } else {
		row.hc   <- NULL 
        row.ord  <- rownames(matrix)
        row.plot <- grid::rectGrob(gp = grid::gpar(lwd = 0, col = 0))
    }
        
    if (dendrograms == "both" || dendrograms == "col") {
        col.hc     <- stats::hclust(stats::dist(t(matrix)), "ward.D")
        col.dendro <- ggdendro::dendro_data(stats::as.dendrogram(col.hc), type = dendrogram.type)
        col.plot   <- ggdendrogram(col.dendro, col = TRUE)
        col.ord    <- match(col.dendro$labels$label, colnames(matrix))
    } else {
		col.hc   <- NULL 
        col.ord  <- colnames(matrix)
        col.plot <- grid::rectGrob(gp = grid::gpar(lwd = 0, col = 0))
    }
    
    mat.ordered <- matrix[row.ord, col.ord]

    data.frame  <- as.data.frame(mat.ordered)

    data.frame[data.frame == "-1"] <- NA
    
    data.frame$markers   <- rownames(mat.ordered)
    data.frame$markers   <- factor(data.frame$markers,levels=unique(data.frame$markers),ordered=TRUE)

    melted.data.frame    <- reshape2::melt(data.frame, id.vars = "markers")
        
    colfunc <- grDevices::colorRampPalette(c("#FFFFFF", "#ECE822", "#F9A22B", "#EE302D", "#A32D33"))

    melted.data.frame$value <- as.factor(melted.data.frame$value)
    
    centre.plot <- ggplot2::ggplot(melted.data.frame, ggplot2::aes_string(x = "variable", y = "markers")) + 
                   ggplot2::geom_tile(ggplot2::aes_string(fill = "value"), colour = tile.color) +
                   ggplot2::scale_fill_manual(values = colfunc(num), na.value = "grey50", guide = ggplot2::guide_legend(title          = legend.title,
                                                                                                                        direction      = "horizontal",
                                                                                                                        ncol           = 5,
                                                                                                                        byrow          = TRUE,
                                                                                                                        label.theme    = ggplot2::element_text(size = 10, angle = 0), 
                                                                                                                        label.position = "bottom",
                                                                                                                        label.hjust    = 0.5,
                                                                                                                        title.position = "top")) +
                   ggplot2::theme(legend.text       = ggplot2::element_text(size = 4),
                                  panel.background  = ggplot2::element_rect("white"),
                                  axis.text.x       = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
                                  legend.position   = c(ifelse(num >= 5, 0.6, 1 - (num * 0.1)), 0.5),
                                  legend.background = ggplot2::element_blank()) +
                   ggplot2::xlab(xlab) +
                   ggplot2::ylab(ylab)
                   

    if (!is.null(clustering.markers)) {
        clustering.markers <- is.element(data.frame$markers, clustering.markers)
        bold.markers       <- ifelse(clustering.markers, "bold", "plain")
        colored.markers    <- ifelse(clustering.markers, "blue", "black")
        centre.plot        <- centre.plot + ggplot2::theme(axis.text.y = ggplot2::element_text(face = bold.markers, color = colored.markers))
    }

    ret <- list(col = col.plot, row = row.plot, centre = centre.plot, row.hc = row.hc, col.hc = col.hc) 

    return(ret)
    
}


# @title Internal - Build a dendrogram plot
#
# @description 
# This function is used internally to generate a 'ggplot' dendrogram.
#
# @details 
# It is to note that 'row' and 'col' are mutuality excluded (both cannot be both TRUE) with priority to row.
# 
# @param dist a numeric matrix containing distances between objects
# @param row a logical value specifying if the horizontal dendrogram must be computed
# @param col a logical value specifying if the vertical dendrogram must be computed
# 
# @return a 'ggplot' dendrogram object
#
#' @import ggplot2 ggdendro
ggdendrogram <- function(dist, row=!col, col=!row) {

    p <- ggplot2::ggplot() +
         ggplot2::geom_segment(data = ggdendro::segment(dist),
                               ggplot2::aes_string(x = "x", y = "y", xend = "xend", yend = "yend")) +
         ggplot2::labs(x = NULL, y = NULL) +
         ggdendro::theme_dendro() +
         ggplot2::theme(axis.line        = ggplot2::element_blank(),
                        axis.text.x      = ggplot2::element_blank(),
                        axis.text.y      = ggplot2::element_blank(),
                        axis.ticks       = ggplot2::element_blank(),
                        axis.title.x     = ggplot2::element_blank(),
                        axis.title.y     = ggplot2::element_blank(),
                        legend.position  = "none",
                        panel.background = ggplot2::element_blank(),
                        panel.border     = ggplot2::element_blank(),
                        panel.grid.major = ggplot2::element_blank(),
                        panel.grid.minor = ggplot2::element_blank(),
                        plot.background  = ggplot2::element_blank(),
                        plot.margin      = grid::unit(c(0,0,0,0), "cm"),
                        panel.margin     = grid::unit(c(0,0,0,0), "cm"))
    
    if (row) {
        p <- p + ggplot2::scale_x_continuous(expand = c(0, 0.55)) +
                 ggplot2::coord_flip()                
    } 
    else {
        p <- p +
                ggplot2::scale_x_continuous(expand = c(0, 0.55))
    }
    return(p)
}


# @title Internal - Extraction of a ggplot element
#
# @description 
# This function is used internally to extract an element from a 'ggplot' objet.
#
# @details 
# Example of valid names are : "guide-box", "axis-b", "xlab", axis-l", "ylab"
#
# @param gplot a 'ggplot' plot
# @param name a character specifying the name of the element to be extracted
# 
# @return a 'ggplot' axis object
#
#' @import ggplot2 gtable
ggextract <- function(gplot, name) {
    built <- ggplot2::ggplot_build(gplot)
    tmp   <- ggplot2::ggplot_gtable(built)
    tmp   <- gtable::gtable_filter(tmp, name)
    return(tmp$grobs[[TRUE]])
}


# @title Internal - Generates an heatmap by assembling elements
#
# @description 
# This function is used internally to displays the heatmap elements build by 'ggheatmap()'
#
# @details
# Example of further parameters passed to arrangeGrob are: "top", "rigth", "left", "bottom" which allow to add text at thoses places.
#
# @param list the list of ggplot object provided by ggheatmap
# @param col.width size of horizontal dendrogram
# @param row.width size of vertical dendrogram
# @param ... further parameters passed to arrangeGrob
# 
# @return a ggplot2 axis
#
#' @import ggplot2 grid gridExtra
ggheatmap.combine <- function(list, col.width=0.15, row.width=0.15, ...) {

    layout <- rbind(c(2, 1, NA),
                    c(5, 3, 4),
                    c(NA, 6, NA))

    legend <- ggextract(list$centre, name = "guide-box")

    x.axis       <- ggextract(list$centre, name = "axis-b")
    x.axis.title <- ggextract(list$centre, name = "xlab")
    x.axis       <- gridExtra::arrangeGrob(x.axis, x.axis.title, nrow = 2)

    y.axis       <- ggextract(list$centre, name = "axis-l")
    y.axis.title <- ggextract(list$centre, name = "ylab")
    y.axis       <- gridExtra::arrangeGrob(y.axis.title, y.axis, ncol = 2)

    center.without_legend = list$centre + ggplot2::theme(axis.line        = ggplot2::element_blank(),
                                                         axis.text.x      = ggplot2::element_blank(),
                                                         axis.text.y      = ggplot2::element_blank(),
                                                         axis.ticks       = ggplot2::element_blank(),
                                                         axis.title.x     = ggplot2::element_blank(),
                                                         axis.title.y     = ggplot2::element_blank(),
                                                         legend.position  = "none",
                                                         panel.background = ggplot2::element_blank(),
                                                         panel.border     = ggplot2::element_blank(),
                                                         panel.grid.major = ggplot2::element_blank(),
                                                         panel.grid.minor = ggplot2::element_blank(),
                                                         plot.background  = ggplot2::element_blank(),
                                                         plot.margin      = grid::unit(c(0, 0, 0, 0), "cm"),
                                                         panel.margin     = grid::unit(c(0, 0, 0, 0), "cm"))

    ret <- gridExtra::arrangeGrob(list$col, #1 on the layout
                                  legend, #2 on the layout
                                  center.without_legend, #3 on the layout
                                  list$row, #4 on the layout
                                  y.axis, #5 on the layout
                                  x.axis, #6 on the layout
                                  layout_matrix = layout,
                                  widths        = grid::unit(c(col.width * 0.5, 1 - (2 * col.width), col.width), "null"),
                                  heights       = grid::unit(c(row.width, 1 - (2 * row.width), row.width * 0.5), "null"),
                                  ...)

    return(ret)
}


# @title Internal - Generate a circle representation
#
# @description 
# This function is used internally to generate a packed circles representation
# 
# @param circles a 2 column dataframe the clusters to be displayed and theirs sizes
# @param class a numeric specifying the class number to be displayed
# @param color a character specifying the colour of the packed circles representation
# @param npoint a numeric specifying the levels of details of circles
# @param limits a numeric specifying the size of the coordinate system centered on (0,0)
# @param maxiter a numeric specifying the maximal number of iterations to perform
#
# @return a ggplot2 object
#
#' @import ggplot2 ggrepel grid gridExtra packcircles
buildCircles <- function(circles,
                         color   = "grey80",
                         class   = NA,
                         npoint  = 100,
                         limits  = 30000,
                         maxiter = 100) {
    set.seed(42)
    xyr <- data.frame(x = stats::runif(nrow(circles), 0, 1),
                      y = stats::runif(nrow(circles), 0, 1),
                      r = sqrt(circles$size) * 50)

    res  <- packcircles::circleLayout(xyr, xlim = c( -limits, limits), ylim = c( -limits, limits), maxiter = 1000, wrap = FALSE)
    data <- packcircles::circlePlotData(layout = res$layout, npoints = npoint)
    text <- cbind(res$layout, cluster = circles$cluster)

    set.seed(42)
    plot <- ggplot2::ggplot(data = data) +
            ggplot2::ggtitle(paste0("Class ", class)) +
            ggplot2::geom_polygon(ggplot2::aes_string(x = "x", y = "y", group = "id"),
                                  fill  = color,
                                  color = "grey40",
                                  alpha = 0.2) +
            ggrepel::geom_text_repel(data          = text, ggplot2::aes_string(x = "x", y = "y", label = "cluster"), size = 3,
                                     box.padding   = grid::unit(0.35, "lines"),
                                     point.padding = grid::unit(0.3, "lines")) +
            ggplot2::coord_equal(xlim = c( -limits, limits), ylim = c( -limits, limits)) +
            ggplot2::theme(axis.line        = ggplot2::element_blank(),
                           axis.text.x      = ggplot2::element_blank(),
                           axis.text.y      = ggplot2::element_blank(),
                           axis.ticks       = ggplot2::element_blank(),
                           axis.title.x     = ggplot2::element_blank(),
                           axis.title.y     = ggplot2::element_blank(),
                           legend.position  = "none",
                           panel.background = ggplot2::element_blank(),
                           panel.border     = ggplot2::element_blank(),
                           panel.grid.major = ggplot2::element_blank(),
                           panel.grid.minor = ggplot2::element_blank(),
                           plot.background  = ggplot2::element_blank(),
                           plot.margin      = grid::unit(c(0, 0, 0, 0), "cm"),
                           panel.margin     = grid::unit(c(0, 0, 0, 0), "cm"))

    return(plot)

}


# @title Internal - Generate a legend for circles representation
#
# @description 
# This function is used internally to generate the legend of a packed circles representation
#
# @param circles a 3 columns data frame with the x, y coordinate of points and their radius
# @param npoint a numeric specifying the levels of details of polygons
# @param limits a numeric specifying the size of the coordinate system centered on (0,0)
# 
# @return a ggplot2 object
#
#' @import ggplot2 ggrepel grid gridExtra packcircles
buildCirclesLegend <- function(circles = data.frame(x = c(-29500, -19000, -8000, 3000, 20000),
                                                    y = c(20000, 20000, 20000, 20000, 20000),
                                                    r = c(500, 1000, 2000, 5000, 10000)),
                               npoint  = 100,
                               limits  = 30000) {
        
    text           <- circles
    colnames(text) <- c("x", "y", "cluster")

    circles$r <- sqrt(circles$r) * 50
    data      <- packcircles::circlePlotData(layout = circles, npoints = npoint)

    plot <- ggplot2::ggplot(data = data) +
            ggplot2::ggtitle("number of cells") +
            ggplot2::geom_polygon(ggplot2::aes_string(x = "x", y = "y", group = "id"),
                                  fill  = "white",
                                  color = "grey40",
                                  alpha = 0.2) +
            ggplot2::geom_text(data = text, ggplot2::aes_string(x = "x", y = "y-20000", label = "cluster"), size = 3) +
            ggplot2::coord_equal(xlim = c( -limits, limits), ylim = c( -limits, limits)) +
            ggplot2::theme(axis.line        = ggplot2::element_blank(),
                           axis.text.x      = ggplot2::element_blank(),
                           axis.text.y      = ggplot2::element_blank(),
                           axis.ticks       = ggplot2::element_blank(),
                           axis.title.x     = ggplot2::element_blank(),
                           axis.title.y     = ggplot2::element_blank(),
                           legend.position  = "none",
                           panel.background = ggplot2::element_blank(),
                           panel.border     = ggplot2::element_blank(),
                           panel.grid.major = ggplot2::element_blank(),
                           panel.grid.minor = ggplot2::element_blank(),
                           plot.background  = ggplot2::element_blank(),
                           plot.margin      = grid::unit(c(0, 0, 0, 0), "cm"),
                           panel.margin     = grid::unit(c(0, 0, 0, 0), "cm"))

    return(plot)
                       
}


# @title Internal - Removing of factors in a dataframe 
#
# @description 
# This function is used internally to remove factors in a dataframe
#
# @param dataframe a dataframe 
#
# @return a dataframe 
removeFactors <- function(dataframe) {
    if (!is.null(dataframe)) {
        for (i in seq_len(ncol(dataframe))) {
            dataframe[, i] <- as.vector(dataframe[, i])
        }
    }
    return(dataframe)
}


# @title Internal - Generation a palette of color using the ggplot2 color style
#
# @description 
# This function is used internally to generate a palette of color using the ggplot2 color style
#
# @param n number of desired colors
#
# @return a character vector of hexademical colors
#
#' @import grDevices
ggcolors <- function(n = 6){
    h = c(0, 360) + 15
    if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
    grDevices::hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

# @title Internal - Returns the mode of a numeric vector
#
# @description 
# This function is used internally to determine the mode of a numeric vector
#
# @param x a numeric vector
#
# @return a list with 2 numeric values specifying the mode "x" and it associated density "y"
#
computemode <- function(x) {
    den <- stats::density(x, kernel = c("gaussian"))
    return(list(x = den$x[den$y == max(den$y)], y = max(den$y)))
}  