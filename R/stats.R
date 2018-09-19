#' @title Identification of the Abundant Clusters
#' 
#' @description 
#' This function is used to identify the Abundant Clusters (AC), that is to say clusters having cell abundance (absolute or relative) statistically greater than a specific threshold.
#' Abundant Clusters are identified using a one-sample test (parametrized or non-parametrized).
#' P-values can be corrected for multiple comparisons.
#' 
#' @param Results a 'Results' object
#' @param samples a character vector providing the sample names to used
#' @param use.percentages a logical specifying if the computations should be performed on percentage
#' @param method a character specifying the statistical method used to identify the Abundant Clusters. The parameter can take the values "t.test" or "wilcox.test"
#' @param method.adjust a character specifying if the p-values should be corrected using multiple correction methods among: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY" and "fdr" (from 'stats::p.adjust' method) 
#' @param mu a numeric specifying the theoretical value (mu) of the one sample statistical test
#' @param th.pvalue a numeric specifying the p-value threshold
#' 
#' @return a S4 object of class 'AC'
#' 
#' @export
identifyAC <- function(Results,
                       samples,
                       use.percentages = TRUE,
                       method          = "t.test",
                       method.adjust   = NULL,
                       mu              = 0,
                       th.pvalue       = 0.05) {
    
    message("[START] - Identification of Abundant Clusters")

    if (is.null(Results)) {
        stop("Error in identifyAC: The 'Results' parameter can be NULL")
    } else if (class(Results)[1] != "Results") {
        stop("Error in identifyAC: The 'Results' parameter must be a 'Results' object")
    }
    
    if (!all(samples %in% Results@sample.names)) {
        stop("Error in identifyAC: The 'samples' parameter must contains only samples names\n Unknown sample names: ",
             paste(setdiff(unique(samples), Results@sample.names), collapse = " "))
    }

    if (!is.logical(use.percentages)) {
        stop("Error in identifyAC: The 'use.percentages' parameter must be a logical")
    }

    if (is.null(method) || !is.element(method, c("t.test", "wilcox.test"))) {
        stop("Error in identifyAC: The 'method' parameter must be a character among : 't.test' or 'wilcox.test'")
    }

    if (!is.null(method.adjust) && !is.element(method.adjust, c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"))) {
        stop("Error in identifyAC: The 'method.adjust' parameter must be a character among : 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr'")
    }

    if (th.pvalue <= 0) {
        stop("Error in identifyAC: The 'th.pvalue' parameter can not be negative")
    }

    if (mu < 0) {
        stop("Error in identifyAC: The 'mu' parameter can not be strictly negative")
    }

    data <- Results@cluster.abundances[, samples, drop = FALSE]
    cluster.size <- apply(data, 1, sum)

    if (use.percentages) {
        data <- prop.table(as.matrix(data), 2)
        data <- data * 100
    } else {
        data <- data
    }

    pv <- apply(data, 1, function(x) {
        return(do.call(method, args = list(x = x, alternative = "greater", mu = mu))$p.value)
    })

    if (!is.null(method.adjust)) {
        pv <- stats::p.adjust(pv, method = method.adjust)
    }
    pv[is.na(pv)] <- 1
    results <- data.frame(cluster = Results@cluster.names,
                          mean    = apply(data, 1, mean),
                          sd      = apply(data, 1, stats::sd),
                          pvalue  = pv)

    results$significant <- ifelse(results$pvalue < th.pvalue, TRUE, FALSE)

    AC <- methods::new("AC",
                       sample.names    = colnames(data),
                       cluster.size    = cluster.size,
                       use.percentages = use.percentages,
                       method          = method,
                       method.adjust   = ifelse(is.null(method.adjust), "none", method.adjust),
                       mu              = mu,
                       th.pvalue       = th.pvalue,
                       results         = results)

    print(AC)
    message("[END] - Identification of Abundant Clusters")

    return(AC)
}


#' @title Identification of the Differentially Abundant Clusters
#' 
#' @description
#' This function is used to identify differentially abundant clusters, that is to say clusters having abundance (absolute or relative) statistically different between two biological conditions.
#' Differentially Abundant Clusters are identified using a two-sample test (parametrized or non-parametrized).
#' P-values can be corrected for multiple comparisons.
#'
#' @param Results a 'Results' object
#' @param condition1 a character vector providing the sample names defined as the first condition
#' @param condition2 a character vector providing the sample names defined as the second condition
#' @param use.percentages a logical specifying if the computations should be performed on percentage
#' @param method a character specifying the name of the statistical test to use "t.test" or "wilcox.test"
#' @param method.adjust a character specifying if the p-values should be corrected using multiple correction methods among: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY" and "fdr" (from 'stats::p.adjust' method) 
#' @param method.paired a logical indicating if the statistical test must be performed in a paired manner
#' @param th.pvalue a numeric specifying the p-value threshold
#' @param th.fc a numeric specifying the fold-change threshold
#'
#' @return a S4 object of class 'DAC'
#' 
#' @export
identifyDAC <- function(Results,
                        condition1,
                        condition2,
                        use.percentages = TRUE,
                        method          = "t.test",
                        method.adjust   = NULL,
                        method.paired   = FALSE,
                        th.pvalue       = 0.05,
                        th.fc           = 2){
    
    message("[START] - Identification of Differentially Abundant Clusters\n")

    if (is.null(Results)) {
        stop("Error in identifyDAC: The 'Results' parameter can be NULL")
    } else if (class(Results)[1] != "Results") {
        stop("Error in identifyDAC: The 'Results' parameter must be a 'Results' object")
    }
    
    if (!all(condition1 %in% Results@sample.names)) {
        stop("Error in identifyDAC: The 'condition1' parameter must contains only samples names\n Unknown sample names: ",
             paste(setdiff(unique(condition1), Results@sample.names), collapse = " "))
    }

    if (!all(condition2 %in% Results@sample.names)) {
        stop("Error in identifyDAC: The 'condition2' parameter must contains only samples names\n Unknown sample names: ",
             paste(setdiff(unique(condition2), Results@sample.names), collapse = " "))
    }

    if (!is.logical(use.percentages)) {
        stop("Error in identifyDAC: The 'use.percentages' parameter must be a logical")
    }

    if (is.null(method) || !is.element(method, c("t.test", "wilcox.test"))) {
        stop("Error in identifyDAC: The 'method' parameter must be a character among : 't.test' or 'wilcox.test'")
    }

    if (!is.null(method.adjust) && !is.element(method.adjust, c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"))) {
        stop("Error in identifyDAC: The 'method.adjust' parameter must be a character among : 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr'")
    }

    if (!is.logical(method.paired)) {
        stop("Error in identifyDAC: The 'method.paired' parameter must be a logical")
    }

    if (th.pvalue <= 0) {
        stop("Error in identifyDAC: The 'th.pvalue' parameter can not be negative")
    }

    if (th.fc < 0) {
        stop("Error in identifyDAC: The 'th.fc' parameter can not be negative")
    }

    data         <- Results@cluster.abundances
    data.cond1   <- data[, condition1, drop = FALSE]
    data.cond2   <- data[, condition2, drop = FALSE]
    data         <- cbind(data.cond1, data.cond2)
    cluster.size <- apply(data, 1, sum)
    
    if (use.percentages) {
        data   <- prop.table(as.matrix(data), 2)
        data   <- data * 100
    }else{
        data   <- data
    }
    
    s1 <- ncol(data.cond1)
    
    pv <- apply(data, 1, function(x) {
                    return(do.call(method, args = list(x = x[1:s1], y = x[ - (1:s1)], paired = method.paired))$p.value)
                })
    pv[is.na(pv)] <- 1
    if (!is.null(method.adjust)) {
        pv <- stats::p.adjust(pv, method = method.adjust)
    }
    
    fc <- apply(data, 1, function(x){
                    fc <- mean(x[1:s1]) / mean(x[-(1:s1)])
                    if (is.na(fc)) {
                        fc <- 1
                    }
                    if (fc < 1) {
                        fc <- (-1 / fc)
                    }
                    return(fc)
               })

    results <- data.frame(cluster     = Results@cluster.names,
                          mean.cond1  = apply(data[, condition1, drop = FALSE], 1, mean),
                          sd.cond1    = apply(data[, condition1, drop = FALSE], 1, stats::sd),
                          mean.cond2  = apply(data[, condition2, drop = FALSE], 1, mean),
                          sd.cond2    = apply(data[, condition2, drop = FALSE], 1, stats::sd),
                          fold.change = fc,
                          pvalue      = pv)
     
    cond1 <- ifelse(length(substitute(condition1)) == 1, deparse(substitute(condition1)), "condition1")
    cond2 <- ifelse(length(substitute(condition2)) == 1, deparse(substitute(condition2)), "condition2")

    colnames(results) <- c("cluster",
                           paste0("mean.", cond1),
                           paste0("sd.", cond1),
                           paste0("mean.", cond2),
                           paste0("sd.", cond2),
                           "fold.change",
                           "pvalue")

    deparse(substitute(condition2), width.cutoff = 40)

    results$significant <- ifelse(results$pvalue < th.pvalue, TRUE, FALSE)
    results$significant <- ifelse(abs(results$fold.change) > th.fc, results$significant, FALSE)       
    
    thresholds <- c(pvalue = th.pvalue, fc = th.fc)
    
    DAC <- methods::new("DAC",
                        sample.cond1    = colnames(data.cond1),
                        sample.cond2    = colnames(data.cond2),
                        cluster.size    = cluster.size, 
                        use.percentages = use.percentages,
                        method          = method,
                        method.paired   = method.paired,
                        method.adjust   = ifelse(is.null(method.adjust), "none", method.adjust),
                        th.fc           = th.fc,
                        th.pvalue       = th.pvalue,
                        results         = results)
    print(DAC)

    message("[END] - Identification of Differentially Abundant Clusters")
    
    return(DAC)
}


#' @title Identification of the correlation of SPADE cluster with a phenotype
#' 
#' @description 
#' This function is used to identify Correlated Clusters, that is to say clusters having an abundance (absolute or relative) statistically correlated with a biological variable.
#' Correlated Clusters are identified using Pearson or Spearman coefficients of correlation.
#' P-values can be corrected for multiple comparisons.
#'
#' @param Results a 'Results' object
#' @param variable a numerical named vector providing the correspondence between a sample name (in rownames) and the specific numerical phenotype
#' @param use.percentages a logical specifying if the computations should be performed on percentage
#' @param method a character indicating the correlation method to use: "pearson", "spearman"
#' @param method.adjust a character specifying if the p-values should be corrected using multiple correction methods among: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY" and "fdr" (from 'stats::p.adjust' method)
#' @param th.correlation a numeric specifying the absolute value of the correlation coefficient threshold 
#' @param th.pvalue a numeric specifying the p-value threshold
#' 
#' @return a S4 object of class 'CC'
#'
#' @export
identifyCC <- function(Results,
                       variable,
                       use.percentages = TRUE,
                       method          = "pearson",
                       method.adjust   = NULL,
                       th.correlation  = 0.75,
                       th.pvalue       = 0.05){
    
    message("[START] - Identification of Correlated Clusters")

    if (is.null(Results)) {
        stop("Error in identifyCC: The 'Results' parameter can be NULL")
    } else if (class(Results)[1] != "Results") {
        stop("Error in identifyCC: The 'Results' parameter must be a 'Results' object")
    }
    
    if (!all(names(variable) %in% Results@sample.names)) {
        stop("Error in identifyCC: The 'variable' parameter must contains only samples names in names\n Unknown sample names: ",
             paste(setdiff(unique(names(variable)), Results@sample.names), collapse = " "))
    }
    
    if (stats::var(variable) == 0) {
        stop("Error in identifyCC: The 'variable' parameter must have a positive variance")
    }
    
    if (!is.logical(use.percentages)) {
        stop("Error in identifyCC: The 'use.percentages' parameter must be a logical")
    }

    if (is.null(method) || !is.element(method, c("pearson", "spearman", "kendall"))) {
        stop("Error in identifyCC: The 'method' parameter must be a character among : 'pearson', 'spearman' or 'kendall'")
    }

    if (!is.null(method.adjust) && !is.element(method.adjust, c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"))) {
        stop("Error in identifyCC: The 'method.adjust' parameter must be a character among : 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr'")
    }

    if (th.pvalue <= 0) {
        stop("Error in identifyCC: The 'th.pvalue' parameter can not be negative")
    }

    if (th.correlation < 0 || th.correlation > 1) {
        stop("Error in identifyCC: The 'th.correlation' parameter must be included in the domain: [0;1]")
    }

    cluster.abundances  <- Results@cluster.abundances
    variable     <- stats::na.omit(variable) 
    data         <- cluster.abundances[, names(variable), drop = FALSE]
    cluster.size <- apply(data, 1, sum)
    
    if (use.percentages) {
        data   <- prop.table(as.matrix(data), 2)
        data   <- data * 100
    }else{
        data   <- as.matrix(data)
    }
    
    n            <- nrow(data)
    cor.estimate <- vector(mode = "numeric", length = n)
    cor.pvalue   <- vector(mode = "numeric", length = n)
    
    for (i in seq_len(n)) {
        suppressWarnings(cor <- stats::cor.test(variable, data[i, 1:ncol(data)], method = method))
        cor.estimate[i] <- ifelse(!is.na(cor$estimate),cor$estimate,0)
        cor.pvalue[i]   <- ifelse(!is.na(cor$p.value),cor$p.value,1)
    }
    
    if (!is.null(method.adjust)) {
        cor.pvalue <- stats::p.adjust(cor.pvalue, method = method.adjust)
    }
    
    results <- data.frame(cluster     = Results@cluster.names,
                          correlation = cor.estimate,
                          pvalue      = cor.pvalue)
    
    results$significant <- ifelse(results$pvalue < th.pvalue, TRUE, FALSE)
    results$significant <- ifelse(abs(results$correlation) > th.correlation, results$significant, FALSE)      
    
    CC <- methods::new("CC",
                       sample.names    = colnames(data),
                       variable        = variable,
                       cluster.size    = cluster.size,
                       use.percentages = use.percentages,
                       method          = method,
                       method.adjust   = ifelse(is.null(method.adjust), "none", method.adjust),
                       th.correlation  = th.correlation,
                       th.pvalue       = th.pvalue,
                       results         = results)

    print(CC)

    message("[END] - Identification of Correlated Clusters")
    
    return(CC)
}


#' @title Classification of clustering results based on the abundance profiles
#' 
#' @description 
#' Classifies clusters based on their abundance profiles (number of cells for each cluster).
#' 
#' @details 
#' The classification is done on cell abundances of each clusters and could be performed using 5 methods:
#' \itemize{
#' \item "hierarchical_k" 
#' This method first compute the Pearson correlation matrix and then use this matrix to performs a hierarchical classification. 
#' The hierarchical classification is cut in order to return the desired number of classes. 
#' This number of classes must be provided as a numeric integer using the 'method.parameter' parameter.
#' It is to note that negative correlations are considered as uncorrelated
#' \item "hierarchical_h" (default method)
#' This method works in the same way than 'hierarchical_k' but the height where the hierarchical tree is specified. 
#' This height is a correlation threshold (a numeric double between 0 and 1 included, default is 0.7) provided using the 'method.parameter' parameter.
#' \item "kmeans"
#' This method works as described in the R stats documentation (?kmeans) using the 'method.parameter' parameter to specify the desired number of classes.
#' \item "eigencell" 
#' This method performs an eigen vector decomposition and then calculate the correlations between cluster values and these vectors.
#' Clusters which correlate above a specific threshold with the same eigen vector are classified together.
#' This correlation threshold (a numeric double between 0 and 1 included, default is 0.8) provided using the 'method.parameter' parameter.
#' \item "clique" 
#' This method first computes the Pearson correlation matrix and then use this matrix to generate an undirected graph.
#' In this graph, an edge is drawn between two nodes if the correlation coefficient in the adjacency matrix is above a specific threshold. 
#' This correlation threshold (a numeric double between 0 and 1 included, default is 0.7) provided using the 'method.parameter' parameter.
#' After building the graph, the method looking for the largest cliques which are considered as classes of nodes. Cliques correspond to subgraph in which every two distinct vertices are adjacent.
#' }
#' 
#' @param Results a 'Results' object
#' @param method a character specifying the clustering method among: "hierarchical_h", "hierarchical_k","k-means","eigencell","clique"
#' @param method.parameter a numeric specifying the numeric value required by the selected method 
#' @param use.percentages a logical specifying if cell cluster abudances must be expressed as percentages
#' 
#' @return a S4 object of class 'AP'
#'
#' @export
classifyAbundanceProfiles <- function(Results,
                                      method           = "hierarchical_h",
                                      method.parameter = NULL,
									  use.percentages  = FALSE){

    default.eigencell.correlation.th    <- 0.8
    default.clique.correlation.th       <- 0.7
    default.hierarchical.correlation.th <- 0.7                           
    
    message("[START] - computing classifyClusteringResults")
    
    if (is.null(Results)) {
        stop("Error in classifyAbundanceProfiles: The 'Results' parameter can be NULL")
    } else if (class(Results)[1] != "Results") {
        stop("Error in classifyAbundanceProfiles: The 'Results' parameter must be a 'Results' object")
    }
    
    
    if (!is.element(method, c("hierarchical_h", "hierarchical_k", "k-means", "eigencell", "clique"))) {
        stop("Error : In classifyAbundanceProfiles, method must be among: 'hierarchical_h','hierarchical_k','k-means','eigencell','clique'")
    }

    data <- Results@cluster.abundances
	
	if (use.percentages) {
        data   <- prop.table(as.matrix(data), 2)
        data   <- data * 100
    }else{
        data   <- as.matrix(data)
    }

    switch(method,
            "hierarchical_h" = {
                if (is.null(method.parameter)) {
                    method.parameter <- default.hierarchical.correlation.th
                    warning("Default hierarchical correlation threshold will be used (", method.parameter, ")")
                }
                classes <- computeHierarchicalClustering(data, class.number = NULL, hierarchical.correlation.th = method.parameter)
            },
            "hierarchical_k" = {
                if (is.null(method.parameter)) {
                    stop("Error in classifyAbundanceProfiles: The parameter 'method.parameter' must contain an integer with 'hierarchical_k' method")
                }
                classes <- computeHierarchicalClustering(data, class.number = method.parameter, hierarchical.correlation.th = NULL)
            },
            "k-means" = {
                if (is.null(method.parameter)) {
                    stop("Error in classifyAbundanceProfiles: The parameter 'method.parameter' must contain an integer with 'k-means' method")
                }
                classes <- computeKmeans(data, method.parameter)
            },
            "eigencell" = {
                if (is.null(method.parameter)) {
                    method.parameter <- default.eigencell.correlation.th
                    warning("Default eigencell correlation threshold will be used (", method.parameter, ")")
                }
                classes <- computeEigenCellClusters(data, method.parameter)
            },
            "clique" = {
                if (is.null(method.parameter)) {
                    method.parameter <- default.clique.correlation.th
                    warning("Default clique correlation threshold will be used (", method.parameter, ")")
                }
                classes <- computeClique(data, method.parameter)
            }, stop("Error in classifyAbundanceProfiles: The parameter 'method' must be among 'hierarchical_h', 'hierarchical_k', 'k-means', 'eigencell' or 'clique'"))

    classes$class <- as.numeric(classes$class)
    classes       <- classes[order(classes$class),]

    cluster.size        <- apply(data, 1, sum)
    names(cluster.size) <- rownames(data)
    
    AP <- methods::new("AP",
                       class.number     = length(unique(classes[!is.na(classes$class), 2])),
                       cluster.size     = cluster.size,
                       method           = method,
                       method.parameter = method.parameter,
                       classes          = classes)
    print(AP)

    message("[END] - computing classifyClusteringResults")
    return(AP)
}


# @title Internal - Hierarchical clustering classification
# 
# @description 
# This function is used internally to classify clusters abundance profiles or phenotype profiles using a hierarchical algorithm. 
# 
# @details 
# This function compute the Pearson correlation matrix associated to the provided matrix. 
# It is to note that negative correlations are considered as uncorrelated.
# This correlation matrix is used to performs a hierarchical classification.
# If 'class.number' parameter is NULL, classification will be determined based on the cut height correlation threshold (i.e. 'hierarchical.correlation.th' parameter)
# 
# @param data a matrix with all clusters in rownames
# @param class.number a numeric specifying the number of classes
# @param hierarchical.correlation.th a numeric value specifying the cut height
# 
# @return a dataframe containing for each cluster, its name and class
computeHierarchicalClustering <- function(data,
                                          class.number                = NULL,
                                          hierarchical.correlation.th = 0.8){
    
    cor.data <- stats::cor(t(data))
    
    cor.data[cor.data < 0] <- 0
    cor.data               <- 1 - cor.data
    
    dist.cor.data <- stats::as.dist(cor.data)
    tree          <- stats::hclust(dist.cor.data)
    
    if (!is.null(class.number)) {
        res <- stats::cutree(tree, k = class.number)
    }else{
        res <- stats::cutree(tree, h = hierarchical.correlation.th)
    }
    
    res <- data.frame(cluster = as.character(names(res)), class = res)
    
    return(res)
}


# @title Internal - Kmeans classification
# 
# @description 
# This function is used internally to classify clusters abundance profiles or phenotype profiles using a k-means algorithm. 
# 
# @details 
# This method works as described in the R stats documentation (?kmeans) using the 'k' parameter to specify the desired number of classes.
#
# @param data a numeric matrix with cluster names in rownames
# @param k a numeric specifying the desired number of classes
# 
# @return a dataframe containing for each cluster, its name and class
computeKmeans <- function(data,
                          k = NULL){
    kmeans  <- stats::kmeans(data, centers = k)    
    results <- data.frame(cluster = as.character(rownames(data)), class = as.numeric(kmeans$cluster))
    
    return(results)
}


# @title Internal - Eigen-vector classification
# 
# @description 
# This function is used internally to classify clusters abundance profiles or phenotype profiles using eigen vector decomposition. 
# 
# @details 
# This method computes an eigen vector decomposition and then calculate the correlations between the matrix rows and these vectors.
# Clusters which correlate above a specific threshold with the same eigen vector are classified together.
# This correlation threshold (a numeric double between 0 and 1 included, default is 0.8) provided using the 'eigencell.correlation.th' parameter.
#
# @param data a numeric matrix with all clusters in rownames
# @param eigencell.correlation.th a numeric value indicating the correlation coefficient threshold
# 
# @return a dataframe containing for each cluster, its name and class
computeEigenCellClusters <- function(data, 
                                     eigencell.correlation.th = 0.80){

    

    svd     <- svd(data)
    eigenCC <- t(svd$v)
    
    res <- data.frame(stringsAsFactors = FALSE)
    for (i in seq_len(eigenCC)) {
        for (j in seq_len(nrow(data))) {
            cor <- stats::cor(as.numeric(eigenCC[i, ]), as.numeric(data[j, ]), method = "pearson")
            if (cor > eigencell.correlation.th) {
                res <- rbind(res, cbind(as.character(rownames(data[j,])), i))
            }
        }
    }

    if (nrow(res) > 0) { 
        colnames(res) <- c("cluster", "class")
        classes.uniq  <- unique(res[, "class"])
        classes       <- data.frame(class = classes.uniq, renumbered = 1:length(classes.uniq))
        joined.class  <- merge(res, classes, by = "class")
        res           <- joined.class[, c("cluster", "renumbered")]
        colnames(res) <- c("cluster", "class")
        unclassified  <- setdiff(rownames(data), res$cluster)
                
    } else {
        unclassified <- rownames(data)
    }

    res <- rbind(res, data.frame(cluster = unclassified, class = rep(NA, length(unclassified))))
    return(res)
}


# @title Internal - Clique percolation classification
# 
# @description 
# This function is used internally to classify clusters abundance profiles or phenotype profiles using a clique percolation algorithm. 
# 
# @details 
# This method first compute the Pearson correlation matrix and then use this matrix to generate an undirected graph.
# In this graph, an edge is drawn between two nodes if the correlation coefficient in the adjacency matrix is above a specific threshold. 
# This correlation threshold (a numeric double between 0 and 1 included, default is 0.7) provided using the 'clique.correlation.th' parameter.
# After building the graph, the method looking for the largest cliques which are considered as classes of nodes. 
# Cliques correspond to subgraph in which every two distinct vertices are adjacent.
# 
# @param data a numeric matrix with all clusters in rownames
# @param clique.correlation.th a numeric value indicating the correlation coefficient threshold
# 
# @return a dataframe containing for each cluster, its name and class
# 
#' @import igraph
computeClique <- function(data,
                          clique.correlation.th = 0.7){
    
    res <- data.frame()
    for (i in seq_len(nrow(data) - 1)) {
        for (j in (i + 1):nrow(data)) {
            cor <- stats::cor(as.numeric(data[i, ]), as.numeric(data[j, ]), method = "pearson")
            if (cor > clique.correlation.th) {
                res <- rbind(res, cbind(j, i, cor))
            }
        }
    }
    
    if (nrow(res) > 0) {
        colnames(res) <- c("cluster", "cluster", "cor")
        res <- data.frame(res)
        graph <- igraph::graph.data.frame(res, directed = FALSE)
        lists <- igraph::largest.cliques(graph)

        res <- data.frame()

        for (i in seq_len(length(lists))) {
            cluster <- as.character(rownames(data[names(lists[[i]]),]))
            res <- rbind(res, cbind(cluster, i))
        }

        colnames(res) <- c("cluster", "class")
        unclassified <- setdiff(rownames(data), res$cluster)
    } else {
        unclassified <- rownames(data)
    }
    
    res <- rbind(res, data.frame(cluster = unclassified, class = rep(NA, length(unclassified))))
    
    return(res)
    
}

