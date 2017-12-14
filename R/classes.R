#' @title Results class definition
#' 
#' @description 
#' The Results object is a S4 object containing cell clustering results. 
#' 
#' This object mainly stores the cluster abundance matrix (i.e. the number of cells associated to each sample for each cluster) and the cluster phenotypes matrix (i.e. the median expressions for each marker of each cluster).
#'  
#' In addition, this object can contain information about clustering results, such a SPADE tree. 
#' 
#' @details 
#' The 'cluster.abundances' dataframe contains the number of cells associated to each sample for each cluster.
#' This dataframe stores the clusters in rows and the samples in columns.
#' 
#' The 'cluster.phenotypes' dataframe stores the median expressions for each marker of each cluster.
#' This dataframe stores in the first column the sample names, in the second column the cluster names, and in the others columns the maker median expressions.
#' 
#' The 'bounds' dataframe contains the marker expressions boundaries (minimum and maximum, or specific percentiles) for each marker.
#' 
#' The 'print()' and 'show()' can be used to display a summary of this object. 
#' 
#' @slot cluster.abundances a dataframe containing the number of cells associated to each sample for each cluster
#' @slot cluster.phenotypes a dataframe containing the median expressions for each marker of each cluster
#' @slot sample.names a character vector containing the sample names
#' @slot cluster.names a character vector containing the cluster names
#' @slot cluster.number a numeric specifying the number of clusters
#' @slot marker.names a character vector containing the marker names
#' @slot clustering.markers a character vector specifying the markers that have been used by the clustering algorithms
#' @slot bounds a numeric data.frame containing the marker expressions boundaries for each marker
#' @slot use.raw.medians a logical specifying if the marker expressions correspond to raw or transformed data
#' @slot flowset a flowSet object containing the imported SPADE FCS files
#' @slot fcs.files a character vector containing the location of the imported FCS files
#' @slot graph a igraph object containing the SPADE tree structure
#' @slot graph.layout a numeric matrix containing the SPADE tree layout
#' @slot assignments a dataframe containing annotations for each sample samples such as a biological condition ("bc"), a timepoint condition ("tp") or an individual ("ind") assignment 
#' @slot th.min_cells a numeric specifying the minimal number of cells that a cluster for a given samples needs to have to be taken into consideration in its phenotypical characterization
#' 
#' @import igraph methods
#' 
#' @name Results-class
#' @rdname Results-class
#' @exportClass Results
Results <- setClass("Results",
    slots = c(cluster.abundances = "data.frame",
        cluster.phenotypes = "data.frame",
        sample.names       = "character",
        cluster.names      = "character",
        cluster.number     = "numeric",
        marker.names       = "character",
        clustering.markers = "character",
        bounds             = "data.frame",
        use.raw.medians    = "logical",
        flowset            = "ANY",
        fcs.files          = "character",
        graph              = "ANY",
        graph.layout       = "ANY",
        assignments        = "ANY",
        th.min_cells       = "numeric"),
    validity = function(object){
        if ((length(object@marker.names) != 0) && (length(object@marker.names) + 2) != ncol(object@cluster.phenotypes)) {
            message(paste0("Error in Results object: marker.names length (",length(object@marker.names)," + 2 for cluster IDs and sample names) are inconsistent with cluster.phenotypes size (number of columns : ",ncol(object@cluster.phenotypes),")"))
            return(FALSE)
        }
        if (nrow(object@cluster.abundances) != object@cluster.number){
            message(paste0("Error in Results object: cluster.number (",object@cluster.number,") is inconsistent with cluster.abundances matrix size (",nrow(object@cluster.abundances),")"))
            return(FALSE)
            }
        if ((length(object@marker.names) != 0) && nrow(object@bounds) != 2) {
            message(paste0("Error in Results object: bounds number of rows (",nrow(object@bounds),") is incorrect (only 2 rows accepted)"))
            return(FALSE)
            }
        if ((length(object@marker.names) != 0) && ncol(object@bounds) != length(object@marker.names)) {
            message(paste0("Error in Results object: bounds number of columns (",ncol(object@bounds),") is inconsistent with marker.names length (",
            length(object@marker.names), ")"))
			message("It is likely that automatic gating results were generated using FCS containing different cell markers.")
			message("Please use FCS files containing identical set of cell markers.")
            return(FALSE)
        }
        if (object@th.min_cells < 0) {
            message(paste0("Error in Results object: th.min_cells must be positif"))
            return(FALSE)
        }
        if (!is.null(object@assignments) && !is.data.frame(object@assignments)) {
            message(paste0("Error in Results object: assignments must be a dataframe"))
            return(FALSE)
        }
        if (!is.null(object@assignments) && (!all(rownames(object@assignments) %in% object@sample.names))) {
            message(paste0("Error in Results object: assignments must contains all samples in rownames.\n These ones are missing:",paste(setdiff(object@sample.names, rownames(object@assignments)), collapse = " ")))
            return(FALSE)
        }
        if (!is.null(object@flowset) && (class(object@flowset)[1] != "flowSet")) {
            message("Error in Results object: flowset must be of class flowSet or null")
            return(FALSE)
        }
        if (!is.null(object@fcs.files) && (class(object@fcs.files)[1] != "character")) {
            message("Error in Results object: fcs.files must be a character vector or null")
            return(FALSE)
        }
        if (!is.null(object@graph) && (class(object@graph)[1] != "igraph")) {
            message("Error in Results object: graph must be of class igraph or null")
            return(FALSE)
        }
        if (!is.null(object@graph.layout) && (class(object@graph.layout)[1] != "matrix")) {
            message("Error in Results object: graph.layout must be a matrix or null")
            return(FALSE)
        }
        if ((length(object@marker.names) != 0) && (length(object@sample.names) * object@cluster.number) != nrow(object@cluster.phenotypes)) {
            message(paste0("Error in Results object: sample.names length (",length(object@sample.names),") and cluster.number (",object@cluster.number,") are inconsistent with cluster.phenotypes size (number of row : ",nrow(object@cluster.phenotypes),")"))
            return(FALSE)
        }
        if ((length(object@marker.names) != 0) && (length(object@marker.names) + 2) != ncol(object@cluster.phenotypes)) {
            message(paste0("Error in Results object: marker.names length (",length(object@marker.names)," + 2 for cluster IDs and sample names) are inconsistent with cluster.phenotypes size (number of columns : ",ncol(object@cluster.phenotypes),")"))
            return(FALSE)
        }
        if (nrow(object@cluster.abundances) != object@cluster.number) {
            message(paste0("Error in Results object: cluster.number (",object@cluster.number,") is inconsistent with cluster.abundances matrix size (",nrow(object@cluster.abundances),")"))
            return(FALSE)
        }
        if (ncol(object@cluster.abundances) != length(object@sample.names)) {
            message(paste0("Error in Results object: number of samples (",length(object@sample.names),") is inconsistent with cluster.abundances matrix size (",ncol(object@cluster.abundances),")"))
            return(FALSE)
        }
        if (length(object@clustering.markers) > length(object@marker.names)) {
            message(paste0("Error in Results object: clustering.markers length (",length(object@clustering.markers),") can not be higher than marker.names length (",length(object@marker.names),")"))
            return(FALSE)
        }
        if (!all(object@clustering.markers %in% object@marker.names)) {
            message(paste0("Error in Results object: clustering.markers must contains markers included in marker.names (",setdiff(object@clustering.markers, object@marker.names),")"))
            return(FALSE)
        }
        if (length(object@clustering.markers) == 0) {
            warning("Warning in Results object: clustering.markers length is 0")
        }
        for (fcs.file in object@fcs.files) {
            if (!file.exists(fcs.file)) {
                message(paste0("Error in Results object: FCS file not exist :", fcs.file))
                return(FALSE)
            }
        }
        return(TRUE)
    }
)

#' @title Abundant Clusters (AC) class definition
#' 
#' @description 
#' The 'AC' object is a S4 object containing the information related to the abundant clusters in a given biological condition. 
#' Moreover, this object contains parameters and results used in the statistical analysis.  
#' 
#' @details 
#' A cluster is considered as a significant abundant cluster if its associated p-value and mean are below the specific thresholds 'th.pvalue'.  
#' 
#' The 'print()' and 'show()' can be used to display a summary of this object. 
#' Moreover all information about this object could be saved as a tab separated file using the 'export()' method.
#' This object is returned by the 'identifyAC()' function. 
#' 
#' @slot sample.names a character vector containing the samples used to compute the abundant clusters
#' @slot cluster.size a numeric vector containing the number of cells for each cluster
#' @slot use.percentages a logical specifying if computation was performed on percentage of cell abundance
#' @slot method a character containing the name of the statistical test used to identify the abundant clusters
#' @slot method.adjust a character containing the name of the multiple correction method used (if any)
#' @slot mu a numeric specifying the theoretical value of the one sample statistical test
#' @slot th.pvalue a numeric value specifying the p-value threshold
#' @slot results a data.frame containing for each cluster (first column): the mean (second column) and the standard deviation (third column) of the biological condition, the associated p-value (fourth column) and a logical (fifth column) specifying if the cluster is significantly abundant.
#' 
#' @name AC-class
#' @rdname AC-class
#' @exportClass AC
AC <- setClass("AC",
    slots = c(sample.names = "character",
        cluster.size    = "numeric",
        use.percentages = "logical",
        method          = "character",
        method.adjust   = "character",
        mu              = "numeric",
        th.pvalue       = "numeric",
        results         = "data.frame"),
    validity = function(object){
        if (length(object@sample.names) == 0) {
            message("Error in AC object: sample.names length can not be equal to 0")
            return(FALSE)
        }
        if (!any(object@method %in% c("t.test","wilcox.test"))) {
            message("Error in AC object: method must match 't.test' or 'wilcox.test'")
            message(paste0("method founded is : ", object@method))
            return(FALSE)
        }
        if (!any(object@method.adjust %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"))) {
            message("Error in AC object: method.adjust must match 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr', 'none')")
            message(paste0("method.adjust founded is : ", object@method.adjust))
            return(FALSE)
        }
        if (object@th.pvalue < 0 || object@th.pvalue > 1) {
            message("Error in AC object: th.pvalue must be include into [0,1] interval")
            message(paste0("th.pvalue founded is : ", object@th.pvalue))
            return(FALSE)
        }                
        if (!identical(colnames(object@results),c("cluster","mean","sd","pvalue","significant"))) {
            print(colnames(object@results))
            message("Error in AC object: results must have this colmuns : 'cluster','mean','sd','pvalue','significant'")
            message("Colmuns founded are : ")
            message(paste(colnames(object@results)))
            return(FALSE)
        }
        if (object@mu < 0) {
            message("Error in AC object: Error in mu slot, mu can not be strictly negative")
            return(FALSE)
        }
        return(TRUE)
    }
)


#' @title Differentially Abundant Clusters (DAC) class definition
#' 
#' @description 
#' The 'DAC' object is a S4 object containing the information related to the differentially abundant clusters between two given biological conditions. 
#' Moreover, this object contains parameters and results used in the statistical analysis.
#' 
#' @details
#' A cluster is considered as a differentially enriched cluster if its associated p-value and fold-change are below the specific thresholds 'th.pvalue' and 'th.fc'.  
#' 
#' The 'print()' and 'show()' can be used to display a summary of this object.
#' Moreover all information about this object could be saved as a tab separated file using the 'export()' method.
#' This object is returned by the 'identifyDAC()' function. 
#' 
#' @slot sample.cond1 a character specifying the names of the samples of the first biological condition
#' @slot sample.cond2 a character specifying the names of the samples of the second biological condition
#' @slot cluster.size a numeric vector containing the number of cells for each cluster
#' @slot use.percentages a logical specifying if computation was performed on percentage of cell abundance
#' @slot method a character containing the name of the statistical test used to identify the DAC
#' @slot method.adjust a character containing the name of the multiple correction method used (if any)
#' @slot method.paired a logical indicating if the statistical test have been performed in a paired manner
#' @slot th.fc a numeric value specifying the fold-change threshold
#' @slot th.pvalue a numeric value specifying the p-value threshold
#' @slot results a data.frame containing for each cluster (first column): the fold-change (second column) and the standard deviation (third column) for the first biological condition, the fold-change (fourth column) and the standard deviation (fifth column) for the second biological condition, the associated p-value (sixth column) and a logical (seventh column) specifying if the cluster is significantly differentially abundant.
#'
#' @name DAC-class
#' @rdname DAC-class
#' @exportClass DAC
DAC <- setClass("DAC",
    slots = c(sample.cond1 = "character",
        sample.cond2    = "character",
        cluster.size    = "numeric",
        use.percentages = "logical",
        method          = "character",
        method.paired   = "logical",
        method.adjust   = "character",
        th.fc           = "numeric",
        th.pvalue       = "numeric",
        results         = "data.frame"),
    validity = function(object){
        if (length(object@sample.cond1) == 0) {
            message("Error in DAC object: sample.cond1 length can not be equal to 0")
            return(FALSE)
        }
        if (length(object@sample.cond2) == 0) {
            message("Error in DAC object: sample.cond2 length can not be equal to 0")
            return(FALSE)
        }
        if (!any(object@method %in% c("t.test","wilcox.test"))) {
            message("Error in DAC object: method must match 't.test' or 'wilcox.test'")
            message(paste0("method founded is : ", object@method))
            return(FALSE)
        }
        if (!any(object@method.adjust %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"))) {
            message("Error in DAC object: method.adjust must match 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr', 'none')")
            message(paste0("method.adjust founded is : ", object@method.adjust))
            return(FALSE)
        }
        if (object@th.pvalue < 0 || object@th.pvalue > 1) {
            message("Error in DAC object: th.pvalue must be include into [0,1] interval")
            message(paste0("th.pvalue founded is : ", object@th.pvalue))
            return(FALSE)
        }
        if (object@th.fc == 0) {
            message("Error in DAC object: th.fc can not be equal to 0")
            message(paste0("th.fc founded is : ", object@th.fc))
            return(FALSE)
        }
        if (ncol(object@results) != 8) {
            print(colnames(object@results))
            message("Error in DAC object: results slot has not the expected number of columns (8)")
            message(paste0("Number of columns found are : ", ncol(object@results)))
            return(FALSE)
        }
        return(TRUE)
    }
)


#' @title Correlated Clusters (CC) class definition
#' 
#' @description 
#' The 'CC' object is a S4 object containing coefficient of correlation associated between each cluster to a phenotypic or functional variable.
#' Moreover, this object contains parameters and results used in the statistical analysis.
#' 
#' @details
#' A cluster is considered as a significant correlated cluster if its associated p-value and correlation threshold are below the specific thresholds 'th.pvalue' and 'th.correlation'.  
#' 
#' The 'print()' and 'show()' can be used to display a summary of this object.
#' Moreover all information about this object could be saved as a tab separated file using the 'export()' method.
#' This object is returned by the 'identifyCC()' function. 
#' 
#' @slot sample.names a character vector containing the samples used to compute correlated clusters
#' @slot variable a numeric vector containing the expression values of the associated variable
#' @slot cluster.size a numeric vector containing the number of cells for each cluster
#' @slot use.percentages a logical specifying if computation was performed on percentage of cell abundance
#' @slot method a character containing the name of the statistical test used to identify the CC
#' @slot method.adjust a character containing the name of the multiple correction method used (if any)
#' @slot th.correlation a numeric value specifying the correlation threshold (R)
#' @slot th.pvalue a numeric value specifying the p-value threshold
#' @slot results a data.frame containing for each cluster (first column): the coefficient of correlation R (second column) , the associated p-value (third column) and a logical (fourth column) specifying if the cluster is significantly correlated.
#' 
#' @name CC-class
#' @rdname CC-class
#' @exportClass CC
CC <- setClass("CC",
        slots = c(sample.names = "character",
        variable        = "numeric",
        cluster.size    = "numeric",
        use.percentages = "logical",
        method          = "character",
        method.adjust   = "character",
        th.correlation  = "numeric",
        th.pvalue       = "numeric",
        results         = "data.frame"),
    validity = function(object){
        if (length(object@variable) != length(object@sample.names)) {
            message(paste0("Error in CC object: variable length (",
            ,length(object@variable),
            ") is inconsistents with sample names length (",
            length(object@sample.names),")"))
            return(FALSE)
        }
        if (length(object@sample.names) == 0) {
            message("Error in CC object: sample.names length can not be equal to 0")
            return(FALSE)
        }
        if (!any(object@method %in% c("pearson", "kendall", "spearman"))) {
            message("Error in CC object: method must match 'pearson', 'kendall', 'spearman'")
            message(paste0("method founded is : ", object@method))
            return(FALSE)
        }
        if (!any(object@method.adjust %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"))) {
            message("Error in CC object: method.adjust must match 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr', 'none')")
            message(paste0("method.adjust founded is : ", object@method))
            return(FALSE)
        }
        if (object@th.pvalue < 0 || object@th.pvalue > 1) {
            message("Error in CC object: th.pvalue must be include into [0,1] interval")
            message(paste0("th.pvalue founded is : ", object@th.pvalue))
            return(FALSE)
        }
        if (object@th.correlation < 0 || object@th.correlation > 1) {
            message("Error in CC object: th.correlation must be include into [0,1] interval")
            message(paste0("th.correlation founded is : ", object@th.pvalue))
            return(FALSE)
        }
        if (!identical(colnames(object@results),c("cluster","correlation","pvalue","significant"))) {
            print(colnames(object@results))
            message("Error in CC object: result slot must have this colmuns : 'cluster','correlation','pvalue','significant'")
            message("Colmuns found are : ")
            message(paste(colnames(object@results)))
            return(FALSE)
        }
        return(TRUE)
    }
)


#' @title Abundance Profiles (AP) class definition
#' 
#' @description 
#' The 'AP' object is a S4 object containing the information related to the cluster classification based on theirs abundance profiles.
#' Moreover, this object contains parameters and results used in the statistical analysis.
#'     
#' @details 
#' Five methods are available to classify cellular clusters: 'hierarchical_k', 'hierarchical_h', 'kmeans', 'eigencell' and 'clique'.
#' Each method can parameterized using the 'method.parameter' parameter.
#'  
#' The 'print()' and 'show()' can be used to display a summary of this object. 
#' Moreover all information about this object could be saved as a tab separated file using the 'export()' method.
#' This object is returned by the 'classifyClusteringResults()' function. 
#'
#' @slot class.number a numeric value specifying the number of clusters
#' @slot cluster.size a numeric vector containing the number of cells for each cluster
#' @slot method a character specifying the method used to classify cluster
#' @slot method.parameter a named list of parameters used by the classification method
#' @slot classes a two column dataframe with the cluster in first column and corresponding class in the second column
#' 
#' @name AP-class
#' @rdname AP-class
#' @exportClass AP
AP <- setClass("AP",
    slots = c(class.number = "numeric",
        cluster.size     = "numeric",
        method           = "character",
        method.parameter = "numeric",
        classes          = "data.frame"),
    validity = function(object) {
        if (is.element(object@method, c("hierarchical_h", "eigencell", "clique")) &&
            (object@method.parameter > 1 || object@method.parameter < 0)) {
            message(paste0("Object AP, Error : with ", objec@tmethod, " method, the method.parameter must be include into [0,1] interval"))
            message(paste0("method.parameter founded is : ", object@th.pvalue))
            return(FALSE)
        }
        if (is.element(object@method, c("hierarchical_k", "kmeans")) &&
            (object@method.parameter != length(unique(object@classes$class)))) {
            message(paste0("Object AP, Error the number of class in the slot classes (", length(unique(object@classes$class)), ") is inconsistent the specified number of class ", method.parameter))
            return(FALSE)
        }
        return(TRUE)
    }
)


#' @title Unload 'flowSet' object from a 'Results' object
#'
#' @description 
#' This function unloads the 'flowSet' object in a 'Results' object.
#'
#' @param Results a 'Results' object
#'
#' @return The new 'Results' object
#'
#' @rdname unload.flowSet-methods
#' @export 
setGeneric("unload.flowSet",  function(Results)
    standardGeneric("unload.flowSet")
)

#' @rdname unload.flowSet-methods
#' @export
setMethod("unload.flowSet","Results",
    definition = function(Results){
        Results@flowset <- NULL
        gc()
        return(Results)
    }
)


#' @title Assignement of a context into a 'Results' object
#'
#' @description 
#' This function allows to assigns three kinds of contextual information to samples: biological condition, timepoints and individuals.
#' 
#' @details
#' The order of rownames (sample names) in this assignments dataframe will be the display order of the elements.
#' 
#' @param Results a Results or Results object
#' @param assignments a data.frame containing all samples names (in rownames) and columns providing contextual associations like "bc" (biological conditions), "tp" (biological conditions) and "ind" (individuals) 
#'
#' @return none
#'
#' @rdname assignContext-methods
#' @export
setGeneric("assignContext",  function(Results, assignments)
    standardGeneric("assignContext")
)

#' @rdname assignContext-methods
#' @export
setMethod("assignContext", "Results",
    definition = function(Results, assignments){
        if (is.null(assignments)) {
            stop("Error in assignContext: The 'assignments' parameter cannot be null")
        }
        if (!is.data.frame(assignments)) {
            stop("Error in assignContext: The 'assignments' parameter must be a dataframe")
        }
        if (!all(colnames(assignments) %in% c("bc", "tp", "ind"))) {
            stop("Error in assignContext: The 'assignments' dataframe can not contain other columns than : \"bc\" for biological conditions, \"tp\" for timepoints and \"ind\" for individuals")
        }

        Results@assignments <- assignments
        validObject(Results)

        message(paste0(ncol(assignments)," column(s) have been found: ", paste(colnames(assignments), collapse = ", ")))
        print(assignments)

        return(Results)
    }
)
 