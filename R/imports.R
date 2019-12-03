#' @title Importation of clustering results from text files
#'
#' @description
#' The 'importResultsFromTables()' function imports cell clustering results from two dataframes ('cluster.abundances' and 'cluster.phenotypes').
#' This function returns a 'Results' object.
#' 
#' @details
#' This function returns a 'Results' object without 'flowset', 'fcs.files', 'graph' and 'graph.layout' slots. 
#' 
#' The 'cluster.abundances' dataframe must be formatted with the cluster names in rownames as following:
#' 
#'     \tabular{ccccc}{
#'    x     \tab sample1 \tab sample2 \tab sample3 \tab ...\cr
#'    cluster1  \tab 749     \tab 5421    \tab   8424  \tab ...\cr
#'    cluster2  \tab 450     \tab 412     \tab   614   \tab ...\cr
#'    cluster3  \tab 288     \tab 782     \tab   478   \tab ...\cr
#'    ...       \tab ...     \tab ...     \tab   ...   \tab ...
#' }
#' 
#' The 'cluster.phenotypes' dataframe must be formatted as following:
#' 
#'    \tabular{cccccc}{
#'    sample     \tab cluster     \tab marker1     \tab marker2  \tab marker3 \tab ...\cr
#'    sample1  \tab cluster1 \tab 0.212    \tab 0.445   \tab 1.756   \tab ...\cr
#'    sample1  \tab cluster3 \tab 0.434    \tab 0.785   \tab 0.654   \tab ...\cr
#'    sample2  \tab cluster1 \tab 0.574    \tab 2.641   \tab 3.854   \tab ...\cr
#'    sample2  \tab cluster2 \tab 1.454    \tab 1.755   \tab -0.568  \tab ...\cr
#'    sample2  \tab cluster3 \tab 1.445    \tab 1.875   \tab -0.786  \tab ...\cr
#'    sample3  \tab cluster1 \tab 0.157    \tab 2.121   \tab 1.648   \tab ...\cr
#'    sample3  \tab cluster2 \tab 1.415    \tab 1.963   \tab 0.786   \tab ...\cr
#'    sample3  \tab cluster3 \tab 1.275    \tab 1.427   \tab 0.754   \tab ...\cr
#'    ...      \tab ...      \tab ...      \tab ...     \tab ...     \tab ...
#' }
#' 
#' @param cluster.abundances a dataframe of cells abundances with clusters in row and samples in column 
#' @param cluster.phenotypes a dataframe containing median marker expression values for each cluster of each sample. In additions of markers, the two first columns are are dedicated to "cluster" and "sample" identifiers.
#' @param clustering.markers a character vector specifying markers that have been used during the clustering procedure (if NULL, all markers will be considered as clustering markers)
#' @param th.min_cells a numeric specifying the minimum number of cell in a cluster of a sample to take its phenotype in account
#' @param assignments a data.frame containing all samples names (in rownames) and columns providing contextual associations like "bc" (biological conditions), "tp" (biological conditions) and "ind" (individuals) 
#'
#' @return a S4 object of class 'Results'
#' 
#' @export 
importResultsFromTables <- function(cluster.abundances,
                                    cluster.phenotypes,
                                    clustering.markers = NULL,
                                    th.min_cells      = 0,
                                    assignments       = NULL) {
    
    message("[START] - importing cell clustering results")
    
    cluster.abundances <- removeFactors(cluster.abundances)
    cluster.phenotypes <- removeFactors(cluster.phenotypes)
    
    if (!is.data.frame(cluster.abundances)) {
        stop("Error in importResultsFromTables(: The 'cluster.abundances' parameter must be a dataframe")
    }
    
    if (!all(apply(cluster.abundances, 2, is.numeric))) {
        stop("Error in importResultsFromTables(: The 'cluster.abundances' parameter must be numeric")
    }
    
    if (!all(apply(cluster.abundances, 2, function(x) {(x >= 0) }))) {
        
        stop("Error in importResultsFromTables(: The 'cluster.abundances' parameter must contain only positive values")
    }
    
    if (!is.data.frame(cluster.phenotypes)) {
        stop("Error in importResultsFromTables(: The 'cluster.phenotypes' parameter must be a dataframe")
    }
    
    if (!is.character(cluster.phenotypes[, 1])) {
        stop("Error in importResultsFromTables(: The 'cluster.phenotypes' parameter must be contain a dataframe with a character in the first column")
    }
    if (!is.character(cluster.phenotypes[, 2])) {
        stop("Error in importResultsFromTables(: The 'cluster.phenotypes' parameter must be contain a dataframe with a character in the second column")
    }
    if (!all(apply(cluster.phenotypes[, -c(1, 2)], 2, is.numeric))) {
        stop("Error in importResultsFromTables(: The 'cluster.phenotypes' parameter must be contain numeric values (with exception for the two first columns")
    }
    
    if (th.min_cells < 0) {
        stop("Error in importResultsFromTables(: The 'th.min_cells' parameter must be a stricly positive integer")
    }
    
    message("\tloading cluster abundances...")
    colnames(cluster.phenotypes)[1] <- "sample"
    message("\tloading cluster phenotypes...")
    colnames(cluster.phenotypes)[2] <- "cluster"
    
    for (cluster in rownames(cluster.abundances)) {
        for (sample in colnames(cluster.abundances)) {
            if (cluster.abundances[cluster, sample] < th.min_cells) {
                cluster.phenotypes[(cluster.phenotypes$cluster == cluster & cluster.phenotypes$sample == sample), 3:ncol(cluster.phenotypes)] <- rep(NA, ncol(cluster.phenotypes) - 2)
            }
        }
    }
    
    message("\tcomputing minimal and maximal bounds...")
    min.bounds <- apply(cluster.phenotypes[, -c(1, 2)], 2, min, na.rm = TRUE)
    max.bounds <- apply(cluster.phenotypes[, -c(1, 2)], 2, max, na.rm = TRUE)
    bounds <- as.data.frame(rbind("0.00" = min.bounds, "1.00" = max.bounds))
    
    marker.names <- colnames(cluster.phenotypes)[3:length(cluster.phenotypes)]
    
    if (is.null(clustering.markers)) {
        clustering.markers <- marker.names
    }
    
    res <- methods::new("Results",
                        cluster.phenotypes = cluster.phenotypes,
                        cluster.abundances = cluster.abundances,
                        cluster.names      = rownames(cluster.abundances),
                        sample.names       = colnames(cluster.abundances),
                        marker.names       = marker.names,
                        cluster.number     = nrow(cluster.abundances),
                        th.min_cells       = th.min_cells,
                        bounds             = bounds,
                        clustering.markers = clustering.markers)

    if (!is.null(assignments)) {
        res <- assignContext(res, assignments)
    }
        
    print(res)
    message("[END] - importing cell clustering results")
    
    return(res)
    
}


#' @title Importation of clustering results from a folder containing FCS files
#'
#' @description 
#' The 'importResultsFromFCS()' function imports cell clustering results from FCS files contained in a specified folder.
#' This function imports the cluster phenotype matrix and count matrix.
#' This function apply an hyperbolic sine transformation to imported FCS data (unless 'trans' = 'logicle' or 'none') and compute the marker range quantiles.
#' 
#' @details 
#' This function returns a 'Results' object including 'flowset', 'fcs.files' slots but not 'graph' and 'graph.layout' slots. 
#' The computation of marker range quantiles can be approximated using 'quantile.approximation' parameter which is more efficient in term of loading time and memory usage.
#' The name of each file will be used as the name of the sample. All FCS files corresponding to all the samples must be contained in a single folder. Importantly, FCS files must contain a channel called "cluster" allowing to identify which cluster is associated to each cell.
#'
#' @param path a character specify the path of the FCS files
#' @param exclude.markers a character vector of markers to exclude (case insensitive)
#' @param clustering.markers a character vector specifying markers that have been used during the clustering procedure (if NULL, all markers will be considered as clustering markers)
#' @param probs a vector of probabilities with 2 values in [0,1] to compute marker range quantiles. First is the lower bound and second is the upper bound.
#' @param trans a character specifying what transformation ("arcsinh", "logicle", or "none") will be applied on the cluster expression matrix ("arcsinh" by default)
#' @param quantile.approximation a logical specifying if marker range quantiles are computed using all cells (FALSE), or is the means of the quantile of each samples (TRUE)
#' @param th.min_cells a numeric specifying the minimum number of cell in a cluster of a sample to take its phenotype in account
#' @param assignments a data.frame containing all samples names (in rownames) and columns providing contextual associations like "bc" (biological conditions), "tp" (biological conditions) and "ind" (individuals)
#'
#' @return a S4 object of class 'Results'
#'
#' @export 
#' 
#' @importFrom flowCore sampleNames
#' @importFrom gtools mixedsort
importResultsFromFCS <- function(path,
                                 exclude.markers        = c("cell_length", "FileNum", "density", "time"),
                                 clustering.markers     = NULL,
                                 probs                  = c(0.05, 0.95),
                                 trans                  = "arcsinh",
                                 quantile.approximation = FALSE,
                                 th.min_cells           = 0,
                                 assignments            = NULL){
    
    message("[START] - importing clustering results from FCS file")
    path <- normalizePath(path, "/", mustWork = TRUE)
    
    message(paste0(basename(path), "\n"))
    
    if (typeof(exclude.markers) != "character") {
        stop("Error in importFromFCS: The 'exclude.markers' parameter must be a character vector")
    }
    
    if (length(probs) != 2) {
        stop("Error in importFromFCS: The 'probs' parameter must only 2 numeric values")
    } else if (probs[1] > probs[2]) {
        stop("Error in importFromFCS: The 'probs' parameter must contain a first value greather than the second value")
    } else if (probs < 0 || probs > 1) {
        stop("Error in importFromFCS: The 'probs' parameter must contain values included in the domain: [0;1]")
    }
    
    if (!is.character(trans)) { stop("Error in importFromFCS: The 'trans' parameter must be a character") }
    if (!is.logical(quantile.approximation)) { stop("Error in importFromFCS: The 'quantile.approximation' parameter must be a logical") }
    
    if (th.min_cells < 0) { stop("Error in importFromFCS: The 'th.min_cells' parameter must be a stricly positive integer")  }
    
    fcs.files  <- dir(path, full.names = TRUE, pattern = ".fcs$")
    list       <- load.flowSet(fcs.files = fcs.files, exclude.markers = exclude.markers, trans = trans, pattern = ".fcs")
    flowset    <- list$flowset
    dictionary <- list$dictionary
    
    message("\tcompute quantiles bounds...")
    
    if (quantile.approximation) {
        quantiles <- computeQuantile.approximation(flowset, probs)
    } else {
        quantiles <- computeQuantile(flowset, probs)
    }
    gc()
    message("\textract results from FCS files...")
    
    samples  <- flowCore::sampleNames(flowset)
    markers  <- setdiff(flowset@colnames, "cluster")
    clusters <- c()
    
    for (sample in samples) { 
        frame    <- flowset[[sample]]@exprs
        clusters <- unique(c(clusters, frame[, "cluster"]))
    }
    
    clusters <- gtools::mixedsort(clusters)
    
    cluster.abundances <- matrix(0, nrow = length(clusters), ncol = length(samples), dimnames = list(clusters, samples))
    nrow               <- length(clusters) * length(samples)
    ncol               <- length(markers) + 2
    cluster.phenotypes <- matrix(NA, nrow = nrow, ncol = ncol, dimnames = list(NULL, c("sample", "cluster", markers)))
    cluster.abundances <- as.data.frame(cluster.abundances, stringsAsFactors = FALSE)
    cluster.phenotypes <- as.data.frame(cluster.phenotypes, stringsAsFactors = FALSE)
    
    cluster.phenotypes$sample  <- rep(samples, each = length(clusters))
    cluster.phenotypes$cluster <- rep(clusters, times = length(samples))
    
    for (sample in samples) {
        frame    <- flowset[[sample]]@exprs
        for (cluster in clusters) {
            frame.cluster <- frame[frame[, "cluster"] == cluster, colnames(frame) != "cluster", drop = FALSE]
            cluster.abundances[cluster, sample] <- nrow(frame.cluster)
            row <- c(sample, cluster, apply(frame.cluster, 2, stats::median, na.rm = TRUE))
            cluster.phenotypes[cluster.phenotypes$sample == sample & cluster.phenotypes$cluster == cluster, ] <- row
        }
    }
    
    cluster.phenotypes$sample  <- as.vector(cluster.phenotypes$sample)
    cluster.phenotypes$cluster <- as.vector(cluster.phenotypes$cluster)
    
    for (i in 3:ncol(cluster.phenotypes)) {
        cluster.phenotypes[, i] <- as.numeric(cluster.phenotypes[, i])
    }
    
    if (is.null(clustering.markers)) {
        clustering.markers <- markers
    }
    
    res <- methods::new("Results",
                        cluster.phenotypes = cluster.phenotypes,
                        trans              = trans,
                        cluster.abundances = cluster.abundances,
                        cluster.names      = rownames(cluster.abundances),
                        sample.names       = flowCore::sampleNames(flowset),
                        marker.names       = markers,
                        clustering.markers = clustering.markers,
                        cluster.number     = nrow(cluster.abundances),
                        flowset            = flowset,
                        fcs.files          = fcs.files,
                        bounds             = quantiles,
                        th.min_cells       = th.min_cells,
                        graph.layout       = NULL,
                        graph              = NULL)
    
    if (!is.null(assignments)) {
        res <- assignContext(res, assignments)
    }
    
    print(res)
    message("[END] - importing clustering results from FCS file")
    
    return(res)
    
}


#' @title Importation of clustering results from a folder containing ACS files
#'
#' @description 
#' The 'importResultsFromCLR.ACS()' function imports cell clustering results from CLR files contained in a specified folder.
#' CLR files must be in ACS format.
#' This function imports the cluster phenotype and abundance.
#' 
#' @details 
#' This function returns a 'Results' object including 'flowset', 'fcs.files' slots but not 'graph' and 'graph.layout' slots. 
#' The computation of marker range quantiles can be approximated using 'quantile.approximation' parameter which is more efficient in term of loading time and memory usage.
#' The name of each file will be used as the name of the sample. All FCS files corresponding to all the samples must be contained in a single folder.
#'
#' For each cell, the selected cluster will be one having the greater probability higher than the 'prob.th' parameter.
#'
#' @param path a character specify the path of the ACS CLR files
#' @param prob.th a numeric specifying the probability threshold to use for considering considering the cluster of each each cells
#' @param exclude.markers a character vector of markers to exclude (case insensitive)
#' @param probs a vector of probabilities with 2 values in [0,1] to compute marker range quantiles. First is the lower bound and second is the upper bound.
#' @param trans a character specifying what transformation ("arcsinh", "logicle", or "none") will be applied on the cluster expression matrix ("arcsinh" by default)
#' @param quantile.approximation a logical specifying if marker range quantiles are computed using all cells (FALSE), or is the means of the quantile of each samples (TRUE)
#' @param th.min_cells a numeric specifying the minimum number of cell in a cluster of a sample to take its phenotype in account
#' @param assignments a data.frame containing all samples names (in rownames) and columns providing contextual associations like "bc" (biological conditions), "tp" (biological conditions) and "ind" (individuals)
#'
#' @return a S4 object of class 'Results'
#'
#' @export 
importResultsFromCLR.ACS <- function(path,
                                     prob.th                = 0.80,
                                     th.min_cells           = 0,
    								 probs                  = c(0.05, 0.95),
                                     trans                  = "arcsinh",
                                     quantile.approximation = FALSE,
    								 exclude.markers        = NULL,
                                     assignments            = NULL){
	
	message("[START] - importing clustering results from a CLR ACS file")
    
    message(paste0(basename(path), "\n"))
      
    if (th.min_cells < 0) { stop("Error in importResultsFromCLR.ACS: The 'th.min_cells' parameter must be a stricly positive integer")  }
    
    path <- normalizePath(path, "/", mustWork = TRUE)
    
    files.asc  <- dir(path, full.names = TRUE, pattern = ".acs$")
    csv.files  <- c()
    fcs.files  <- c()
   
    
    for(file in files.asc){
        temp <- paste0(dirname(file),"/",".extract-",basename(file))

        if(!file.exists(temp)){
            dir.create(temp)
        }
        utils::unzip(file,exdir=temp)
        csv.file  <- dir(temp, full.names = TRUE, pattern = ".csv$")
        fcs.file  <- dir(temp, full.names = TRUE, pattern = ".fcs$")
        csv.files <- c(csv.files,csv.file)
        fcs.files <- c(fcs.files,fcs.file)
        
    }
    
	# count
	cluster.abundances <- c()
	samples            <- c()
	clusters           <- NA
	
	associated.clusters <- vector("list", length(csv.files))
	cluster.abundances  <- vector("list", length(csv.files))
	
	for(file in csv.files){
	    
	    sample                          <- gsub(".csv","",basename(file))
	    samples                         <- c(samples,sample)
		all.cells                       <- utils::read.delim(file,sep=",",quote="")

		if (any(is.na(clusters))){
		    clusters <- colnames(all.cells)
		}
		
		if(!identical(colnames(all.cells), clusters)){
		    stop("Error in importResultsFromCLR.ACS: the header of the csv files must be identical")
		}
		
		cells.associated.clusters <- apply(all.cells, 1, FUN = function (x){
		   cluster <- NA
		   if (any(!is.na(x))) {
		       x <- x[which.max(x)]
		       if (x > prob.th){
		        cluster <- ifelse((length(x)>0) & (x > prob.th), names(x), NA)
		       }
		   }
		})
		
		nb_cells.by.cluster <- as.data.frame(table(cells.associated.clusters))
		sum                 <- stats::setNames(nb_cells.by.cluster$Freq, nb_cells.by.cluster$cells.associated.clusters)
		
		represented.clusters <- names(sum)

		missing              <- setdiff(clusters,represented.clusters)
		
		sum                  <- c(sum, rep(0,length(missing)))
		names(sum)           <- c(represented.clusters,missing)
		
		associated.clusters[[sample]] <- cells.associated.clusters
		cluster.abundances[[sample]]  <- sum
		
	}
	
	cluster.abundances <- as.data.frame(do.call(cbind, cluster.abundances))

	list       <- load.flowSet(fcs.files = fcs.files, exclude.markers = exclude.markers, trans = trans, pattern = ".fcs")
	flowset    <- list$flowset
    dictionary <- list$dictionary

	if (quantile.approximation) {
        quantiles <- computeQuantile.approximation(flowset, probs)
    } else {
        quantiles <- computeQuantile(flowset, probs)
    }
    
	markers  <- colnames(quantiles)
	samples  <- flowCore::sampleNames(flowset)
	
	nrow               <- length(clusters) * length(samples)
    ncol               <- length(markers) + 2
    cluster.phenotypes <- matrix(NA, nrow = nrow, ncol = ncol, dimnames = list(NULL, c("sample", "cluster", markers)))
    cluster.phenotypes <- as.data.frame(cluster.phenotypes, stringsAsFactors = FALSE)
    
    cluster.phenotypes$sample  <- rep(samples, each = length(clusters))
    cluster.phenotypes$cluster <- rep(clusters, times = length(samples))
 
    for (sample in samples) {
    	    
        for (cluster in clusters){
    	        
    	   flowframe <- flowset[[sample]]@exprs
            
           marker.median.exprs <- apply(flowframe[associated.clusters[[sample]] %in% cluster,], 2, stats::median)
                
           cluster.phenotypes[cluster.phenotypes$sample == sample & cluster.phenotypes$cluster == cluster, 3:ncol(cluster.phenotypes)] <- marker.median.exprs

        }
    	   
    }
	
	res <- methods::new("Results",
                        cluster.phenotypes = cluster.phenotypes,
                        trans              = trans,
                        cluster.abundances = cluster.abundances,
                        cluster.names      = rownames(cluster.abundances),
                        sample.names       = colnames(cluster.abundances),
                        marker.names       = markers,
                        clustering.markers = markers,
                        cluster.number     = nrow(cluster.abundances),
                        flowset            = flowset,
                        fcs.files          = fcs.files,
                        bounds             = quantiles,
                        th.min_cells       = th.min_cells,
                        graph.layout       = NULL,
                        graph              = NULL)
	 
    if (!is.null(assignments)) {
        res <- assignContext(res, assignments)
    }
    message("[END] - importing clustering results from CLR ASC files")
    
    unlink(paste0(dirname(files.asc),"/",".extract-",basename(files.asc)),recursive = TRUE)

    print(res)
    return(res)							 
}

#' @title Importation of clustering results from a folder containing CVS files
#'
#' @description 
#' The 'importResultsFromCLR.CSV()' function imports cell clustering results from CLR files contained in a specified folder.
#' CLR files must be in CSV format ().
#' This function imports only the cluster abundance.
#' 
#' @details 
#' This function returns a 'Results' object with only cluster abundance information (no phenotype). 
#' The name of each file will be used as the name of the sample.
#' 
#' For each cell, the selected cluster will be one having the greater probability higher than the 'prob.th' parameter.
#'
#' @param path a character specify the path of the CVS CLR files
#' @param prob.th a numeric specifying the probability threshold to use for considering considering the cluster of each each cells
#' @param th.min_cells a numeric specifying the minimum number of cell in a cluster of a sample to take its phenotype in account
#' @param assignments a data.frame containing all samples names (in rownames) and columns providing contextual associations like "bc" (biological conditions), "tp" (biological conditions) and "ind" (individuals)
#'
#' @return a S4 object of class 'Results'
#'
#' @export 
importResultsFromCLR.CSV <- function(path,
                                     prob.th                = 0.80,
                                     th.min_cells           = 0,
                                     assignments            = NULL){
    
    message("[START] - importing clustering results from a CLR CSV file")
    
    message(paste0(basename(path), "\n"))
    
    if (th.min_cells < 0) { stop("Error in importResultsFromCLR CVS: The 'th.min_cells' parameter must be a stricly positive integer")  }
    
    file <- normalizePath(path, "/", mustWork = TRUE)
    
    # count
    cluster.abundances <- c()
    samples  <- c()
    clusters <- c()

    # count
    files.csv  <- dir(path, full.names = TRUE, pattern = ".csv$")
    cluster.abundances <- c()
    for(file in files.csv){
        sample  <- gsub(".csv","",basename(file))
        samples <- c(samples, sample)
        abundances <- utils::read.delim(file,sep=",",quote="")
        apply(abundances, 1, FUN = function(x){
            if (!any(is.na(x))){
                idx.max    <- which.max(x)
                x[]        <- 0
                if(x[idx.max] >= prob.th){
                    x[idx.max] <- 1
                }
            }
            return(x)
        })
        sum                <- apply(abundances, 2, sum, na.rm = TRUE)
        cluster.abundances <- cbind(cluster.abundances, sum)
    }

    cluster.abundances           <- as.data.frame(cluster.abundances)
    colnames(cluster.abundances) <- samples
    cluster.phenotypes           <- data.frame()
    markers                      <- character(0)
    flowset                      <- NULL
    quantiles                    <- data.frame()

    res <- methods::new("Results",
                        cluster.phenotypes = cluster.phenotypes,
                        trans              = NA,
                        cluster.abundances = cluster.abundances,
                        cluster.names      = rownames(cluster.abundances),
                        sample.names       = colnames(cluster.abundances),
                        marker.names       = markers,
                        clustering.markers = markers,
                        cluster.number     = nrow(cluster.abundances),
                        flowset            = flowset,
                        fcs.files          = character(0),
                        bounds             = quantiles,
                        th.min_cells       = th.min_cells,
                        graph.layout       = NULL,
                        graph              = NULL)
    
    if (!is.null(assignments)) {
        res <- assignContext(res, assignments)
    }
    message("[END] - importing clustering results from CLR CSV files")
    print(res)
    return(res)							 
}





#' @title Importation of clustering results generated by SPADE
#'
#' @description 
#' The 'importResultsFromSPADE()' function imports SPADE cell clustering results from a specified path.
#' This function imports the cluster phenotype matrix and count matrix as well as the SPADE tree.
#' This function apply an hyperbolic sine transformation to imported FCS data (unless 'trans' = log or none) and compute the marker range quantiles.
#' 
#' @details
#' This function returns a 'Results' object including 'flowset', 'fcs.files', 'graph' and 'graph.layout' slots. 
#' The computation of marker range quantiles can be approximated using 'quantile.approximation' parameter which is more efficient in term of loading time and memory usage.
#'  
#' @param path a character specify the path of SPADE results folder
#' @param exclude.markers a character vector of markers to exclude (case insensitive)
#' @param probs a vector of probabilities with 2 values in [0,1] to compute marker range quantiles. First is the lower bound and second is the upper bound.
#' @param trans a character specifying what transformation ("arcsinh", "logicle", or "none") will be applied on the cluster expression matrix ("arcsinh" by default)
#' @param quantile.approximation a logical specifying if marker range quantiles are computed using all cells (FALSE), or is the means of the quantile of each samples (TRUE)
#' @param th.min_cells a numeric specifying the minimum number of cell in a cluster of a sample to take its phenotype in account
#' @param load.phenotype a logical specifying if the phenotype matrix and fcs file will be loaded
#' @param assignments a data.frame containing all samples names (in rownames) and columns providing contextual associations like "bc" (biological conditions), "tp" (biological conditions) and "ind" (individuals)
#'
#' @return a S4 object of class 'Results'
#'
#' @export 
#' 
#' @import igraph
importResultsFromSPADE <- function (path, exclude.markers = c("cell_length", "FileNum", 
    "density", "time"), probs = c(0.05, 0.95), trans = "arcsinh", 
    quantile.approximation = FALSE, th.min_cells = 0, load.phenotype = TRUE, 
    assignments = NULL) 
{
    message("[START] - importing SPADE clustering results")
    path <- normalizePath(path, "/", mustWork = TRUE)
    message(paste0(basename(path), "\n"))
    if (typeof(exclude.markers) != "character") {
        stop("Error in importResultsFromSPADE: The 'exclude.markers' parameter must be a character vector")
    }
    if (length(probs) != 2) {
        stop("Error in importResultsFromSPADE: The 'probs' parameter must only 2 numeric values")
    }
    else if (probs[1] > probs[2]) {
        stop("Error in importResultsFromSPADE: The 'probs' parameter must contain a first value greather than the second value")
    }
    else if (probs < 0 || probs > 1) {
        stop("Error in importResultsFromSPADE: The 'probs' parameter must contain values included in the domain: [0;1]")
    }
    if (is.na(match(trans, c("arcsinh", "logicle", "none")))) {
        stop("Error in importResultsFromSPADE: The 'trans' parameter must be arcsinh, log or none")
    }
    if (!is.logical(quantile.approximation)) {
        stop("Error in importResultsFromSPADE: The 'quantile.approximation' parameter must be a logical")
    }
    if (th.min_cells < 0) {
        stop("Error in importResultsFromSPADE: The 'th.min_cells' parameter must be a stricly positive integer")
    }
    if (load.phenotype) {
        fcs.files <- dir(path, full.names = TRUE, pattern = ".fcs.density.fcs.cluster.fcs$")
        list <- load.flowSet(fcs.files = fcs.files, exclude.markers = exclude.markers, 
            trans = trans)
        flowset <- list$flowset
        dictionary <- list$dictionary
        message("\tcompute quantiles bounds...")
        if (quantile.approximation) {
            quantiles <- computeQuantile.approximation(flowset, 
                probs)
        }
        else {
            quantiles <- computeQuantile(flowset, probs)
        }
        gc()
    }
    else {
        fcs.files <- character(0)
        flowset <- NULL
        quantiles <- data.frame()
    }
    message("\treading SPADE results...")
    files <- dir(paste(path, "/tables/bySample/", sep = ""), 
        full.names = TRUE)
    if (length(files) == 0) 
        stop("Error when importing cell cluster abundances. Please check that subfolder \"./tables/bySample/\" is well existing")
    cluster.phenotypes <- data.frame(stringsAsFactors = FALSE)
    cluster.abundances <- data.frame()
    sample.names <- c()
    for (file in files) {
        name <- gsub(".fcs.density.fcs.cluster.fcs.anno.Rsave_table.csv$", 
            "", basename(file))
        sample.names <- c(sample.names, name)
        SPADE.matrix <- utils::read.table(file, sep = ",", header = TRUE, 
            stringsAsFactors = FALSE, check.names = FALSE)
        SPADE.matrix[, "ID"] <- as.character(SPADE.matrix[, "ID"])
        cluster.abundances.sample <- SPADE.matrix[, "count"]
        if (nrow(cluster.abundances)) {
            cluster.abundances <- cbind(cluster.abundances, cluster.abundances.sample)
            samples.headers <- append(samples.headers, name)
        }
        else {
            cluster.abundances <- data.frame(row.names = SPADE.matrix[, 
                "ID"], cluster.abundances.sample)
            samples.headers <- name
        }
        if (load.phenotype) {
            cluster.phenotypes.sample <- SPADE.matrix[, grep("count|percenttotal", 
                colnames(SPADE.matrix), invert = TRUE)]
            cluster.phenotypes.sample[cluster.abundances.sample < 
                th.min_cells, 2:ncol(cluster.phenotypes.sample)] <- rep(NA, 
                ncol(cluster.phenotypes.sample) - 1)
            cluster.phenotypes.sample <- cbind(name = rep(name, 
                nrow(cluster.phenotypes.sample)), cluster.phenotypes.sample)
            cluster.phenotypes <- rbind(cluster.phenotypes, cluster.phenotypes.sample)
        }
    }
    colnames(cluster.abundances) <- samples.headers
    if (load.phenotype) {
        cluster.phenotypes.header <- colnames(cluster.phenotypes)
        cluster.phenotypes.header[1] <- "sample"
        cluster.phenotypes.header[2] <- "cluster"
        colnames(cluster.phenotypes) <- cluster.phenotypes.header
        cluster.phenotypes <- filter.medians(cluster.phenotypes, 
            trans)
        cluster.phenotypes.header <- colnames(cluster.phenotypes)
        clustering.markers.index <- grep("_clust", cluster.phenotypes.header)
        cluster.phenotypes.header <- gsub("_clust", "", cluster.phenotypes.header)
        colnames(cluster.phenotypes) <- rename.markers(cluster.phenotypes.header, 
            dictionary = dictionary)
        clustering.markers <- colnames(cluster.phenotypes)[clustering.markers.index]
        if (!is.null(exclude.markers)) {
            cluster.phenotypes <- exclude.markers(cluster.phenotypes, 
                exclude.markers)
            clustering.markers <- setdiff(clustering.markers, 
                exclude.markers)
        }
        graph <- igraph::read.graph(paste(path, "/mst.gml", sep = ""), 
            format = "gml")
        graph.layout <- as.matrix(utils::read.table(paste0(path, 
            "/layout.table"), sep = " ", quote = "", stringsAsFactors = FALSE))
        markers.names <- colnames(cluster.phenotypes[, -c(1, 
            2)])
    }
    else {
        clustering.markers <- character(0)
        markers.names <- character(0)
        graph.layout <- NULL
        graph <- NULL
    }
    res <- methods::new("Results", cluster.phenotypes = cluster.phenotypes, 
        trans = trans, cluster.abundances = cluster.abundances, 
        cluster.names = rownames(cluster.abundances), sample.names = sample.names, 
        marker.names = markers.names, clustering.markers = clustering.markers, 
        cluster.number = nrow(cluster.abundances), flowset = flowset, 
        fcs.files = fcs.files, bounds = quantiles, th.min_cells = th.min_cells, 
        graph.layout = graph.layout, graph = graph)
    if (!is.null(assignments)) {
        res <- assignContext(res, assignments)
    }
    message("[END] - importing SPADE clustering results")
    return(res)
}

# @title Internal - Renaming of cell markers
# 
# @description 
# This function is used internally to rename the cell markers based on a dictionary.
#
# @details 
# Dictionary is a data.frame used to rename the marker names. 
# The first column must correspond to the original marker names, the second column must correspond to the new marker names. 
# In case of a redondance in the second column, the marker will be renamed in this way oldmarkernames-().
#
# @param header a character vector containing the original marker names
# @param dictionary a character vector containing a correspondence between the original and the new marker names
#
# @return a character vector containing the renamed marker names
rename.markers <- function(header, dictionary) {
    
    header <- make.names(header)
    
    dictionary[, 1] <- as.vector(dictionary[, 1])
    dictionary[, 2] <- as.vector(dictionary[, 2])
	
	occurences <- table(dictionary[, 2])
    occurences <- occurences[occurences > 1]
    redondant.names <- names(occurences)
    
    for (i in 1:nrow(dictionary)) {
        if (any(dictionary[i, 2] %in% redondant.names)) {
            temp             <- gsub("X\\.(.*)\\.Di(_clust$|$)","(\\1)Di\\2", dictionary[i, 1])
            dictionary[i, 2] <- paste0(dictionary[i, 2], "-", temp)
        }
        header[which(header == dictionary[i, 1])[1]] <- dictionary[i, 2]
    }
	
	header[is.na(header)] <- dictionary[, 1][is.na(header)]
   
    return(header)
}

# @title Internal - Removing of cell markers to exclude from a matrix
#
# @description 
# This function is used internally to remove one or several cell markers.
# 
# @details 
# If the data parameter is a dataframe the colnames.FCS parameter is ignored but if the data parameter is a flowset, the colnames.FCS parameter is required.
# 
# @param data a numeric matrix or flowset
# @param exclude a character vector containing the cell markers to be excluded (case intensive)
# @param colnames.FCS a character vector containing column names if data is a FCS flowset
# 
# @return a numeric matrix without the cell markers to exclude
exclude.markers <- function(data, exclude, colnames.FCS = NULL){
    
    if (!is.null(colnames.FCS)) {
        column <- colnames.FCS
    } else {
        column <- colnames(data) 
    }
    
    exclude.flags <- toupper(exclude) %in% toupper(column)
    
    if (any(!(exclude.flags))) {
        warning(paste0("Unknown marker to exclude: ", paste(exclude[!exclude.flags], collapse = ", ")))
    }
    
    data    <- data[ , -which(toupper(column) %in% toupper(exclude))]
    
    return(data)
}

# @title Internal - Filtering of medians from a SPADE result matrix
#
# @description 
# This function is used internally to remove raw or transform medians from a SPADE result matrix. CVS medians are always removed.
# 
# @param data a SPADE matrix
# @param trans a character specifying what transformation ("arcsinh", "logicle", or "none") will be applied on the cluster expression matrix ("arcsinh" by default)# 
# @return a numeric matrix without the cell markers to exclude
filter.medians <- function (data, trans = "arcsinh") {
    if (trans=="none") {
        exclude <- "^medians|^cvs"
    }
    else {
        exclude <- "^raw_medians|^cvs"
    }
    data <- data[, grep(exclude, colnames(data), invert = TRUE, 
        ignore.case = TRUE)]
    colnames(data) <- gsub("^medians|^cvs|^raw_medians", "", 
        colnames(data))
    return(data)
}

# @title Internal - Computation of quantile with FCS flowset marker by marker 
#
# @description 
# This function is used internally to compute the marker range quantiles.
# 
# @details 
# This function performs the exact calculation of quantiles with all cells but needs more resources (time and memory usage) than 'computeQuantile.approximation'.
# 
# @param flowset a flowCore flowset
# @param probs a numeric vector of 2 values specifying the quantiles to compute
# 
# @return a numeric matrix containing the quantiles of each marker
computeQuantile <- function(flowset, probs = c(0.05,0.95)){

    bounds  <- data.frame()
    markers <- flowset@colnames
    markers <- setdiff(markers, "cluster")
    
    for (marker in markers) {
        temp <- c()
        for (sample in 1:length(flowset)) {    
            frame <- flowset[[sample]]@exprs
            temp  <- c(temp, frame[, marker])
        }
        if (nrow(bounds) > 0) {
            bounds <- cbind(bounds, stats::quantile(temp, probs = probs))
        }else{
            bounds <- as.data.frame(stats::quantile(temp, probs = probs))
        }
    }
    
    colnames(bounds) <- markers
    rownames(bounds) <- probs
    
    return(bounds)
}


# @title Internal - Computation of quantiles with FCS flowset sample by sample
#
# @description 
# This function is used internally to compute bound using the 2 selected percentiles from each sample. 
# For the first percentile, the mean of this percentile from all samples is compute (same for the second percentiles). 
# This approach avoid to compute percentile based the whole dataset and then seed up computation.
# 
# @details 
# This function performs an approximate calculation of quantiles using less memory than computeQuantile.
# 
# @param flowset a flowCore flowset 
# @param probs a numeric vector of 2 values specifying the quantiles to compute
# 
# @return a numeric matrix containing the quantiles of each marker
# 
#' @importFrom flowCore fsApply
computeQuantile.approximation <- function(flowset,probs = c(0.05,0.95)){
    
    bounds.by.sample <- flowCore::fsApply(flowset[, flowset@colnames != "cluster"], flowCore::each_col, stats::quantile, probs = probs)
    
    lower.bounds <- bounds.by.sample[seq(from = 1, to = nrow(bounds.by.sample), by = 2), ]
    upper.bounds <- bounds.by.sample[seq(from = 2, to = nrow(bounds.by.sample), by = 2), ]
    
    lower.bounds <- apply(lower.bounds, 2, mean)
    upper.bounds <- apply(upper.bounds, 2, mean)
    
    bounds <- data.frame(row.names = probs, rbind(lower.bounds, upper.bounds), check.names = FALSE)
    
    return(bounds)
    
}


#' @title Loading of FCS files object into a 'Results' object
#'
#' @description 
#' This function loads the FCS files to the 'flowset' slot of the 'Results' object.
#' 
#' @details
#' If a 'Results' object is provided, others parameters ('fcs.files', 'exclude.markers' and 'trans') will be ignored (they will be retrieved from the 'Results' object). 
#' 
#' @param Results a Results object (with 'fcs.files' slot not null) (optional)
#' @param fcs.files a character vector containing the absolute path of the original FCS files
#' @param exclude.markers a character vector of markers to exclude (case insensitive)
#' @param trans a character specifying what transformation ("arcsinh", "logicle", or "none") will be applied on the cluster expression matrix ("arcsinh" by default)
#' @param pattern a character specifying the pattern of the FCS file 
#' 
#' @return a S4 'flowSet' object
#'
#' @export 
#' 
#' @importFrom flowCore read.flowSet arcsinhTransform transformList transform
load.flowSet <- function (Results = NULL, fcs.files, exclude.markers, trans, 
    pattern = ".fcs.density.fcs.cluster.fcs") 
{
    message("FCS files loading:")
    if (!is.null(Results)) {
        if (!is.null(Results@fcs.files)) {
            fcs.files <- Results@fcs.files
        }
        else {
            stop("Error in load.flowSet: The 'Results' parameter required a 'Results' object with a not null 'fcs.files' slot")
        }
    }
    flowset <- flowCore::read.flowSet(fcs.files, emptyValue = TRUE)
    samples.names <- gsub(pattern, "", basename(fcs.files))
    flowCore::sampleNames(flowset) <- samples.names
    dictionary <- flowset[[1]]@parameters@data[, c(1, 2)]
    dictionary[, 1] <- make.names(dictionary[, 1])
    dictionary[is.na(dictionary[, 2]), 2] <- dictionary[is.na(dictionary[, 
        2]), 1]
    colnames(flowset) <- rename.markers(colnames(flowset), dictionary = dictionary)
    if (!is.null(Results)) {
        exclude.markers <- setdiff(colnames(flowset), c(Results@marker.names, 
            "cluster"))
    }
    if (!is.null(exclude.markers)) {
        flowset <- exclude.markers(flowset, exclude.markers, 
            colnames.FCS = colnames(flowset))
    }
    if ((is.null(Results) && trans=="arcsinh") || ((!is.null(Results)) && 
        !Results@trans=="arcsinh")) {
        message("\tarchsin transform...")
        transform.arcsinh <- flowCore::arcsinhTransform(a = 0, 
            b = 0.2)
        marker.toTransform <- setdiff(colnames(flowset), "cluster")
        transformations <- flowCore::transformList(marker.toTransform, 
            transform.arcsinh)
        flowset <- flowCore::transform(flowset, transformations)
    }
   if ((is.null(Results) && trans=="logicle") || ((!is.null(Results)) && 
        !Results@trans=="logicle")) {
        message("\tlog transform...")
        transform.logicle <- flowCore::logicleTransform()
        marker.toTransform <- setdiff(colnames(flowset), "cluster")
        transformations <- flowCore::transformList(marker.toTransform, 
            transform.logicle)
        flowset <- flowCore::transform(flowset, transformations)
    }
    return(list(flowset = flowset, dictionary = dictionary))
}
