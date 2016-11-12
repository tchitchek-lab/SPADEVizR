#' @title Removing of cell clusters from a Results object
#' 
#' @description 
#' This function is used to remove one or more clusters from a Results object.
#' 
#' @param Results a Results object
#' @param clusters a character vector containing the names of the clusters to remove
#'  
#' @return a Results object
#'
#' @export
removeClusters <- function(Results, clusters){
	
	newResults <- Results
	
	if(!is.character(clusters)){
		stop(paste("error in removeCluster: clusters is not character"))
	}
	
	if(any(!clusters %in% newResults@cluster.names)){
		unknown <- clusters[!clusters %in% newResults@cluster.names]
		stop(paste("error in removeCluster: clusters ",unknown," are unknown"))
	}
	
	tokeep <- newResults@cluster.names[!(newResults@cluster.names %in% clusters)]
	
	newResults@cluster.names      <- tokeep
	newResults@cluster.abundances <- newResults@cluster.abundances[tokeep,]
	rownames(newResults@cluster.abundances) <- newResults@cluster.names
	newResults@cluster.phenotypes <- newResults@cluster.phenotypes[newResults@cluster.phenotypes$cluster %in% tokeep,]
	newResults@cluster.number     <- apply(newResults@cluster.abundances,1,sum)
	
	return(newResults)
}


#' @title Merging of cell clusters from a Results object
#' 
#' @description 
#' This function is used to merge one or more clusters in a Results object.
#' 
#' @details
#' The function merges the abundances and the phenotypes of clusters into a new one.
#' Clusters to merge are removed from the Results object.
#' 
#' @param Results a Results object
#' @param clusters a character vector containing the names of the clusters to remove
#' @param name a character specifiyng the name of the new cluster to create
#'  
#' @return a Results object
#'
#' @export
#'
#' @import plyr
mergeClusters <- function(Results, clusters, name){

	newResults <- Results
	
	if(!is.character(clusters)){
		stop(paste("error in mergeClusters: clusters is not character"))
	}
	
	if(!is.character(name)){
		stop(paste("error in mergeClusters: name is not character"))
	}
	
	if(any(!clusters %in% newResults@cluster.names)){
		stop(paste("error in mergeClusters: clusters ",clusters," are unknown"))
	}
	
	tomerge <- newResults@cluster.names[newResults@cluster.names %in% clusters]
	
	newResults@cluster.names      <- c(newResults@cluster.names,name)
	abundances_merged             <- apply(newResults@cluster.abundances[tomerge,],2,sum)
	newResults@cluster.abundances <- rbind(newResults@cluster.abundances,abundances_merged)
	rownames(newResults@cluster.abundances) <- newResults@cluster.names
	phenotypes_merged             <- newResults@cluster.phenotypes[newResults@cluster.phenotypes$cluster %in% tomerge,]
	phenotypes_merged             <- plyr::ddply(phenotypes_merged,"cluster") 
	newResults@cluster.phenotypes <- rbind(newResults@cluster.phenotypes,phenotypes_merged)
	newResults@cluster.number     <- c(newResults@cluster.number,sum(abundances_merged))
	
	newResults <- removeCluster(newResults,clusters)
	
	return(newResults)
}
