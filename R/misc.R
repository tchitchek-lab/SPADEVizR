#' @title Removing of cell clusters from a Results object
#' 
#' @description 
#' This function is used to remove one or more cell clusters from a Results object.
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
		stop(paste("error in removeClusters: clusters is not character"))
	}
	
	if(any(!clusters %in% newResults@cluster.names)){
		unknown <- clusters[!clusters %in% newResults@cluster.names]
		stop(paste("error in removeClusters: clusters ",unknown," are unknown"))
	}
	
	tokeep <- newResults@cluster.names[!(newResults@cluster.names %in% clusters)]
	
	newResults@cluster.names      <- tokeep
	newResults@cluster.abundances <- newResults@cluster.abundances[tokeep,]
	rownames(newResults@cluster.abundances) <- newResults@cluster.names
	newResults@cluster.phenotypes <- newResults@cluster.phenotypes[newResults@cluster.phenotypes$cluster %in% tokeep,]
	newResults@cluster.number     <- apply(newResults@cluster.abundances,1,sum)
	
	invisible(newResults)
}


#' @title Merging of cell clusters from a Results object
#' 
#' @description 
#' This function is used to merge one or more clusters in a Results object.
#' 
#' @details
#' The function merges the abundances and the phenotypes of clusters into a new one.
#' Clusters to merge are not removed from the Results object after the merging.
#' 
#' @param Results a Results object
#' @param clusters a character vector containing the names of the clusters to merge
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
	
	invisible(newResults)
}


#' @title Annotating cell clusters
#' 
#' @description 
#' This function is used to annotate cell clusters based on their marker expression categories.
#' The annotations of cell clusters must be specifyed using a character dataframe indicating for each marker of each population the matching categorial values.
#' 
#' @details
#' The annotation data.frame must have the annotations in rows and the markers in colomuns.
#' Matching values must be specified as a character list of acceptable values. 
#' 
#' @param Results a Results object
#' @param annotations a character dataframe specifying the annotations 
#' @param num a numeric value specifying the number of markers expression categories to be used
#'
#' @return a Results object
#'
#' @export
annotateClusters <- function(Results, annotations, num=5){
	
	newResults <- Results
	
	if(!is.data.frame(annotations)){
		stop(paste("error in annotateClusters: annotations is not a dataframe"))
	}
	
	if(any(!colnames(annotations) %in% newResults@marker.names)){
		unknown <- colnames(annotations)[!colnames(annotations) %in% newResults@marker.names]
		stop(paste("error in annotateClusters: markers in annotations dataframe (",unknown,") are not in markers of the Results object"))
	}
	
	hm <- computePhenoTable(newResults@cluster.phenotypes,newResults@bounds,num=num)
	hm <- reshape2::dcast(hm,marker~cluster)
	rownames(hm) <- hm$marker
	hm$marker    <- NULL
	
	for(cluster in colnames(hm)){
		pheno_hm <- hm[,cluster]
		names(pheno_hm) <- rownames(hm)
		
		for(annotation in rownames(annotations)){
			pheno_annot <- annotations[annotation,]
		
			chk = TRUE
			for(marker in colnames(pheno_annot)[!is.na(pheno_annot)]){
				toeval = paste0(pheno_hm[marker], " %in% ", pheno_annot[marker])
				if(!eval(parse(text=toeval))){
					chk = FALSE
				}
			}
			
			if(chk==TRUE){
				newname <- paste0(cluster,":",annotation) 
				
				if(!is.null(newname)){
					newResults@cluster.names[newResults@cluster.names==cluster] <- newname
					newResults@cluster.phenotypes$cluster[newResults@cluster.phenotypes$cluster==cluster] <- newname
				}
			}
		}
		
	}
	
	invisible(newResults)
}
