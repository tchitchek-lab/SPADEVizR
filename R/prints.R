#' @title Textual previews for all SPADEVizR objects
#'
#' @description 
#' Prints a preview for a SPADEVizR object.
#'
#' @param x a SPADEVizR object
#' 
#' @return none
#'  
#' @name print
#' @rdname print-methods
NULL

#' @rdname print-methods
#' @export
setMethod("print","Results",
    function(x){
        cat("Object class: Results\n")
        cat("Markers: \n")
        cat(paste0(x@marker.names, collapse = "\n"))
        cat("\n")
        if (is.null(x@assignments)) {
            cat("Samples: ")
            cat("\n")
            cat(paste0(x@sample.names, collapse = "\n"))
        }else{
            cat("Samples and theirs contextual assignments: ")
            cat("\n")
            print(x@assignments)
        }
		cat("\n")
    }
)

#' @rdname print-methods
#' @export
setMethod("print","AC",
    function(x){
        cat("Object class: Abundant Clusters (AC)\n")
        cat("Samples: \n")
        cat(paste0(x@sample.names, collapse = "\n"))
        cat("\n")
        cat(paste0("Used: ", ifelse(x@use.percentages, "relative abundances", "absolute abundances")))
        cat("\n")
        cat(paste0("Number of identified AC: ", nrow(x@results[x@results$significant == TRUE, ])))
        cat("\n")
        cat(paste0("Statistical test used: ", x@method))
        cat("\n")
        cat(paste0("Adjusted method used: ", x@method.adjust))
        cat("\n")
        cat(paste0("Theorical value (mu): ", x@mu))
        cat("\n")
        cat(paste0("P-value threshold: ", x@th.pvalue))
        cat("\n")    
    }
)

#' @rdname print-methods
#' @export
setMethod("print","DAC",
    function(x){
        cat("Object class: Differentially Abundant Clusters (DAC)\n")
        cat("Samples of condition 1: \n")
        cat(paste0(x@sample.cond1, collapse = "\n"))
        cat("\n")
        cat("Samples of condition 2: \n")
        cat(paste0(x@sample.cond2, collapse = "\n"))
        cat("\n")
        cat(paste0("Used: ", ifelse(x@use.percentages, "relative abundances", "absolute abundances")))
        cat("\n")
        cat(paste0("Number of identified DAC: ", nrow(x@results[x@results$significant == TRUE, ])))
        cat("\n")
        cat(paste0("Statistical test used: ", x@method))
        cat("\n")
        cat(paste0("Paired: ", x@method.paired))
        cat("\n")
        cat(paste0("Adjusted method used: ", x@method.adjust))
        cat("\n")
        cat(paste0("P-value threshold: ", x@th.pvalue))
        cat("\n")
        cat(paste0("Fold-change threshold: ", x@th.fc))
        cat("\n")
    }
)

#' @rdname print-methods
#' @export
setMethod("print","CC",
    function(x){
        cat("Object class: Correlated Clusters (CC)\n")
        cat("Samples = variables :\n")
        for (sample in x@sample.names) {
            cat(paste0(sample, " = ", x@variable[sample], "\n"))
        }
        cat("\n")
        cat(paste0("Used: ", ifelse(x@use.percentages, "relative abundances", "absolute abundances")))
        cat("\n")
        cat(paste0("Number of identified CC: ", nrow(x@results[x@results$significant == TRUE, ])))
        cat("\n")
        cat(paste0("Statistical test used: ", x@method))
        cat("\n")
        cat(paste0("Adjusted method used: ", x@method.adjust))
        cat("\n")
        cat(paste0("P-value threshold: "))
        cat("\n")
        cat(paste0(" ", x@th.pvalue, collapse = "\n")) ###
        cat("\n")
        cat(paste0("Correlation threshold: "))
        cat("\n")
        cat(paste0(" ", x@th.correlation, collapse = "\n"))
        cat("\n")
    }
)

#' @rdname print-methods
#' @export
setMethod("print","AP",
    function(x){
        cat("Object class: AP\n")
        cat(paste0("Number of class: "))
        cat(paste0(x@class.number, collapse = "; "))
        cat("\n")
        cat(paste0("Classification method used: "))
        cat(paste0(x@method, collapse = "; "))
        cat("\n")
        cat(paste0("Parameter used"))
        cat(paste0(names(x@method.parameter), " = ", x@method.parameter, collapse = "; "))
        cat("\n")
    }
)


#' @title Textual previews for SPADEVizR objects
#'
#' @description 
#' Show a preview for a SPADEVizR object.
#'
#' @param object a SPADEVizR object
#' 
#' @return none
#'  
#' @name show
#' @rdname show-methods
NULL

#' @rdname show-methods
#' @export
setMethod("show","Results",
    definition = function(object){print(object)}
)

#' @rdname show-methods
#' @export
setMethod("show","AC",
    definition = function(object){print(object)}
)

#' @rdname show-methods
#' @export
setMethod("show","DAC",
    definition = function(object){print(object)}
)

#' @rdname show-methods
#' @export
setMethod("show","CC",
    definition = function(object){print(object)}
)

#' @rdname show-methods
#' @export
setMethod("show","AP",
    definition = function(object){print(object)}
)
