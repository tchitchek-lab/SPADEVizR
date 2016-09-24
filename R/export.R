#' @title Exportation of SPADEVizR objects
#'
#' @description Exports a SPADEVizR object into a tab separated file.
#'
#' @param object a SPADEVizR object
#' @param filename a character indicating the location of the output file
#'
#' @return none
#'
#' @name export
#' @rdname export-methods
setGeneric("export", function(object, filename = "export.txt") { standardGeneric("export") })

#' @rdname export-methods
#' @export
setMethod("export",c("Results"),
    function(object,filename){
        cat(file = filename, paste0("#Results object: ", names(object),"\n"))
        cat(file = filename, "#sample.names: ", paste0("\"", object@sample.names, "\"", collapse = "\t"), sep = "", "\n", append = TRUE)
        cat(file = filename, "#cluster.names: ", paste0(object@cluster.names, collapse = "\t"), sep = "", "\n", append = TRUE)
        cat(file = filename, "#marker.names: ", paste0("\"", object@marker.names, "\"", collapse = "\t"), sep = "", "\n", append = TRUE)
        cat(file = filename, "#fcs.files:\n", paste0("\"", object@fcs.files, "\"", collapse = "\t"), sep = "", "\n", append = TRUE)
        cat(file = filename, "#cluster.abundances:\n", append = TRUE)
        cat(file = filename, colnames(object@cluster.abundances), append = TRUE)
        cat("\n")
        utils::write.table(object@cluster.abundances, file = filename, append = TRUE, sep = "\t", col.names = NA)
        cat(file = filename, "#cluster.phenotypes:\n", append = TRUE)
        cat(file = filename, colnames(object@cluster.phenotypes), append = TRUE)
        cat("\n")
        utils::write.table(object@cluster.phenotypes, file = filename, append = TRUE, sep = "\t", col.names = NA)
        cat(file = filename, "#quantiles:\n", append = TRUE)
        cat(file = filename, colnames(object@bounds), append = TRUE)
        cat("\n")
        utils::write.table(object@bounds, file = filename, append = TRUE, sep = "\td", col.names = NA)
    }
)

#' @rdname export-methods
#' @export
setMethod("export",c("AC"),
    function(object,filename){
        cat(file = filename, "#AC object (Abundant Clusters)\n")
        cat(file = filename, "#sample names: ", paste0("\"", object@sample.names, "\"", collapse = "\t"), sep = "", "\n", append = TRUE)
        cat(file = filename, "#cluster sizes: ", paste0(object@cluster.size, collapse = "\t"), "\n", sep = "", append = TRUE)
        cat(file = filename, "#use.percentages: ", paste0(object@use.percentages, collapse = "\t"), "\n", sep = "", append = TRUE)
        cat(file = filename, "#mu: ", paste0(object@mu, collapse = "\t"), "\n", sep = "", append = TRUE) 
        cat(file = filename, "#method: ", paste0(object@method, collapse = "\t"), "\n", sep = "", append = TRUE)
        cat(file = filename, "#method.adjust: ", paste0(object@method.adjust, collapse = "\t"), "\n", sep = "", append = TRUE)
        cat(file = filename, "#p-value threshold: ", paste0(object@th.pvalue, collapse = "\t"), "\n", sep = "", append = TRUE)
        cat(file = filename, "#results :\n", append = TRUE)
        cat(file = filename, colnames(object@results), append = TRUE)
        cat(file = filename, "\n", append = TRUE)
        utils::write.table(object@results, file = filename, append = TRUE, sep = "\t", row.names = FALSE, col.names = FALSE)
    }
)

#' @rdname export-methods
#' @export
setMethod("export",c("DAC"),
    function(object,filename){
        cat(file = filename, "#DAC object (Differentially Abundant Clusters)\n")
        cat(file = filename, "#Samples of condition 1: ", paste0("\"", object@sample.cond1, "\"", collapse = "\t"), "\n", sep = "", append = TRUE)
        cat(file = filename, "#Samples of condition 2: ", paste0("\"", object@sample.cond2, "\"", collapse = "\t"), "\n", sep = "", append = TRUE)
        cat(file = filename, "#cluster.size: ", paste0(object@cluster.size, collapse = "\t"), "\n", sep = "", append = TRUE)
        cat(file = filename, "#use.percentages: ", paste0(object@use.percentages, collapse = "\t"), "\n", sep = "", append = TRUE)
        cat(file = filename, "#method: ", paste0(object@method, collapse = "\t"), "\n", sep = "", append = TRUE)
        cat(file = filename, "#method.adjust: ", paste0(object@method.adjust, collapse = "\t"), "\n", sep = "", append = TRUE)
        cat(file = filename, "#method.paired: ", paste0(object@method.paired, collapse = "\t"), "\n", sep = "", append = TRUE)
        cat(file = filename, "#th.pvalue: ", paste0(object@th.pvalue, collapse = "\t"), "\n", sep = "", append = TRUE)
        cat(file = filename, "#th.fc: ", paste0(object@th.fc, collapse = "\t"), "\n", sep = "", append = TRUE)
        cat(file = filename, "#results :\n", append = TRUE)
        cat(file = filename, colnames(object@results), append = TRUE)
        cat(file = filename, "\n", append = TRUE)
        utils::write.table(object@results, file = filename, append = TRUE, sep = "\t", row.names = FALSE, col.names = FALSE)
    }
)

#' @rdname export-methods
#' @export
setMethod("export",c("CC"),
    function(object,filename){
        cat(file = filename, "#CC object (Correlated Clusters)\n")
        cat(file = filename, "#Sample names: ", paste0("\"", object@sample.names,"\"", collapse = "    "), "\n", sep = "", append = TRUE)
        cat(file = filename, "#variable: ", paste0(object@variable, collapse = "\t"), "\n", sep = "", append = TRUE)
        cat(file = filename, "#cluster.size: ", paste0(object@cluster.size, collapse = "\t"), "\n", sep = "", append = TRUE)
        cat(file = filename, "#use.percentages: ", paste0(object@use.percentages, collapse = "\t"), "\n", sep = "", append = TRUE)
        cat(file = filename, "#method: ", paste0(object@method, collapse = "\t"), "\n", sep = "", append = TRUE)
        cat(file = filename, "#method.adjust: ", paste0(object@method.adjust, collapse = "\t"), "\n", sep = "", append = TRUE)
        cat(file = filename, "#th.pvalue: ", paste0(object@th.pvalue, collapse = "\t"), "\n", sep = "", append = TRUE)
        cat(file = filename, "#th.correlation: ", paste0(object@th.correlation, collapse = "\t"), "\n", sep = "", append = TRUE)
        cat(file = filename, "#results :\n", append = TRUE)
        cat(file = filename, colnames(object@results), append = TRUE)
        cat(file = filename, "\n", append = TRUE)
        utils::write.table(object@results, file = filename, append = TRUE, sep = "\t", row.names = FALSE, col.names = FALSE)
    }
)

#' @rdname export-methods
#' @export
setMethod("export",c("AP"),
    function(object,filename){
        cat(file = filename, "#AP object (Abundance Profiles)\n")
        cat(file = filename, "#Number of identified classes: ", paste0(object@class.number, collapse = "\t"), "\n", sep = "", append = TRUE)
        cat(file = filename, "#Method used: ", paste0("\"", object@method, "\"", collapse = "\t"), "\n", sep = "", append = TRUE)
        cat(file = filename, "#Method parameters: ", paste0(object@method.parameter, collapse = "\t"), "\n", sep = "", append = TRUE)
        cat(file = filename, "#Classes:\n", append = TRUE)
        cat(file = filename, colnames(object@classes), append = TRUE)
        cat(file = filename, "\n", append = TRUE)
        utils::write.table(object@classes, file = filename, append = TRUE, sep = "\t", row.names = FALSE, col.names = FALSE)
    }
)
