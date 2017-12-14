#' @title Visualization for all SPADEVizR objects
#'
#' @description 
#' Generates a graphical representation for all SPADEVizR objects ('Results', 'AC', 'DAC', 'CC', and 'AP' objects).
#'
#' @details 
#' 'Results' objects are represented using the 'heatmapViewer()' function.
#' 'AC' objects are represented using the 'abundantClustersViewer()' function. 
#' 'DAC' objects are represented using the 'volcanoViewer()' function.
#' 'CC' objects are represented using the 'correlatedClustersViewer()' function.
#' 'AP' objects are represented using the 'circlesPackingViewer()' function.
#'
#' @param x a 'Results', 'AC', 'DAC', 'CC' and 'AP' object
#' @param y unused parameter
#' @param ... supplementary parameters transmitted to the 'heatmapViewer()', 'abundantClustersViewer()', 'volcanoViewer()', 'correlatedClustersViewer()' or 'circlesPackingViewer()' functions
#' 
#' @return a 'ggplot' object
#'  
#' @name plot
#' @rdname plot-methods
setGeneric("plot", function(x, y = NULL, ...){ standardGeneric("plot") })

#' @rdname plot-methods
#' @export
setMethod("plot", c("Results", "missing"),
    function(x, ...) {
        return(heatmapViewer(x, ...))
    }
)

#' @rdname plot-methods
#' @export
setMethod("plot", c("AC", "missing"),
    function(x, y, ...){
        return(abundantClustersViewer(x, ...))
    }
)

#' @rdname plot-methods
#' @export
setMethod("plot", c("DAC", "missing"),
    function(x, ...){
        return(volcanoViewer(x, ...))
    }
)

#' @rdname plot-methods
#' @export
setMethod("plot", c("CC", "missing"),
    function(x, y, ...){
        return(correlatedClustersViewer(x, ...))
    }
)

#' @rdname plot-methods
#' @export
setMethod("plot", c("AP", "missing"),
    function(x, ...){
        return(circlesPackingViewer(x, ...))
    }
)

