#' @title Generating SPADEVizR reports
#'
#' @description 
#' Generates a customizable PDF report based on SPADEVizR visualization features.
#' 
#' @details 
#' Available plots are :
#' \itemize{
#' \item {"count" (included by default):}{Displays a Count Viewer representation;}
#' \item {"tree":}{Displays a SPADE Tree Viewer representation;}
#' \item {"heatmap" (included by default):}{Displays an Heatmap Viewer representation;}
#' \item {"boxplot":}{Displays a Boxplot Viewer representation;}
#' \item {"kinetics":}{Displays a Kinetic Viewer representation;}
#' \item {"streamgraph":}{Display a Streamgraph Viewer representation;}
#' \item {"pheno" (included by default):}{Displays a Pheno Viewer representation;}
#' \item {"MDSclusters" (included by default):}{Displays a MDS Viewer representation at the cluster level}
#' \item {"MDSsamples":}{Displays a MDS Viewer representation at the sample level}
#' \item {"distogram":}{Displays a Distogram Viewer representation;}
#' \item {"kinetics_pheno":}{Displays two "kinetics" and "cluster" representations juxtaposed (one on the side of the other) for each cluster;}
#' \item {"boxplot_pheno":}{Displays two "boxplot" and "cluster" representations juxtaposed (one on the side of the other) for each cluster;}
#' \item {AC, DAC, CC and AP objects}
#' }
#' 
#' @param Results a 'Results' object
#' @param PDFfile a character specifying the output path
#' @param select.plots a vector combining character and stat objects ('AC', 'DAC', 'CC' and 'AP') specifying the order of the desired plots (see details) 
#' @param clusters a character vector of clusters to include in the report (all will be included by default)
#' @param markers a character vector of markers to include in the report (all will be included by default)
#' @param samples a character vector providing the sample names to used (all samples by default)
#' @param stat.objects a list containing one or several AC, DEC, CC or AP objects to plot in the report
#' @param width a numeric specifying the width of the PDF file
#' @param height a numeric specifying the height of the PDF file
#' @param verbose a boolean specifying if some messages must be displayed during the generation of the report
#'
#' @return none
#'
#' @export
#' 
#' @import gridExtra
createReport <- function(Results,
                           PDFfile      = "report.pdf",
                           select.plots = c("count", "heatmap", "MDSclusters", "pheno"),
                           clusters     = NULL,
                           markers      = NULL,
                           samples      = NULL,
                           stat.objects = list(),
                           width        = 30,
                           height       = 10,
                           verbose      = TRUE) {

    message("[BEGIN] - generating report")
    on.exit(dev.off())

    if (is.null(clusters)) {
        clusters <- Results@cluster.names
    }

    nb.plot <- length(select.plots)
    i <- 1
    pdf(PDFfile, width = width, height = height)
    for (current.plot in select.plots) {

        current.name <- ifelse(typeof(current.plot) == "character", current.plot, "object")

        if (verbose) {
            message(paste0("\tGenerate: ",
                               ifelse(current.name == "object", class(current.plot)[1], current.name),
                               ", ", i, " on ", nb.plot))
        }
        i <- i + 1

        switch(current.name,
               "count" = {
            countViewer(Results, samples = samples, clusters = clusters, show.on_device = TRUE)
        },
               "MDSclusters" = {
            MDSViewer(Results, space = "clusters", samples = samples, clusters = clusters, show.on_device = TRUE)
        },
               "MDSsamples" = {
            MDSViewer(Results, space = "samples", samples = samples, clusters = clusters, show.on_device = TRUE)
        },
               "heatmap" = {
            heatmapViewer(Results, clusters = clusters, markers = markers, show.on_device = TRUE)
        },
               "tree" = {
            treeViewer(Results, samples = samples, show.on_device = TRUE)
        },
               "kinetics_pheno" = {

            for (j in clusters) {
                if (verbose) {
                    message(paste0("\t\tCluster ", j, " on ", length(clusters)))
                }
                kinetics.plot <- kineticsViewer(Results, samples = samples, clusters = j, show.on_device = FALSE)
                cluster.plot  <- phenoViewer(Results, samples = samples, clusters = j, markers = markers, show.on_device = FALSE)
                gridExtra::grid.arrange(kinetics.plot, cluster.plot, ncol = 2)
            }
        },
               "boxplot_pheno" = {

            for (j in clusters) {
                if (verbose) {
                    message(paste0("\t\tCluster ", j, " on ", length(clusters)))
                }
                boxplot.plot <- boxplotViewer(Results, samples = samples, clusters = j, show.on_device = FALSE)
                cluster.plot <- phenoViewer(Results, samples = samples, clusters = j, markers = markers, show.on_device = FALSE)
                gridExtra::grid.arrange(boxplot.plot, cluster.plot, ncol = 2)
            }
        },
               "boxplot" = {
            for (j in clusters) {
                if (verbose) {
                    message(paste0("\t\tCluster ", j, " on ", length(clusters)))
                }
                boxplotViewer(Results, samples = samples, clusters = j, show.on_device = TRUE)
            }
        },
               "kinetics" = {
            for (j in clusters) {
                if (verbose) {
                    message(paste0("\t\tCluster ", j, " on ", length(clusters)))
                }
                kineticsViewer(Results, samples = samples, clusters = j, show.on_device = TRUE)
            }
        },
                "pheno" = {
            for (j in clusters) {
                if (verbose) {
                    message(paste0("\t\tCluster ", j, " on ", length(clusters)))
                }
                phenoViewer(Results, samples = samples, clusters = j, markers = markers, show.on_device = TRUE)
            }

        },
              "distogram" = {
            distogramViewer(Results, samples = samples, show.on_device = TRUE)
        },
              "streamgraph" = {
            streamgraphViewer(Results, samples = samples, clusters = clusters, show.on_device = TRUE)
        },
              "object" = {
            plot(current.plot, show.on_device = TRUE)
        }, warning(paste0("Unknown plots names ignored:", current.name)))

    }

    message("[END] - generating report")

}

