################################################################################
# feature plots
################################################################################

#' @title nice Seurat feature plot
#' @description create and save a nice Seurat feature plot in folder `featureplot`
#' @param object Seurat object
#' @param par column name in markers.csv
#' @param width width of output plot (default: 16)
#' @param height height of output plot (default: length of genes divided by two)
#' @return save feature plot to folder `featureplot`
#' @importFrom ggplot2 theme element_blank element_rect ggsave
#' @examples \dontrun{fPlot(sc_merge, par = "main", width = 5, height = 5)}
#' @export

fPlot <- function(object, par, width = 16, height = ceiling(length(genes)/2)) {
    if(!file.exists("markers.csv")) {
        stop("Please make sure that markers.csv file exists")
    }
    if(!methods::is(object) == "Seurat") {
        stop("Object must be a Seurat object")
    }
    dir.create(here::here("featureplot"), showWarnings = FALSE)
    markers <- readr::read_csv(here::here("markers.csv")) |>
    as.list(markers) |>
    lapply(function(x) x[!is.na(x)])
    genes <- markers[[par]]
    if(is.null(genes)) {
        stop("No genes were found. Make sure that `par` exists in markers.csv")
}
    object_parse <- deparse(substitute(object))
    fp <- Seurat::FeaturePlot(object = object, features = unique(genes), cols = viridis::viridis(100), reduction = "umap", pt.size = .1, order = TRUE, coord.fixed = TRUE) +
        theme(axis.text = element_blank(),
                       axis.ticks = element_blank(),
                       panel.border = element_rect(color = "black", size = 1, fill = NA))
    ggsave(here::here("featureplot", glue::glue("fp_{object_parse}_{par}.png")), width = width, height = height, limitsize = FALSE)
}


################################################################################
# dot plots
################################################################################

#' @title nice Seurat dot plot
#' @description create and save a nice Seurat dot plot in folder `featureplot`
#' @param object Seurat object
#' @param par column name in markers.csv
#' @param dot_min minimal dot size
#' @param width width of output plot (default: 10)
#' @param height height of output plot (default: 10)
#' @param ortho boolean; convert mouse to human?
#' @return save feature plot to folder `featureplot`
#' @importFrom ggplot2 ggplot scale_size theme xlab ylab element_text ggsave
#' @examples
#' \dontrun{
#' dotPlot(object = sc_merge, par = "cellmarkers_covid", dot_min = 0.1)
#' }
#' @export

dotPlot <- function(object, par, dot_min, width = 10, height = 10, ortho = FALSE) {
    if(!file.exists("markers.csv")) {
        stop("Please make sure that markers.csv file exists")
    }
    if(!methods::is(object) == "Seurat") {
        stop("Object must be a Seurat object")
    }
    dir.create(here::here("dotplot"), showWarnings = FALSE)
    markers <- readr::read_csv(here::here("markers.csv")) |>
        as.list(markers) |>
        lapply(function(x) x[!is.na(x)])
    genes <- markers[[par]]
    if(is.null(genes)) {
        stop("No genes were found. Make sure that `par` exists in markers.csv")
}
    if (ortho) {
        genes <- homologene::mouse2human(genes, db = homologeneData2)$humanGene
    }
    object_parse <- deparse(substitute(object))
    Seurat::DotPlot(object, features = unique(genes), dot.scale = 10, scale.by = "size", dot.min = dot_min) +
        viridis::scale_color_viridis(option = "viridis") +
        scale_size(range=c(0,10))+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "italic"))+
        xlab(NULL)+
        ylab(NULL)
    ggsave(here::here("dotplot", glue::glue("dp_{object_parse}_{par}.pdf")), width = width, height = height, limitsize = FALSE)
}


################################################################################
# pheatmap
################################################################################

#' @title nice pheatmap
#' @description create and save a nice pheatmap in folder `heatmapplot` using color breakes
#' @param matrix numeric matrix of the values to be plotted
#' @param scale should the values be centered and scaled in row, column direction? Allowed: `row`, `column`, `none` (default: `none`)
#' @param width width of output plot (default: number of columns divided by two)
#' @param height height of output plot (default: number of rows divided by three)
#' @param cellwidth individual cell width (default: 10)
#' @param cellheight individual cell height (default: 10)
#' @param treeheight_row height of a tree for rows (default: 10)
#' @param treeheight_col height of a tree for columns (default: 10)
#' @param fontsize fontsize (default: 10)
#' @param cluster_rows cluster rows? (default: true)
#' @param cluster_cols cluster columns? (default: true)
#' @param annotation_row data frame that contains the annotations. Rows in the data and in the annotation are matched using row names. (default: NA)
#' @return save heatmap to folder `featureplot`
#' @examples \dontrun{pHeatmap(szabo_tc_tc_avg, scale = "row", cluster_cols = FALSE)}
#' @export

pHeatmap <- function(matrix, scale = "none", height = ceiling(nrow(matrix)/3), width = ceiling(ncol(matrix)/2), cellwidth = 10, cellheight = 10, treeheight_row = 10, treeheight_col = 10, fontsize = 10, cluster_rows = TRUE, cluster_cols = TRUE, annotation_row = NA) {
    dir.create(here::here("heatmap"), showWarnings = FALSE)
    matrix_parse <- deparse(substitute(matrix))
    matrix <- matrix[!rowSums(matrix) == 0,] # filter rows with only zeros
    break_max <- round(max(abs(c(max(pheatmap:::scale_mat(matrix, scale = scale)), min(pheatmap:::scale_mat(matrix, scale = scale)))))-0.1,1) #use internal function to get scaled matrix and max value for color legend
    break_min <- -break_max
    phmap <- pheatmap::pheatmap(matrix,
                                color = viridis::viridis(100),
                                scale = scale,
                                cellwidth = cellwidth,
                                cellheight = cellheight,
                                fontsize = fontsize,
                                treeheight_row = treeheight_row,
                                treeheight_col = treeheight_col,
                                cluster_rows = cluster_rows,
                                cluster_cols = cluster_cols,
                                clustering_method = "complete",
                                border_color = NA,
                                legend_breaks = seq(break_min, break_max, length = 3),
                                annotation_row = annotation_row
                                )
    pdf(here::here("heatmap", glue::glue("hm_{matrix_parse}.pdf")), width = width, height = height)
    print(phmap)
    dev.off()
}
