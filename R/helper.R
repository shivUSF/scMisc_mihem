################################################################################
# average expression
################################################################################

#' @title wrapper function for Seurats AverageExpression
#' @description create average expression matrix with Seurat based on `markers.csv`
#' @param par column name in markers.csv
#' @param object Seurat object
#' @param assay which assay to use
#' @param slot which slot to use
#' @param ortho convert to orthologues? Allowed values: `none`, `mouse2human` or `human2mouse`
#' @return matrix with genes as rows and identity clases as columns
#' @examples \dontrun{szabo_tc_tc_avg <- avgExp("szabo_tc", object = sc_tc_fil, assay = "RNA", slot = "data")}
#' @export 

avgExp <- function(par, object, assay, slot, ortho = "none") {
    if(!file.exists("markers.csv")) {
        stop("Please make sure that markers.csv file exists")
    }
    if(!methods::is(object) == "Seurat") {
        stop("Object must be a Seurat object")
    }
    markers <- readr::read_csv(here::here("markers.csv")) |>
    as.list(markers) |>
    lapply(function(x) x[!is.na(x)])
    genes <- markers[[par]]
    if(is.null(genes)) {
        stop("No genes were found. Make sure that `par` exists in `markers.csv`")
}

    if(!(ortho %in% c("none", "mouse2human", "human2mouse"))) {
        stop("ortho must take values: `none`, `mouse2human` or `human2mouse`")
    }
    if (ortho == "mouse2human") {
        genes <- homologene::mouse2human(genes, db = homologene::homologeneData2)$humanGene
    }
    if (ortho == "human2mouse") {
        genes <- homologene::human2mouse(genes, db = homologene::homologeneData2)$mouseGene
    }
    if (ortho == "none") {
    genes <- genes
    }

    Seurat::AverageExpression(object, features = genes, slot = slot, assay = assay)[[assay]]
}

################################################################################
# scale rows from pheatmap
################################################################################

#' @title scale rows from pheatmap only internal use
#' @description copied from pheatmap for internal use
#' @param x matrix input
#' @return matrix scaled rows

scale_rows <- function(x) {
    m <- apply(x, 1, mean, na.rm = T)
    s <- apply(x, 1, sd, na.rm = T)
    return((x - m) / s)
}


################################################################################
# scale mat from pheatmap
################################################################################

#' @title scale matrix from pheatmap only internal use
#' @description copied from pheatmap for internal use
#' @param mat matrix input
#' @param scale should the values be centered and scaled in row, column direction? Allowed: `row`, `column`, `none`
#' @return matrix scaled 

scale_mat <- function(mat, scale) {
    if(!(scale %in% c("none", "row", "column"))) {
        stop("scale argument shoud take values: 'none', 'row' or 'column'")
    }
    mat <- switch(scale, none = mat, row = scale_rows(mat), column = t(scale_rows(t(mat))))
    return(mat)
}


