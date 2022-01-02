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
# find markers presto
################################################################################

#' @title wrapper function for presto find markers
#' @description find significant DE genes using presto
#' @param ident1 cell population 1
#' @param ident2 cell population 2, if NULL use all (default: NULL)
#' @param object Seurat object
#' @param only_pos only return positive markers? (default: FALSE)
#' @return data frame with significant DE genes arranged by log2FC
#' @examples \dontrun{findMarkersPresto(ident1 = "biopsy", ident2 = "blood", object = sc_tc_fil)}
#' @export 

findMarkersPresto <- function(ident1, ident2 = NULL, object, only_pos = FALSE) {
    if(!methods::is(object) == "Seurat") {
        stop("Object must be a Seurat object")
    }
    result <- SeuratWrappers::RunPresto(object, ident.1 = ident1, ident.2 = ident2, min.pct = 0.1, logfc.threshold = 0.25, only.pos = only_pos) |>
    rownames_to_column("gene") |>
    filter(p_val_adj < 0.05) |>
    relocate(gene, avg_log2FC, p_val, p_val_adj) |>
    arrange(desc(avg_log2FC))
return(result)
}
