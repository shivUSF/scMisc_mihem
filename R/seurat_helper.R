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
    markers <- readr::read_csv("markers.csv") |>
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
#' @param min_pct minimum fraction of cells in either two of the populations (default: 0.1)
#' @param logfc_threshold minimum x-fold difference (default: 0.25)
#' @param assay which assay to use in DE testing (e.g. RNA or SCT)
#' @return data frame with significant DE genes arranged by log2FC
#' @examples \dontrun{findMarkersPresto(ident1 = "biopsy", ident2 = "blood", object = sc_tc_fil, assay = "RNA")}
#' @export 

findMarkersPresto <- function(ident1, ident2 = NULL, object, only_pos = FALSE, min_pct = 0.1, logfc_threshold = 0.25, assay = assay) {
    if(!methods::is(object) == "Seurat") {
        stop("Object must be a Seurat object")
    }
    if(is.null(assay)) {
    stop("Please provie assay information")
    }
    result <- SeuratWrappers::RunPresto(object, ident.1 = ident1, ident.2 = ident2, min.pct = min_pct, logfc.threshold = logfc_threshold, only.pos = only_pos, assay = assay) |>
    tibble::rownames_to_column("gene") |>
    dplyr::filter(p_val_adj < 0.05) |>
    dplyr::relocate(gene, avg_log2FC, p_val, p_val_adj) |>
    dplyr::arrange(desc(avg_log2FC))
return(result)
}

################################################################################
# table abundance
################################################################################

#' @title wrapper function to calculate absolute and relative abundance
#' @description wrapper function to calculate absolute and relative abundance
#' @param object Seurat object
#' @param row_var variable in meta data that will represent the rows
#' @param col_var variable in meta data that will represent the columns
#' @examples \dontrun{abundanceTbl(object = aie, row_var = "predicted.id", col_var = "AIE_type")}
#' @export 

abundanceTbl <- function(object, row_var, col_var) {
    if(!methods::is(object) == "Seurat") {
        stop("Object must be a Seurat object")
    }
    if(!dir.exists(file.path("results", "abundance"))) {
        stop("Directory `/results/abundance/` must exist")
    }
    object_parse <- deparse(substitute(object))
    result_abs <- as.data.frame.matrix(table(object@meta.data[[row_var]], object@meta.data[[col_var]])) |>
        rownames_to_column("cell") 

    result_pct <- result_abs |>
        mutate(across(where(is.numeric), function(x) x/sum(x)*100)) |>
        mutate(across(where(is.numeric), function(x) round(x, 2)))

    writexl::write_xlsx(list("absolute" = result_abs, "percentage" = result_pct), file.path("results", "abundance", glue::glue("abundance_tbl_{object_parse}_{col_var}.xlsx")))
}
