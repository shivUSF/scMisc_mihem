################################################################################
# clustifyr UMAP and heatmap
################################################################################

#' @title clustifyr UMAP heatmap wrap
#' @description create and save clustifyr UMAP and heatmap
#' @param my_object query input object
#' @param my_annotations annotations of query input object
#' @param ref_object reference object
#' @param ref_annotations reference annotations
#' @param format format of reference object: possible values are seurat, sc or matrix
#' @param ortho convert to orthologues? Allowed values: `none`, `mouse2human` or `human2mouse` (default: `none`)
#' @param filter should the reference be filtered? Allowed values: TRUE, FALSE (default: FALSE)
#' @param regex if filter: specify filter term
#' @param umap also create UMAP? Allowed values: TRUE, FALSE (default FALSE)
#' @param width width of heatmap plot (default: 10)
#' @param height height of heatmap plot (default: 7)
#' @param cluster_rows cluster rows? (default: true)
#' @param cluster_cols cluster columns? (default: true)
#' @param cellwidth individual cell width (default: 10)
#' @param cellheight individual cell height (default: 10)
#' @param fontsize fontsize (default: 10)
#' @param color color palette
#' @param name fontsize (default: 10)
#' @return save umap/heatmap plot to folder `/results/clustifyr/`
#' @examples
#' \dontrun{clustifyFun(my_object = tcells,
#'            my_annotations = "cluster",
#'            ref_object = projectil_til,
#'            ref_annotations = "functional.cluster",
#'            format = "seurat",
#'            ortho = "none",
#'            height = 25,
#'            width = 10,
#'            umap = TRUE,
#'            color = col_vector,
#'            name = "projectil_til")
#' }
#' @export

clustifyFun <- function(my_object,
                        my_annotations,
                        ref_object,
                        ref_annotations,
                        format,
                        ortho = "none",
                        filter = FALSE,
                        regex = "",
                        umap = FALSE,
                        width = 10,
                        height = 7,
                        cluster_rows = TRUE,
                        cluster_cols = TRUE,
                        cellwidth = 10,
                        cellheight = 10,
                        fontsize = 10,
                        color,
                        name) {
    if(!(ortho %in% c("none", "mouse2human", "human2mouse"))) {
        stop("ortho must take values: `none`, `mouse2human` or `human2mouse`")
        }
    dir.create(file.path("results", "clustifyr"), showWarnings = FALSE)
    if(format == "sce") {
        ref <- clustifyr::object_ref(
            input = ref_object,               # SCE object
            cluster_col = ref_annotations,       # name of column in colData containing cell identities
        )
    } else if (format == "seurat") {
        ref <- clustifyr::seurat_ref(
        seurat_object = ref_object,
        cluster_col = ref_annotations)

    } else if (format == "matrix") {
        ref <- ref_object
    } else {
    stop("please specify input format: sce, seurat or matrix")
    }
    if(filter) {
        ref <- ref[,str_detect(colnames(ref), regex)]
    }
    if(ortho == "mouse2human") {
        lookup <- homologene::mouse2human(rownames(ref), db = homologene::homologeneData2)
        rownames(ref) <- lookup$humanGene[match(rownames(ref), lookup$mouseGene)]
        message("genes converted from mouse to human")
    } else if(ortho == "human2mouse") {
        lookup <- homologene::human2mouse(rownames(ref), db = homologene::homologeneData2)
        rownames(ref) <- lookup$mouseGene[match(rownames(ref), lookup$humanGene)]
        message("genes converted from human to mouse")
    } else if(ortho == "none") {
        message("no genes were converted")
    }
    res_seurat <-clustifyr::clustify(
                                input = my_object,       # a Seurat object
                                ref_mat = ref,    # matrix of RNA-seq expression data for each cell type
                                cluster_col = my_annotations, # name of column in meta.data containing cell clusters
                                compute_method = "spearman",
                                obj_out = TRUE      # output Seurat object with cell type inserted as "type" column
                            )
        Idents(res_seurat) <- res_seurat$type
        if(umap){
            DimPlot(res_seurat, label = TRUE, cols = color)
            ggsave(file.path("results", "clustifyr", glue::glue("clustifyr_{name}_umap.pdf")), width = 10, height = 7)
        }
        res_matrix <-clustifyr::clustify(
                                    input = my_object,       # a Seurat object
                                    ref_mat = ref,    # matrix of RNA-seq expression data for each cell type
                                    cluster_col = my_annotations, # name of column in meta.data containing cell clusters
                                    obj_out = FALSE      # output Seurat object with cell type inserted as "type" column
                                )
        min_break <- round(min(res_matrix), 2)
        max_break <- round(max(res_matrix), 2)
        phmap <- pheatmap::pheatmap(res_matrix,
                                    color = viridis::viridis(100),
                                    scale = "none",
                                    cluster_rows = cluster_rows,
                                    cluster_cols = cluster_cols,
                                    clustering_method = "complete",
                                    border_color = NA,
                                    treeheight_row = 10,
                                    treeheight_col = 10,
                                    cellwidth = cellwidth,
                                    cellheight = cellheight,
                                    fontsize = fontsize,
                                    breaks = seq(min_break, max_break, length = 100),
                                    legend_breaks = round(seq(min_break, max_break, length = 3), 2))
        pdf(file.path("results", "clustifyr", glue::glue("clustifyr_{name}_hmap.pdf")), width = width, height = height)
        print(phmap)
        dev.off()
}
