theme_rect <-function() {
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.border = element_rect(color = "black", size = 1, fill = NA))
}


##############################################################################
##############################################################################
filterSeurat <- function(x, rna, mt) {
subset(x, subset = nFeature_RNA > 200 & nFeature_RNA < rna & percent_mt < mt & scDblFinder.class == "singlet")
}

##############################################################################
##############################################################################

FPlot <- function(object, par, width = 16, height = ceiling(length(genes)/2)) {
markers <- readxl::read_excel("./markers.xlsx") %>%
    as.list(markers) %>%
    lapply(function(x) x[!is.na(x)])
object_parse <- deparse(substitute(object))
genes <- markers[[par]]
fp <- Seurat::FeaturePlot(object = object, features = unique(genes), cols = viridis::viridis(100), reduction = "umap", pt.size = .1, order = TRUE, coord.fixed = TRUE) &
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.border = element_rect(color = "black", size = 1, fill = NA))
ggsave(file.path(project_path, "featureplot", glue::glue("{object_parse}_{par}.png")), width = width, height = height, limitsize = FALSE)
}


##############################################################################
##############################################################################


lss <- function () {
data.frame('object' = ls(".GlobalEnv")) %>% 
      dplyr::mutate(size_unit = object %>%sapply(. %>% get() %>% object.size %>% format(., unit = 'auto')),
                    size = sapply(strsplit(size_unit, split = ' '), FUN = function(x) x[1]),
                    unit = factor(sapply(strsplit(size_unit, split = ' '), FUN = function(x) x[2]), levels = c('Gb', 'Mb', 'Kb', 'bytes'))) %>% 
      dplyr::arrange(unit, dplyr::desc(size)) %>%
      dplyr::select(-size_unit) %>%
    as_tibble()

}

##############################################################################
##############################################################################


dirCreate <- function(folder_name) {
    fs::dir_create(project_path, folder_name)
}

##############################################################################
##############################################################################


dotPlot <- function(object, par, dot.min, width = 10, height = 10,ortho = FALSE) {
markers <- readxl::read_excel("./markers.xlsx") %>%
    as.list(markers) %>%
    lapply(function(x) x[!is.na(x)])
object_parse <- deparse(substitute(object))
genes <- markers[[par]]
if (ortho) {
genes <- homologene::mouse2human(genes, db = homologeneData2)$humanGene
}
Seurat::DotPlot(object, features = unique(genes), dot.scale = 10, scale.by = "size", dot.min = dot.min)+ viridis::scale_color_viridis(option = "viridis")+scale_size(range=c(0,10))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "italic"))+
    xlab(NULL)+
    ylab(NULL)
ggsave(file.path(project_path, "dotplot", glue::glue("dp_{object_parse}_{par}.pdf")), width = width, height = height, limitsize = FALSE)
}

##############################################################################
##############################################################################

avgExp <- function(par, object, assay, slot, ortho = FALSE) {
markers <- readxl::read_excel("./markers.xlsx") %>%
    as.list(markers) %>%
    lapply(function(x) x[!is.na(x)])
genes <- markers[[par]]
if (ortho) {
genes <- homolegene::mouse2human(genes, db = homologeneData2)$humanGene
}
AverageExpression(object, features = genes, slot = slot, assay = assay)[[assay]]
}

##############################################################################
##############################################################################


avgExpScore <- function(par, object, assay, slot) {
markers <- readxl::read_excel("./markers.xlsx") %>%
    as.list(markers) %>%
    lapply(function(x) x[!is.na(x)])
genes <- markers[[par]]
avg_exp <- AverageExpression(object, features = genes, slot = slot, assay = assay)[[assay]]
result <- data.frame(colMeans(avg_exp))
colnames(result) <- par
return(result)
 }

##############################################################################
##############################################################################

phmap <- function(matrix, scale = "none", height = ceiling(nrow(matrix)/3), width = ceiling(ncol(matrix)/2), cellwidth = 10, cellheight = 10, treeheight_row = 10, treeheight_col = 10, fontsize = 10, cluster_rows = TRUE, cluster_cols = TRUE, annotation_row = NA) {
matrix_parse <- deparse(substitute(matrix))
matrix <- matrix[!rowSums(matrix) == 0,] # filter rows with only zeros
break_max <- round(max(abs(c(max(pheatmap:::scale_mat(matrix, scale = scale)), min(pheatmap:::scale_mat(matrix, scale = scale)))))-0.1,1) #use internal function to get scaled matrix and max value for color legend
break_min <- -break_max
phmap <- pheatmap::pheatmap(matrix,
        color = viridis::viridis(100),
#        breaks = c(0.5,0.6),
         scale = scale,
         cellwidth = cellwidth,
         cellheight = cellheight,
        fontsize = fontsize,
        treeheight_row = treeheight_row,
        treeheight_col = treeheight_col,
        cluster_rows = cluster_rows,
        cluster_cols = cluster_cols,
        # clustering_distance_rows = "euclidean",
         #clustering_distance_cols = "euclidean",
         #clustering_distance_rows = dist((1-cor(facs.mean.matrix, method = "spearman"))),
         clustering_method = "complete",
        border_color = NA,
       legend_breaks = seq(break_min, break_max, length = 3),
        annotation_row = annotation_row
        #legend_breaks = pheatmap:::generate_breaks(matrix, 2, center = TRUE)
        # filename = "pheatmap_mean_euclidean.png"
         )
pdf(file.path(project_path, "heatmap", glue::glue("hm_{matrix_parse}.pdf")), width = width, height = height)
print(phmap)
dev.off()
}

##############################################################################
##############################################################################
findMarkersPresto <- function(ident1, ident2 = NULL, object, only_pos = FALSE) {
object_parse <- deparse(substitute(object))
result <-
#    FindMarkers(object, ident.1 = ident1, ident.2 = ident2, min.pct = 0.1, logfc.threshold = 0.25) %>%
SeuratWrappers::RunPresto(object, ident.1 = ident1, ident.2 = ident2, min.pct = 0.1, logfc.threshold = 0.25, only.pos = only_pos) %>%
    rownames_to_column("gene") %>%
    filter(p_val_adj < 0.05) %>%
    relocate(gene, avg_log2FC, p_val, p_val_adj) %>%
    arrange(desc(avg_log2FC))
return(result)
}

##############################################################################
##############################################################################
sbpFun <- function(object, color_by, x_axis, color_order, x_order, width, height, color = color) {
object_parse <- deparse(substitute(object))
result_wide <- as.data.frame.matrix(table(object@meta.data[[color_by]], object@meta.data[[x_axis]])) %>%
    rownames_to_column("cell") %>%
    mutate(across(where(is.numeric), function(x) x/sum(x)*100))
result_long <- result_wide%>%
    pivot_longer(!cell, names_to = "type", values_to = "count") %>%
    mutate(cell = factor(cell, levels = color_order)) %>%
    mutate(type = factor(type, levels = x_order)) %>%
    filter(count != 0)
sbp <- ggplot(data = result_long)+
    geom_col(aes(x = type, y = count, fill = cell), color = "black", size = 0.1, position = "fill")+
    scale_fill_manual(values = color)+
    guides(fill = guide_legend(title = NULL))+ #remove guide label
    theme_classic()+ #remove background
    ylab("Proportion of cells")+
    xlab("")+
    theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.3))
ggsave(file.path(project_path, "abundance", glue::glue("stacked_barplot_{object_parse}_{x_axis}.pdf")), sbp, width = width, height = height)
return(sbp)
}

##############################################################################
##############################################################################

abundanceTbl <- function(object, row_var, col_var) {
object_parse <- deparse(substitute(object))
result_abs <- as.data.frame.matrix(table(object@meta.data[[row_var]], object@meta.data[[col_var]])) |>
    rownames_to_column("cell") 

result_pct <- result_abs |>
    mutate(across(where(is.numeric), function(x) x/sum(x)*100)) |>
    mutate(across(where(is.numeric), function(x) round(x, 2)))

writexl::write_xlsx(list("absolute" = result_abs, "percentage" = result_pct), file.path(project_path, "abundance", glue::glue("abundance_tbl_{object_parse}_{col_var}.xlsx")))
}
##############################################################################
##############################################################################
VolPlotFun <- function(object, cluster_idents, group_by, group1, group2, cluster_order, lookup, my_color) {
object_parse <- deparse(substitute(object))

cl_size_ind <- as.data.frame.matrix(table(object@meta.data[[cluster_idents]], object@meta.data$sample)) |>
    rownames_to_column("cluster") |>
    mutate(across(where(is.numeric), function(x) x/sum(x)*100)) |>
    pivot_longer(!cluster, names_to = "sample", values_to = "count") |>
#    mutate(cell = factor(cell, levels = cluster_idents_order)) |>
 #   mutate(type = factor(type, levels = group_by_order)) |>
#    dplyr::filter(count != 0) |>
    left_join(lookup, by = "sample") |>
#    dplyr::filter(pool_names == group1 | pool_names == group2)
    dplyr::filter(!!sym(group_by) == group1 | !!sym(group_by) == group2) |>
    dplyr::rename(col_interest = !!sym(group_by))

pvalue_res <- vector("double") # define output

for (i in cluster_order) {
    out1 <- cl_size_ind[cl_size_ind$cluster == i,]
    pvalue_res[i] <- wilcox.test(count ~ col_interest, data = out1)$p.value
}

#wilcox_res <- p.adjust(wilcox_res, "BH")
pvalue_cl <- data.frame(cluster = cluster_order, pvalue = pvalue_res)

cl_size <- as.data.frame.matrix(table(object@meta.data[[cluster_idents]], object@meta.data[[group_by]])) |>
    rownames_to_column("cluster") |>
    mutate(across(where(is.numeric), function(x) x/sum(x)*100)) |>
    mutate(logratio = log2(!!sym(group1)/!!sym(group2))) |>
    left_join(pvalue_cl, by = "cluster") |>
    mutate(log_pvalue = -log10(pvalue))|>
    mutate(cluster = factor(cluster, levels = cluster_order))

ggplot(cl_size, aes(x = logratio, y = log_pvalue, color = cluster, size = 3, label = cluster))+
    geom_point()+
    scale_color_manual(values = my_color)+
    theme_classic()+
#    geom_text(vjust = 0, nudge_y = 0.05)+
    ggrepel::geom_text_repel(nudge_y = 0.07)+
    geom_hline(yintercept = -log10(0.05), color = "blue", linetype = "dashed")+ #horizontal line p unadjusted
    #geom_hline(yintercept = -log10(0.05/nrow(out)), color = "blue")+
    geom_vline(xintercept = 0, color = "red", linetype = "dashed")+ #vertical line
    xlab(bquote(~Log[2]~ 'fold change'))+
    ylab(bquote(-Log[10]~ "p value"))+
    theme(legend.position = "none") #remove guide
ggsave(file.path(project_path, "abundance", glue::glue("volcano_{object_parse}_{group1}_{group2}.pdf")), width = 5, height = 5)
}


##############################################################################
##############################################################################
DotClPlot <- function(object, color_by, x_axis, cluster1, cluster2, color_order, color, width = 4, height = 4) {
object_parse <- deparse(substitute(object))
out <- as.data.frame.matrix(table(object@meta.data[[color_by]], object@meta.data[[x_axis]])) %>%
    rownames_to_column("cell") %>%
    mutate(across(where(is.numeric), function(x) x/sum(x)*100)) %>%
    mutate(logratio = log2(!!sym(cluster1)/!!sym(cluster2))) %>%
    mutate(cell = factor(cell, levels = color_order))
ggplot(out, aes(x = logratio, y = fct_reorder(cell, logratio), color = cell))+
    geom_point(size = 5)+
    theme_classic()+
    geom_vline(xintercept = 0, color = "red", linetype = "dashed")+ #vertical line
    scale_color_manual(values = color)+
    xlab("Log2 fold change")+
    ylab(NULL)+
    theme(legend.position="none") # remove legend
ggsave(file.path(project_path, "abundance", glue::glue("dotplot_{object_parse}_{cluster1}_{cluster2}.pdf")), width = width, height = height)
}

##############################################################################
##############################################################################
barAbundance <- function(object, color_by, x_axis, cluster1, cluster2, color_order, color, width = 4, height = 4) {
object_parse <- deparse(substitute(object))
out <- as.data.frame.matrix(table(object@meta.data[[color_by]], object@meta.data[[x_axis]])) %>%
    rownames_to_column("cell") %>%
    mutate(across(where(is.numeric), function(x) x/sum(x)*100)) %>%
    #mutate(logratio = log2(!!sym(cluster1)/!!sym(cluster2))) %>%
    mutate(difference = !!sym(cluster1)-!!sym(cluster2)) %>%
    mutate(cell = factor(cell, levels = color_order))
ggplot(out, aes(x = logratio, y = fct_reorder(cell, logratio), fill = cell))+
#ggplot(out, aes(x = difference, y = fct_reorder(cell, difference), fill = cell))+
    geom_col()+
    theme(panel.border = element_rect(color = "black", size = 1, fill = NA),
          panel.background = element_blank(),
          legend.position = "none")+
    scale_fill_manual(values = color)+
   xlab("Log2 fold change")+
    #xlab("Difference (%g)")+
    ylab(NULL)
ggsave(file.path(project_path, "abundance", glue::glue("barplot_{object_parse}_{cluster1}_{cluster2}_diff.pdf")), width = width, height = height)
}

##############################################################################
##############################################################################
moduleScore <- function(module, object, assay) {
object_parse <- deparse(substitute(object))
markers <- readxl::read_excel("./markers.xlsx") %>%
    as.list(markers) %>%
    lapply(function(x) x[!is.na(x)])
object_new <- AddModuleScore(object, features = markers[module], assay = assay, name = module)
colnames(object_new@meta.data) <- stringr::str_replace(string = colnames(object_new@meta.data), pattern = glue::glue("{module}1"), replacement = module)
FeaturePlot(object = object_new, features = module, cols = viridis::viridis(100), reduction = "umap", pt.size = .1, order = TRUE, coord.fixed = TRUE, min.cutoff = 0)+ ggtitle(module) + 
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.border = element_rect(color = "black", size = 1, fill = NA))
ggsave(file.path(project_path, "featureplot", glue::glue("module_score_{module}.pdf")), width = 6, height = 6)
VlnPlot(object = object_new, features = module)+ ggtitle(module)
ggsave(file.path(project_path, "violin", glue::glue("module_score_{object_parse}_{module}.pdf")), width = 10, height = 10)
}


##############################################################################
##############################################################################
clustifyList <- function(object, marker, seurat_col = "seurat_clusters", name, cluster_rows = TRUE, cluster_cols = TRUE, width = 10, height = 7, filter = FALSE, regex){
if(filter) {
        marker <- marker[,str_detect(colnames(marker), regex)]
}
res_seurat <- clustifyr::clustify_lists(
  input = object,             # matrix of normalized single-cell RNA-seq counts
  marker = marker,                 # list of known marker genes
  cluster_col = seurat_col,
 metric = "jaccard",                   # test to use for assigning cell types
  seurat_out = TRUE)
Idents(res_seurat) <- res_seurat$type

DimPlot(res_seurat, label = TRUE)
ggsave(fs::path(project_path, "clustifyr", glue::glue("clustifyr_{name}_umap"), ext = "pdf"), width = 15, height = 10)

res_list <- clustifyr::clustify_lists(
  input = object,             # matrix of normalized single-cell RNA-seq counts
  marker = marker,                 # list of known marker genes
  cluster_col = seurat_col,
  seurat_out = FALSE,
  metric = "jaccard")
min_break <- round(min(res_list), 2)
max_break <- round(max(res_list), 2)
    pdf(file.path(project_path, "clustifyr", glue::glue("clustifyr_{name}_hmap.pdf")), width = width, height = height)
print(pheatmap::pheatmap(res_list,color = viridis::viridis(100), scale = "none", cluster_rows = cluster_rows, cluster_cols = cluster_cols, clustering_method = "complete", border_color = NA, treeheight_row = 10, treeheight_col = 10, cellwidth = 50, cellheight = 50, fontsize = 20, breaks = seq(min_break, max_break, length = 100), legend_breaks = round(seq(min_break, max_break, length = 3), 2)))
    dev.off()

}

##############################################################################
##############################################################################
clustifyFun <- function(ref_object, ref_annotations, my_object, my_annotations, name, width = 10, height = 7, format, ortho = FALSE,  filter = FALSE, regex, umap = FALSE, cluster_rows = TRUE, cluster_cols = TRUE, cellwidth = 50, fontsize = 20, color = col_vector){
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
    if(ortho) {
        lookup <- homologene::mouse2human(rownames(ref), db = homologene::homologeneData2)
    rownames(ref) <- lookup$humanGene[match(rownames(ref), lookup$mouseGene)]
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
    ggsave(fs::path(project_path, "clustifyr", glue::glue("clustifyr_{name}_umap"), ext = "pdf"), width = 15, height = 10)
    }
    res_matrix <-clustifyr::clustify(
                                input = my_object,       # a Seurat object
                                ref_mat = ref,    # matrix of RNA-seq expression data for each cell type
                                cluster_col = my_annotations, # name of column in meta.data containing cell clusters
                                obj_out = FALSE      # output Seurat object with cell type inserted as "type" column
                            )
min_break <- round(min(res_matrix), 2)
max_break <- round(max(res_matrix), 2)
phmap <- pheatmap::pheatmap(res_matrix,color = viridis::viridis(100), scale = "none", cluster_rows = cluster_rows, cluster_cols = cluster_cols, clustering_method = "complete", border_color = NA, treeheight_row = 10, treeheight_col = 10, cellwidth = cellwidth, fontsize = fontsize, breaks = seq(min_break, max_break, length = 100), legend_breaks = round(seq(min_break, max_break, length = 3), 2))
    pdf(file.path(project_path, "clustifyr", glue::glue("clustifyr_{name}_hmap.pdf")), width = width, height = height)
    #print(clustifyr::plot_cor_heatmap(cor_mat = res_matrix, col = viridis::viridis(100), column_names_gp = grid::gpar(fontsize = col_fontsize), row_dend_width = unit(3, "mm"), column_dend_height = unit(3, "mm")))
print(phmap)
    dev.off()
#return(ggplotify::as.ggplot(phmap))
}


##############################################################################
##############################################################################

volcanoPlot <- function(de_list, fc_cutoff = 0.5, p_cutoff = 1e-30, selectLab = NULL, drawConnectors = TRUE, width, height) {
stopifnot(is.character(de_list))
sheets <- readxl::excel_sheets(file.path(project_path, "de", glue::glue("de_{de_list}.xlsx")))
input_list <- list()
for (x in sheets) {
    input_list[[x]] <- readxl::read_excel(file.path(project_path, "de", glue::glue("de_{de_list}.xlsx")), sheet = x)
    if(nrow(input_list[[x]]) != 0){ #only if tibble is not empty
        volcano <- EnhancedVolcano::EnhancedVolcano(data.frame(input_list[[x]]),
                                   lab = input_list[[x]][[1]],
                                   x = 'avg_log2FC',
                                   y = 'p_val_adj',
                                   pCutoff = p_cutoff,
                                   FCcutoff = fc_cutoff,
                                   axisLabSize = 20,
                                   pointSize = 5,
                                   labSize = 7,
                                   legendLabSize = 20,
                                   legendIconSize = 10,
                                   subtitle = NULL,
                                   caption = NULL,
                                   drawConnectors = TRUE,
                                   title = NULL,
                                        #    title = paste(x, input),
                                   boxedLabels = TRUE,
                                   selectLab = selectLab[[x]],
                                        #xlim=c(0, 2),
                                        #  ylim =c(0,50),
                                   xlab = bquote(~Log[2]~ 'fold change'),
                                   ylab = bquote(-Log[10]~ "adjusted p value"),
                                  legendLabels = c('NS', "FC", "p val", "FC + p val"),
                             #                       'p-value', "p-value and avg logFC"),
                                            legendPosition = "right")
                pdf(file.path(project_path, "de", glue::glue("volcano_{de_list}.pdf")), width = width, height = height)
        print(volcano)
        dev.off()
        return(volcano)
    }
}
}
##############################################################################
##############################################################################

enrichrFun <- function(de_list, threshold, sheet) {
sheets <- readxl::excel_sheets(file.path(project_path, "de", glue::glue("de_{de_list}.xlsx")))
input <- readxl::read_excel(file.path(project_path, "de", glue::glue("de_{de_list}.xlsx")), sheet = sheet)
input_pos <- filter(input, avg_log2FC > threshold)
input_neg <- filter(input, avg_log2FC < -threshold)
input_merge <- list(pos = input_pos, neg = input_neg)
result <- list()
for (i in names(input_merge)) {
    result[[i]] <- enrichR::enrichr(input_merge[[i]]$gene, dbs)
    for (j in seq_along(result[[i]])){
        result[[i]][[j]][c("Old.Adjusted.P.value", "Old.P.value")] <- NULL
    }#remove old p value
    names(result[[i]])[1] <- "TF_Pertubations"
    names(result[[i]])[3] <- "Enrichr_Submissions_TF"
    writexl::write_xlsx(result[[i]], file.path(project_path, "enrichr", glue::glue("enrichr_{de_list}_{sheet}_{i}.xlsx")))
}
}
