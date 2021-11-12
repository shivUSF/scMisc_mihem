#' @title nice Seurat feature plot
#' @description create and save a nice Seurat feature plot in folder `featureplot`
#' @param object Seurat object
#' @param par column name in markers.xlsx
#' @param width width of output plot (default: 16)
#' @param height height of output plot (default: length of genes divided by two)
#' @return save feature plot to folder `featureplot`
#' @examples \dontrun{fPlot(sc_merge, par = "main", width = 5, height = 5)}
#' @export

fPlot <- function(object, par, width = 16, height = ceiling(length(genes)/2)) {
    dir.create(here::here("featureplot"), showWarnings = FALSE)
    stopifnot(file.exists("markers.csv"))
    markers <- readr::read_csv(here::here("markers.csv")) |>
    as.list(markers) |>
    lapply(function(x) x[!is.na(x)])
object_parse <- deparse(substitute(object))
genes <- markers[[par]]
fp <- Seurat::FeaturePlot(object = object, features = unique(genes), cols = viridis::viridis(100), reduction = "umap", pt.size = .1, order = TRUE, coord.fixed = TRUE) &
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.border = element_rect(color = "black", size = 1, fill = NA))
ggsave(here::here(project_path, "featureplot", glue::glue("{object_parse}_{par}.png")), width = width, height = height, limitsize = FALSE)
}
