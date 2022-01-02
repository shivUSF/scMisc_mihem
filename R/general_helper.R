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

################################################################################
# list environment objects
################################################################################

#' @title pretty function to list environment objects
#' @description list objects from the environment, changes to human-readable units and arranges them by size
#' @return tibble with environment objects arranged by size
#' @export

lss <- function () {
data.frame('object' = ls(".GlobalEnv")) %>% 
      dplyr::mutate(size_unit = object %>%sapply(. %>% get() %>% object.size %>% format(., unit = 'auto')),
                    size = sapply(strsplit(size_unit, split = ' '), FUN = function(x) x[1]),
                    unit = factor(sapply(strsplit(size_unit, split = ' '), FUN = function(x) x[2]), levels = c('Gb', 'Mb', 'Kb', 'bytes'))) %>% 
      dplyr::arrange(unit, dplyr::desc(size)) %>%
      dplyr::select(-size_unit) %>%
    tibble::as_tibble()
}
