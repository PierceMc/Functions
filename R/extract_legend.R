#' Extract legends from ggplot
#'
#' @param my_ggp GGPlot object
#' @return grob containing the legend from my_gpp
#' @export
extract_legend <- function(my_ggp) {
  step1 <- ggplot_gtable(ggplot_build(my_ggp))
  step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
  step3 <- step1$grobs[[step2]]
  return(step3)
}
