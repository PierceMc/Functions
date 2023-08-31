#' Load ggplot themes on library load
#'
#' @export
.onLoad <- function(libname = find.package("phd"), pkgname = "phd") {
	textsize=15
	theme_set(theme_classic())
	theme_update(legend.position='none', legend.title=element_text(size=textsize),legend.text=element_text(size=textsize),axis.title=element_text(size=textsize))
	ggplot <- function(...) ggplot2::ggplot(...) + scale_fill_brewer(palette="Dark2") + scale_color_brewer(palette="Dark2")
}
