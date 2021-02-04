#' Make raster plot of an OutbreakLandscape file
#'
#' @param d Directory containing data
#' @param V Variable of interest
#' @param v value of interest
#' @param s landscape size. Default=200
#' @param p Species. Default=SBW
#' @return Histogram of patch sizes
#' @export
makeimage <- function(d, V, v, s=200, p='SBW'){
	Out <- list()
	for(i in v){
		tmp <- read.csv(paste0(d, 'OutbreakLandscape', V, i, 'run4.csv', p, '.csv'))
		landscapesquare <- square(tmp[100,])
		landscapesmall <- as.matrix(landscapesquare[1:s,1:s])
		colnames(landscapesmall) <- paste0("P", c(1:s))
		rownames(landscapesmall) <- paste0("p", c(1:s))
		landscapeplot <- as.data.frame(reshape2::melt(landscapesmall))
		Out[[which(v == i)]] <- ggplot(landscapeplot, aes(x = Var2, y = Var1)) +   geom_raster(aes(fill=value))+ scale_fill_gradient(low="white", high="black") + labs(x="", y="") + theme_bw()+ theme(legend.position='none', axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank()) + ggtitle(paste(V, i))

	}
	return(Out)
}
