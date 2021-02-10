#' Make raster plot of an OutbreakLandscape file
#'
#' @param d Directory containing data.
#' @param V Variable of interest.
#' @param v vector of values. 
#' @param s landscape size. Default=200.
#' @param p Species. Default=SBW.
#' @param t Optional. Title for the plot. Default: False.
#' @param timestep Optional. Timestep of the data file to plot. Default: 50
#' @return Raster image of OutbreakLandscape file.
#' @export
makeimage <- function(d, V, v, s=200, p='SBW', t=F, timestep=50){
	require(ggplot2)
	require(reshape2)
	Out <- list()
	for(i in v){
		tmp <- read.csv(paste0(d, 'OutbreakLandscape', V, i, 'run4.csv', p, '.csv'))
		landscapesquare <- square(tmp[timestep,])
		landscapesmall <- as.matrix(landscapesquare[1:s,1:s])
		colnames(landscapesmall) <- paste0("P", c(1:s))
		rownames(landscapesmall) <- paste0("p", c(1:s))
		landscapeplot <- as.data.frame(reshape2::melt(landscapesmall))
		landscapeplot[,3] <- landscapeplot[,3]+1
		if(t == T){ 
		Out[[which(v == i)]] <- ggplot(landscapeplot, aes(x = Var2, y = Var1)) +   geom_raster(aes(fill=log(value)))+ scale_fill_gradient(low="white", high="green") + labs(x="", y="") + theme_minimal()+ theme(legend.position='none', axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank()) + coord_equal() + ggtitle(paste(V, i)) 
		} else if(is.character(t)){
		Out[[which(v == i)]] <- ggplot(landscapeplot, aes(x = Var2, y = Var1)) +   geom_raster(aes(fill=log(value)))+ scale_fill_gradient(low="white", high="green") + labs(x="", y="") + theme_minimal()+ theme(legend.position='none', axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank()) + ggtitle(t) + coord_equal()
		} else { 
		Out[[which(v == i)]] <- ggplot(landscapeplot, aes(x = Var2, y = Var1)) +   geom_raster(aes(fill=log(value)))+ scale_fill_gradient(low="white", high="green") + labs(x="", y="") + theme_minimal()+ theme(legend.position='none', axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank()) + coord_equal()
		}
	}
	return(Out)
}

