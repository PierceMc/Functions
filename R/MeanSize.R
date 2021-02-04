#' Generate summary statistics for OutbreakLandscape file
#'
#' @param  spacedata OutbreakLandscape file
#' @return Summary Statistics table
#' @export
MeanSize <- function(spacedata, s=0){
	cluster = makeCluster(3, type = "SOCK")
	registerDoSNOW(cluster)
	Output <- foreach(i=c((nrow(spacedata)-100):nrow(spacedata)), .combine='rbind') %dopar% {
		tmpsquare <- square(spacedata[i,])
		if(s==0){
			tmphosen <- hosen(tmpsquare)
		} else {
			tmphosen <- hosen(tmpsquare, s)
		}
		tmpsizes <- sizes(tmphosen)
		c(mean(tmpsizes[,2]), max(tmpsizes[,2]), sd(tmpsizes[,2]), length(which(tmpsizes[,2] > 2)))
	}
	stopCluster(cluster)
	colnames(Output) <- c('Size', 'Max', 'SD', 'No.')
	return(Output)
}

