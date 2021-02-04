
#' @export
singlecellplots <- function(tmpdata1, tmpdata2, tmpdata3, ymax="NT", ymin="NT", values=c(1, 3, 5), fun="mean"){

rowmeans1 = apply(tmpdata1, 1, mean)
rowsds1 = apply(tmpdata1, 1, sd)
plotdata1 =data.frame("Iteration"=c(31500:(31499+nrow(tmpdata1))), "Mean"=rowmeans1, "SD"=rowsds1)

rowmeans2 = apply(tmpdata2, 1, mean)
rowsds2 = apply(tmpdata2, 1, sd)
plotdata2 =data.frame("Iteration"=c(31500:(31499+nrow(tmpdata2))), "Mean"=rowmeans2, "SD"=rowsds2)

rowmeans3 = apply(tmpdata3, 1, mean)
rowsds3 = apply(tmpdata3, 1, sd)
plotdata3 =data.frame("Iteration"=c(31500:(31499+nrow(tmpdata3))), "Mean"=rowmeans3, "SD"=rowsds3)


if(fun == "mean"){

	if(ymax=="NT"){
		max1=max(rowmeans1)+max(rowsds1)
		max2=max(rowmeans2)+max(rowsds2)
		max3=max(rowmeans3)+max(rowsds3)
		ymax=max(max1, max2, max3)
	}
	if(ymin=="NT"){
		min1=min(rowmeans1)-max(rowsds1)
		min2=min(rowmeans2)-max(rowsds2)
		min3=min(rowmeans3)-max(rowsds3)
		ymin=min(min1, min2, min3)
	}


	runplot1 <- ggplot(data=plotdata1, aes(x=Iteration, y=Mean, ymin=Mean-SD, ymax=Mean+SD))+geom_ribbon()+geom_point()+coord_cartesian(y=c(ymin,ymax)) + ggtitle(values[1])
	runplot2 <- ggplot(data=plotdata2, aes(x=Iteration, y=Mean, ymin=Mean-SD, ymax=Mean+SD))+geom_ribbon()+geom_point()+coord_cartesian(y=c(ymin,ymax)) + ggtitle(values[2])
	runplot3 <- ggplot(data=plotdata3, aes(x=Iteration, y=Mean, ymin=Mean-SD, ymax=Mean+SD))+geom_ribbon()+geom_point()+coord_cartesian(y=c(ymin,ymax)) + ggtitle(values[3])

	return(list(runplot1, runplot2, runplot3))
} else if(fun == "sd"){

	if(ymax=="NT"){
		max1=max(rowsds1)
		max2=max(rowsds2)
		max3=max(rowsds3)
		ymax=max(max1, max2, max3)
	}
	if(ymin=="NT"){
		min1=min(rowsds1)
		min2=min(rowsds2)
		min3=min(rowsds3)
		ymin=min(min1, min2, min3)
	}
	runplot1 <- ggplot(data=plotdata1, aes(x=Iteration, y=SD))+geom_line()+coord_cartesian(y=c(ymin,ymax)) + ggtitle(values[1])
	runplot2 <- ggplot(data=plotdata2, aes(x=Iteration, y=SD))+geom_line()+coord_cartesian(y=c(ymin,ymax)) + ggtitle(values[2])
	runplot3 <- ggplot(data=plotdata3, aes(x=Iteration, y=SD))+geom_line()+coord_cartesian(y=c(ymin,ymax)) + ggtitle(values[3])

	return(list(runplot1, runplot2, runplot3))
}
}

#' @export
meancellsperiodogram <- function(x1, x2, x3, value=c(1,3,5)){

	rowmeans1 = apply(x1, 1, mean)
	periodogram1 <- periodogramplotdata(rowmeans1, 1, 2)
	plotdata1 =data.frame("Freq"=periodogram1[[1]], "Power"=periodogram1[[2]])

	rowmeans2 = apply(x2, 1, mean)
	periodogram2 <- periodogramplotdata(rowmeans2, 1, 2)
	plotdata2 =data.frame("Freq"=periodogram2[[1]], "Power"=periodogram2[[2]])
	
	rowmeans3 = apply(x3, 1, mean)
	periodogram3 <- periodogramplotdata(rowmeans3, 1, 2)
	plotdata3 =data.frame("Freq"=periodogram3[[1]], "Power"=periodogram3[[2]])

	ymax=max(max(periodogram1[[2]]), max(periodogram2[[2]]), max(periodogram3[[2]]))
	ymin=0

	runplot1 <- ggplot(data=plotdata1, aes(x=Freq, y=Power))+geom_line()+coord_cartesian(y=c(ymin,ymax)) + ggtitle(value[1]) + ggtitle(value[1])
	runplot2 <- ggplot(data=plotdata2, aes(x=Freq, y=Power))+geom_line()+coord_cartesian(y=c(ymin,ymax)) + ggtitle(value[2]) + ggtitle(value[2])
	runplot3 <- ggplot(data=plotdata3, aes(x=Freq, y=Power))+geom_line()+coord_cartesian(y=c(ymin,ymax)) + ggtitle(value[3]) + ggtitle(value[3])

	return(list(runplot1, runplot2, runplot3))

}


