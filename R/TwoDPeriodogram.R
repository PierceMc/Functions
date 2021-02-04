#' @export
ThreePlots <- function(d,r, m, RN=F){
	SBWabundance <- sum(d[r,])
	plots <- list()
	snap <- as.matrix(square(d[r,]))
	row.names(snap) <- c(1:nrow(snap))
	N=nrow(snap)
	#image(log(snap), col=hcl.colors(5, "Grays"))
	msnap <- reshape2::melt(snap)
	plots[[1]] <- ggplot(data=as.data.frame(msnap)) + geom_tile(aes(Var2,Var1, fill=value)) + scale_fill_distiller(palette = "YlGnBu", direction = 1) + theme_bw()  + theme(legend.position="none") + xlab("") + ylab("")

	p=TwoDimPeriodogram(snap)
	x=dim(p)[1]
	y=dim(p)[2]
	xover2=x/2
	yover2=y/2

	preturn=fftshift(p)
	mp <- reshape2::melt(fftshift(p))
	mp[,1] <- rep((-xover2):(xover2-1),x)
	mp[,2] <- rep((-yover2):(yover2-1), each=y )
	plots[[2]] <- ggplot(data=as.data.frame( mp )) + geom_tile(aes(Var2,Var1, fill=value)) + scale_fill_distiller(palette = "YlGnBu", direction = 1) + theme_bw() + ggtitle(paste(r, ", SBW abundance =", SBWabundance))+ theme(legend.position="none") + xlab("Cycles/Cell") + ylab("Cycles/Cell")
	 #plot(mp, col=hcl.colors(100, "Grays"))+ theme(legend.position="none")

	p=continue(p, N)
	pr=radial.psd(snap, plot=F)
	p=list("X"=pr$wavenumber, "Y"=pr$r_spectrum)
	plots[[3]] <- ggplot(as.data.frame(p)) + geom_line(aes(x=p[["X"]], y=p[["Y"]])) + scale_fill_distiller(palette = "YlGnBu", direction = 1) + theme_bw()+ theme(legend.position="none") + ylab("Power") + xlab("Cycles/Cell") #+ coord_cartesian(ylim=c(0,2.5e23))
	#plot(p[["X"]], p[["Y"]], col=hcl.colors(5, "Grays"), type='l')

#	squared <- square(d)
#	hist0 <- SpatialHistogramOnSquare(squared, "qplot")
#	plots[[4]] <- qplot(histo[,2], geom="histogram", binwidth = 10)

	if(RN==F){
		return(plots)
	} else if(RN=="2D"){
		return(snap)
	}else {
		return(p)
	}
}


#' @export
fftshift <- function(input_matrix, dim = -1) {

    rows <- dim(input_matrix)[1]    
    cols <- dim(input_matrix)[2]    

    swap_up_down <- function(input_matrix) {
        rows_half <- ceiling(rows/2)
        return(rbind(input_matrix[((rows_half+1):rows), (1:cols)], input_matrix[(1:rows_half), (1:cols)]))
    }

    swap_left_right <- function(input_matrix) {
        cols_half <- ceiling(cols/2)
        return(cbind(input_matrix[1:rows, ((cols_half+1):cols)], input_matrix[1:rows, 1:cols_half]))
    }

    if (dim == -1) {
        input_matrix <- swap_up_down(input_matrix)
        return(swap_left_right(input_matrix))
    }
    else if (dim == 1) {
        return(swap_up_down(input_matrix))
    }
    else if (dim == 2) {
        return(swap_left_right(input_matrix))
    }
    else {
        stop("Invalid dimension parameter")
    }
}

#' @export
ifftshift <- function(input_matrix, dim = -1) {

    rows <- dim(input_matrix)[1]    
    cols <- dim(input_matrix)[2]    

    swap_up_down <- function(input_matrix) {
        rows_half <- floor(rows/2)
        return(rbind(input_matrix[((rows_half+1):rows), (1:cols)], input_matrix[(1:rows_half), (1:cols)]))
    }

    swap_left_right <- function(input_matrix) {
        cols_half <- floor(cols/2)
        return(cbind(input_matrix[1:rows, ((cols_half+1):cols)], input_matrix[1:rows, 1:cols_half]))
    }

    if (dim == -1) {
        input_matrix <- swap_left_right(input_matrix)
        return(swap_up_down(input_matrix))
    }
    else if (dim == 1) {
        return(swap_up_down(input_matrix))
    }
    else if (dim == 2) {
        return(swap_left_right(input_matrix))
    }
    else {
        stop("Invalid dimension parameter")
    }
}
#sbw=d

#' @export
TwoDimPeriodogram<- function(sbw, spli=2, val="", averageline=0){
	sbw = sbw-mean(sbw)
	N <- nrow(sbw)
	xdft <- (1/(N^(1/2)))*fft2d(sbw)
	p <- abs(xdft)^2
	return(p)
}
#' @export
continue <- function(p, N, val = ""){
	Fs <- 1 #Sampling Rate
	spli=2
	p <- p[1:(N/spli+1)]
	p[2:(N/spli)] = spli*p[2:(N/spli)];
	freq<-seq(0,(Fs/spli),by=(Fs/N))
	return(list("X"=freq[2:(N/spli)], "Y"=p[2:(N/spli)]) )
}

#' @export
SpatialPeriodogram <- function(d, r){
squareexls <- as.matrix(square(d[r,]))
p=TwoDimPeriodogram(squareexls)
x=dim(p)[1]
y=dim(p)[2]
xover2=x/2
yover2=y/2
mp <- reshape2::melt(fftshift(p))
mp[,1] <- rep((-xover2):(xover2-1),x)
mp[,2] <- rep((-yover2):(yover2-1), each=y )
Plot <- ggplot(data=mp) + geom_tile(aes(Var2,Var1, fill=value)) + scale_fill_distiller(palette = "YlGnBu", direction = 1) + theme_bw() + ggtitle("")+ theme(legend.position="none") +ylab("Distance (Cells)") + xlab("Distance (Cells)")
return(Plot)
}



