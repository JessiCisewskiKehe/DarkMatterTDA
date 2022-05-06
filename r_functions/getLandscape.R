###Compute landscape function
#Define the standard landscape function kernel
kern_land_original <- function(x0, y0, x){
	###  kernel function (usual landscape)
	# (x0, y0) are the observations
	# x is the location to compute the kernel
	out <- ifelse(x-x0 < y0 - x, x-x0, y0-x)
	out[out<0] <- 0
	return(out)
}

#Compute the landscape function matrix
getLandscape <- function(subdiag, tseq){
	#subdiag = birth and death times (for given homology group)
	#tseq = grid on which to compute landscape
	kernel_matrix <- matrix(0, nrow = length(tseq), ncol = nrow(subdiag))
		for(ii in 1:nrow(subdiag)){
			kernel_matrix[,ii] <- kern_land_original(subdiag[ii, 1], subdiag[ii, 2], tseq)
	}
	lambda <- sapply(seq(along = tseq), FUN = function(i) {
        sort(kernel_matrix[i, ], decreasing = TRUE)
    })
	return(t(lambda))
}

