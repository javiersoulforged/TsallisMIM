Renyi.mi <- function(freqs2d){

  	freqs2d = as.matrix(freqs2d/sum(freqs2d)) # just to make sure ...

  	freqs.x = rowSums(freqs2d) # marginal frequencies of x
  	freqs.y = colSums(freqs2d) # marginal frequencies of y
  	freqs.null = freqs.x %o% freqs.y # independence null model

	R2 = log(sum( ifelse(freqs2d > 0, freqs2d^2/freqs.null, 0) ))

  	return(R2)	
}



