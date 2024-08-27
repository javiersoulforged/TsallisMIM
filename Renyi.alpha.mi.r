Renyi.alpha.mi <- function(freqs2d,alpha){

  	freqs2d = as.matrix(freqs2d/sum(freqs2d)) # just to make sure ...

  	freqs.x = rowSums(freqs2d) # marginal frequencies of x
  	freqs.y = colSums(freqs2d) # marginal frequencies of y
  	freqs.null = freqs.x %o% freqs.y # independence null model

	R.alpha = log(sum( ifelse(freqs2d > 0, freqs2d^alpha/(freqs.null)^(alpha-1), 0) ))

  	return(R.alpha)	
}



