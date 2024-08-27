
Global.Renyi.alpha <- function(A,alpha){
	ev <- eigen(A)
	E=ev$values
	k = length(E)
	ST = sum(E)
	R = 1 + log(sum((E/ST)^alpha)) / ((alpha-1)*log(k))
	return(R)
}