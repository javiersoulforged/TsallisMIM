
Global.Renyi <- function(A){
	ev <- eigen(A)
	E=ev$values
	k = length(E)
	ST = sum(E)
	R = 1 + log(sum((E/ST)^2)) / log(k)
	return(R)
}