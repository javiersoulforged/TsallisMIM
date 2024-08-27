
Global.Shannon <- function(A){
	ev <- eigen(A)
	E=abs(ev$values)
	k = length(E)
	ST = sum(E)
	H = 1 + sum(E*log(E/ST)/ST) / log(k)
	return(H)
}