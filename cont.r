
cont <- function(x){
	co=0
	n <- length(x)
	for (i in 1:n){
		if(is.na(x[i])==TRUE) co=co+1 }
	pos=NA
	w=1
	for (i in 1:n){
		if(is.na(x[i])==TRUE) pos[w]=i
		w=w+1 	}
	pos <- pos[is.na(pos)==FALSE]
	list(cuenta=co,pos=pos)
}