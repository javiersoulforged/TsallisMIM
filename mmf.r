
mmf <- function(X,s=1,mm=TRUE,p=0.75){

#mm: TRUE: media movil
#    FALSE: media fija

	m=dim(X)[2]
	t=dim(X)[1]
	n=dim(X)[1]/s
	u=(1-p)*s

	if(mm==TRUE){
		datos.a=matrix(NA,t,m)
		for(i in (s+1):t) for(j in 1:m){ 
			aux=X[(i-s+1):i,j]
			if(cont(aux)$cuenta<=u){ 
				datos.a[i,j]=mean(aux[is.na(aux)==FALSE])
			}
		}
	}

	if(mm==FALSE){
		datos.a = matrix(NA,n,m)
		for(i in 1:n) for(j in 1:m){ 
			aux=X[(s*(i-1)+1):(s*i),j]
			if(cont(aux)$cuenta<=u){ 
				datos.a[i,j]=mean(aux[is.na(aux)==FALSE])
			}
		}				
	}
	datos.a
}
