
Tsallis.poisson <- function(lambda,q,x.lim=100){
	#px <- dpois(1:x.lim, lambda, log = FALSE)
	px <- exp(-lambda)*lambda^(1:x.lim) / factorial(1:x.lim)
      if(q == 1) T = -sum(px*log(px))
	if(q != 1) T = (1 - sum(px^q)) / (q - 1)
	return(T)
}

lambda = 15
L = 5
q <- c(seq(0.1,0.9,0.1),seq(1.1,L,0.1)) 
T <- rep(NA,length(q))
for(i in 1:length(q)) T[i] = Tsallis.poisson(lambda,q[i])

plot(q,T,type="l",lwd=2,ylim=c(0,25),ylab="Tsallis entropy")
lines(q,T,col="blue",lwd=2)
lines(q,T,col="red",lwd=2)
legend("topright",c(expression(gamma~"= 5"),expression(gamma~"= 10"),
expression(gamma~"= 15")), col=c("black","blue","red"), lwd=2, bty="n")
