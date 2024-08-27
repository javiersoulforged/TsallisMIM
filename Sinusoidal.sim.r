library(LongMemoryTS)
library(entropy)
library(ggplot2)
library(plotrix)
library(GA)

source("Global.Shannon.r")
source("Renyi.empirical.r")
source("Renyi.mi.r")
source("Global.Renyi.r")
source("Renyi.alpha.empirical.r")
source("Renyi.alpha.mi.r")
source("Global.Renyi.alpha.r")

Tsallis <- function(R,q) (exp((1-q)*R) - 1) / (1-q) 

Global.Tsallis <- function(Gr,q,k){
	d = log(k)
	A = 1 - exp((1-q)*(1-Gr)*d)
	B = (1-q)*d
	Gt = 1 + A/B
	return(Gt)
}

k <- seq(10,200,10)
L=100 # 100, 500, 1000
H <- matrix(NA,length(k),3)

for(p in 1:length(k)){

	X <- matrix(NA, L, k[p])

	for(i in 1:k[p]){
		T = 5*i
		eps <- rnorm(L)
		for(t in 1:L){
		beta = 10
		if(i > 20) beta = 0
		X[t,i] = beta*sin(2*pi*t/T) + eps[t]		
		}
	}

		
A <- B <- matrix(NA,k[p],k[p])

for(i in 1:k[p]) {
	nbx = round(diff(range(X[,i]))/(2*IQR(X[,i])/length(X[,i])^(1/3)),0)
	y1d = discretize(X[,i], numBins=nbx)
	A[i,i] = entropy(y1d)
	B[i,i] = Renyi.empirical(y1d)

	for(j in 1:k[p]){
	if(i > j){
		nby = round(diff(range(X[,j]))/(2*IQR(X[,j])/length(X[,j])^(1/3)),0)
		y2d = discretize2d(X[,i], X[,j], numBins1=nbx, numBins2=nby)
		A[i,j] = A[j,i] = mi.empirical(y2d)
		B[i,j] = B[j,i] = Renyi.mi(y2d)
	}
	}
}

E <- eigen(A)$values
ST = sum(E)
H[p,1] = 1 + sum(E*log(E/ST)/ST) / log(k[p])

E <- eigen(B)$values
ST = sum(E)
H[p,2] = 1 + log(sum((E/ST)^2)) / log(k[p])

H[p,3] = Global.Tsallis(H[p,2], 2, k[p])
}

plot(k, H[,1], col="blue", xlab=expression(k), ylim=c(0,0.9), pch=1,
     ylab="Global measures", main=paste("N=",L,sep=""))
points(k, H[,2], col="red", pch=2)
lo <- loess(H[,1]~k, na.action = na.exclude)
lines(k, predict(lo), col='blue', lwd=2)
lo <- loess(H[,2]~k)
lines(k, predict(lo), col='red', lwd=2)
points(k, H[,3], col="violet", pch=6)
lo <- loess(H[,3]~k, na.action = na.exclude)
lines(k, predict(lo), col='violet', lwd=2)

legend("right", c("Shannon", "Renyi", "Tsallis"), bty="n",
col=c("blue","red","violet"), pch=c(1,2,6), lty=c(NA,NA,NA), lwd=2)
