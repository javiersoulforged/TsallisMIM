#install.packages("GA")
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

rho <- seq(0,0.99,0.01)
k = 6  # 2, 4, 6
L=500  # 100, 250, 500

d <- rep(0, 10) # c(0.4, 0.3, 0.2, 0.4, 0.3, 0.3, -0.7, -0.3, 0.3, -0.8)
phi = 0 # 0.4
THETA=0  # -0.2

H <- matrix(NA,length(rho),3)

for(p in 1:length(rho)){

sigma2 = rho[p] # correlacion
Sigma <- matrix(sigma2, k, k)
diag(Sigma) = 1
X <- VARFIMA.sim(phi=phi, THETA=THETA, Sigma=Sigma, d.vec=d[1:k], T=L)
		
A <- B <- matrix(NA,k,k)

for(i in 1:k) {
	nbx = round(diff(range(X[,i]))/(2*IQR(X[,i])/length(X[,i])^(1/3)),0)
	y1d = discretize(X[,i], numBins=nbx)
	A[i,i] = entropy(y1d)
	B[i,i] = Renyi.empirical(y1d)

	for(j in 1:k){
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
H[p,1] = 1 + sum(E*log(E/ST)/ST) / log(k)

E <- eigen(B)$values
ST = sum(E)
H[p,2] = 1 + log(sum((E/ST)^2)) / log(k)

H[p,3] = Global.Tsallis(H[p,2], 2, k)
}

plot(rho, H[,1], col="blue", xlab=expression(rho), ylim=c(0,1), pch=1,
     ylab="Global measures", main=paste("k=",k ,", N=",L,sep=""))
points(rho, H[,2], col="red", pch=2)
lo <- loess(H[,1]~rho, na.action = na.exclude)
lines(rho, predict(lo), col='blue', lwd=2)
lo <- loess(H[,2]~rho)
lines(rho, predict(lo), col='red', lwd=2)
points(rho, H[,3], col="violet", pch=6)
lo <- loess(H[,3]~rho, na.action = na.exclude)
lines(rho, predict(lo), col='violet', lwd=2)

legend("topleft", c("Shannon", "Renyi", "Tsallis"),
col=c("blue","red","violet"), pch=c(1,2,6), lty=c(NA,NA,NA), lwd=2, 
bty="n")


