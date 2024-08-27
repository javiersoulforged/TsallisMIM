#install.packages("plot3D")
library(entropy)
library(ggplot2)
library(plotrix)
library(plot3D)

source("Global.Shannon.r")
source("Renyi.empirical.r")
source("Renyi.mi.r")
source("Global.Renyi.r")
source("Renyi.alpha.empirical.r")
source("Renyi.alpha.mi.r")
source("Global.Renyi.alpha.r")
source("cont.r")
source("mmf.r")

Tsallis <- function(R,q) (exp((1-q)*R) - 1) / (1-q) 

Global.Tsallis <- function(Gr,q,k){
	d = log(k)
	A = 1 - exp((1-q)*(1-Gr)*d)
	B = (1-q)*d
	Gt = 1 + A/B
	return(Gt)
}

datosO3 <- read.table("data-O3.txt",header=TRUE)
names(datosO3)
attach(datosO3)

dim(datosO3)
d03=datosO3[20060300<datosO3$date & datosO3$date<20060400,] #March,2006
#d03=datosO3
names(d03)
data1=d03[,4:10]
names(data1)

s <- seq(4,48,1)
k=dim(data1)[2] # Number of Stations

H <- matrix(NA,length(s),3)

for(p in 1:length(s)){

	daux=mmf(data1,s=s[p],mm=TRUE)[-c(1:s[p]),]

	X <- daux
		
	A <- B <- matrix(NA,k,k)

for(i in 1:k) {
	nbx = round(diff(range(X[,i]))/(2*IQR(X[,i], na.rm = TRUE)/length(X[,i])^(1/3)),0)
	y1d = discretize(X[,i], numBins=nbx)
	A[i,i] = entropy(y1d)
	B[i,i] = Renyi.empirical(y1d)

	for(j in 1:k){
	if(i > j){
		nby = round(diff(range(X[,j]))/(2*IQR(X[,j], na.rm = TRUE)/length(X[,j])^(1/3)),0)
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

round(H, 4)

plot(s, H[,1], col="blue", ylim=c(0.15,1.25), xlab=expression(s), 
     ylab="Global measures", main="", type="l", lwd=2, cex=1.5, 
	cex.lab=1.5, cex.axis=1.5)
lines(s, H[,2], col="red", pch=2, lwd=2)
lines(s, H[,3], col="violet", pch=2, lwd=2)

#lo <- loess(H[,1]~s, na.action = na.exclude)
#lines(s,predict(lo), col='blue', lwd=2)
#lo <- loess(H[,2]~s)
#lines(s,predict(lo), col='red', lwd=2)

arrows(4, 1.06, 8, 1.175, col = "black",lwd=2)
text(x=11, y=1.18, label=paste("s=4"), cex=1.5)
arrows(21, 0.6, 23, 0.45, col = "black",lwd=2)
text(x=30, y=1.05, label=paste("s=28"), cex=1.5)
arrows(28, 0.85, 30, 1, col = "black",lwd=2)
text(x=23, y=0.4, label=paste("s=21"), cex=1.5)

#arrows(22, 0.197, 27, 0.025, col = "black",lwd=2)
#  text(x=32, y=0.015, label=paste("s=22"), cex=1.5)

legend("topright", c("Shannon","Renyi","Tsallis"), 
col=c("blue","red","violet"), lwd=2, bty="n")

max(H[,1])
s[5]

max(H[,2])
s[1]

min(H[,1])
s[19]

min(H[,2])
s[18]

### Completa data

N=dim(data1)[1]
k=dim(data1)[2] # Number of Stations
data2 <- data1

for(i in 1:N) for(j in 1:k){
	if(is.na(data1[i,j])) data2[i,j] = round(mean(data1[c(i-1,i+1),j]),0)
}


### Usando ventanas

s = 24
N=dim(data2)[1]
f=N/s
k=dim(data2)[2] # Number of Stations
alpha <- seq(1,15,1)

R <- T <- matrix(NA,f,length(alpha))

for(p in 1:f) for(a in 1:length(alpha)) {
	X <- data2[(s*(p-1)+1):(s*p),]
	A <- matrix(NA,k,k)

for(i in 1:k) {
	nbx = round(diff(range(X[,i], na.rm = TRUE))/(2*IQR(X[,i], na.rm = TRUE)/length(X[,i])^(1/3)),0)
	y1d = discretize(X[,i], numBins=nbx)
	if(alpha[a] == 1) A[i,i] = entropy(y1d)
	if(alpha[a] > 1) A[i,i] = Renyi.alpha.empirical(y1d, alpha[a])

	for(j in 1:k){
	if(i > j){
		nby = round(diff(range(X[,j], na.rm = TRUE))/(2*IQR(X[,j], na.rm = TRUE)/length(X[,j])^(1/3)),0)
		y2d = discretize2d(X[,i], X[,j], numBins1=nbx, numBins2=nby)
		if(alpha[a] == 1) A[i,j] = A[j,i] = mi.empirical(y2d)
		if(alpha[a] > 1) A[i,j] = A[j,i] = Renyi.alpha.mi(y2d, alpha[a])
	}
	}
}

if(alpha[a] == 1) R[p,a] = T[p,a] = Global.Shannon(A)
if(alpha[a] > 1){ 
	R[p,a] = Global.Renyi.alpha(A,alpha[a])
	T[p,a] = Global.Tsallis(R[p,a],alpha[a],k)
}

}

R
dim(R)
Global.Tsallis(2,4,k)

R[6,1] = mean(c(R[5,1],R[7,1]))
R[20,1] = mean(c(R[19,1],R[21,1]))
image2D(z=R, y=alpha, x=c(1:31), ylab=expression(alpha), xlab="days")

T[6,1] = mean(c(T[5,1],T[7,1]))
T[20,1] = mean(c(T[19,1],T[21,1]))
image2D(z=log(T), y=alpha, x=c(1:31), ylab=expression(q), xlab="days")

max(T)

