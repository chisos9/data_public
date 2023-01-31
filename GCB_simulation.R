rm(list=ls())
cat("\014")

set.seed(139)

cT<-62
ax<-(1959:(1959+cT-1))
beta1<-0.01
beta2<-0.01
d<-0.1
alpha1<- -4.13
alpha2<- -5.11
sigmaC<-sqrt(1.42)
sigmaLND<-sqrt(0.42)
sigmaOCN<-sqrt(0.015)
sigmaE<-sqrt(0.01)

E<-matrix(0,cT,1)
S_LND<-matrix(0,cT,1)
S_OCN<-matrix(0,cT,1)
C<-matrix(0,cT,1)

E[1]<-4.249917
C[1]<-2.127*315.39
S_LND[1]<-alpha1+beta1*C[1]
S_OCN[1]<-alpha2+beta2*C[1]

eps1<-sigmaLND*rnorm(cT)
eps2<-sigmaOCN*rnorm(cT)
eps3<-sigmaE*rnorm(cT)
eps4<-sigmaC*rnorm(cT)

c<-1+beta1+beta2

for (t in 2:cT){
  E[t]<-E[t-1]+d+eps3[t]
  S_LND[t]<-alpha1*(1-beta1/c)-alpha2*beta1/c+d*beta1/c + beta1/c*(E[t-1]+C[t-1]) + (1-beta1/c)*eps1[t]+beta1/c*(eps3[t]+eps4[t]-eps2[t])
  S_OCN[t]<-alpha2*(1-beta2/c)-alpha2*beta2/c+d*beta2/c + beta2/c*(E[t-1]+C[t-1]) + (1-beta2/c)*eps2[t]+beta2/c*(eps3[t]+eps4[t]-eps1[t])
  C[t]<-(d-alpha1-alpha2)/c + 1/c*(E[t-1]+C[t-1]) + 1/c*(eps3[t]+eps4[t]-eps1[t]-eps2[t])
}

par(mfrow=c(2,2))
plot(ax,E,type='l',xlab='Year')
plot(ax,S_LND,type='l',xlab='Year')
plot(ax,S_OCN,type='l',xlab='Year')
plot(ax[2:cT],diff(C),type='l',xlab='Year',ylab='G_ATM')

