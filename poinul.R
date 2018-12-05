#N(0,1) vs. N(mu,1)
#with a N(0,1) prior

n=10^2
xobs=rnorm(n,mean=0,sd=1)
dnormxobs=dnorm(xobs)
meanobs=sum(xobs)/(n+1)
sdobs=1/sqrt(n+1)

#Dirichlet weights
a0=c(1/n,0.1,.2,.3,.4,.5,1.0)
#Gibbs sampler
T=10^5
alfs=matrix(.5,7,T)
thes=rep(meanobs,T)

for (i in 1:7){

  z=1+(runif(n)<1/(1+alfs[i,1]*dnorm(xobs,mean=thes[1])/((1-alfs[i,1])*dnormxobs)))
  for (t in 1:(T-1)){

  #Hasting independent proposal:
  proplik=curlik=sum(dnorm(xobs[z==2],log=TRUE))
  curlik=curlik+sum(dnorm(xobs[z==1],mean=thes[t],log=TRUE))
  prop=rnorm(1,meanobs,sdobs)
  proplik=proplik+sum(dnorm(xobs[z==1],mean=prop,log=TRUE))
  thes[t+1]=thes[t]
  if (log(runif(1))<proplik-curlik+dnorm(prop,log=TRUE)+dnorm(thes[t],meanobs,sdobs,log=TRUE)-
	dnorm(thes[t],log=TRUE)-dnorm(prop,meanobs,sdobs,log=TRUE)) thes[t+1]=prop
  z=1+(runif(n)<1/(1+alfs[i,t]*dnorm(xobs,mean=thes[t+1])/((1-alfs[i,t])*dnormxobs)))
  alfs[i,t+1]=1/(1+rgamma(1,a0[i]+sum((z==2)))/rgamma(1,a0[i]+sum((z==1))))
  }}

Alfs=alfs[,-1]

boxplot(t(alfs),names=c("1/n",".1",".2",".3",".4",".5","1"),col="wheat2",outline=FALSE)
par(new=TRUE)
plot(alfs[1,],col="grey10",type="l",ylim=c(0,1),axes=FALSE)
for (i in 2:6) lines(alfs[i,],col=paste("grey",15*i,sep=""))

#observing the evolución with sample size
T=10^5
n=10^3
sqrtt=sqrt(2)
robs=rnorm(n,mean=0)#(1-2*(runif(n)<.5))*rexp(n,rat=sqrtt) #rnorm(n)
size=c(5,10,50,100,500,n)
a0=1/(sqrt(size)) #rep(.1,6)
befs=matrix(.5,6,T)
meanbefs=rep(0,6)
meanobs=mean(robs)
dnormrobs=dnorm(robs)
thes=rep(meanobs,T)

for (i in 1:6){
  ni=size[i]
  xobs=robs[1:ni]
  meanobs=mean(xobs)
  dnormxobs=dnormrobs[1:ni]
  sdobs=1/sqrt(ni)

  if (i>1) befs[i,1]=meanbefs[i-1]
  z=1+(runif(ni)<1/(1+befs[i,1]*dnorm(xobs,mean=thes[1])/((1-befs[i,1])*dnormxobs)))

  for (t in 1:(T-1)){

  #Hasting independent proposal:
  curlik=sum(dnorm(xobs[z==2],log=TRUE))+sum(dnorm(xobs[z==1],mean=thes[t],log=TRUE))
  prop=rnorm(1,meanobs,sdobs)
  proplik=sum(dnorm(xobs[z==2],log=TRUE))+sum(dnorm(xobs[z==1],mean=prop,log=TRUE))
  thes[t+1]=thes[t]
  if (runif(1)<exp(proplik-curlik)*dnorm(prop)*dnorm(thes[t],meanobs,sdobs)/
	(dnorm(thes[t])*dnorm(prop,meanobs,sdobs))) thes[t+1]=prop
  z=1+(runif(ni)<1/(1+befs[i,t]*dnorm(xobs,mean=thes[t])/((1-befs[i,t])*dnormxobs)))
  meanbefs[i]=meanbefs[i]+(a0[i]+sum((z==1)))/(2*a0[i]+ni)
  befs[i,t+1]=1/(1+rgamma(1,a0[i]+sum((z==2)))/rgamma(1,a0[i]+sum((z==1))))
  }}

befs=befs[,-1]
meanbefs=meanbefs/(T-1)

boxplot(t(befs),names=c("5","10","50","100","500","1000"),col="wheat2",outline=FALSE)
par(new=TRUE)
plot(apply(befs,1,mean),col="tomato",axes=FALSE,ylim=c(0,1))
lines(apply(befs,1,median),col="sienna")
X11()
par(mfrow=c(6,1),mar=c(0,0,0,0))
for (i in 1:6) plot(befs[i,],type="l",col="sienna",axes=FALSE)

#off with Gibbs
T=10^6
n=10^3
sqrtt=sqrt(2)
robs=rnorm(n)#(1-2*(runif(n)<.5))*rexp(n,rat=sqrtt) #rnorm(n)
size=c(5,10,50,100,500,n)
a0=rep(.1,6) #1/((size)) #rep(.1,6)
meanbefs=rep(0,6)
meanobs=mean(robs)
dnormrobs=dnorm(robs)
thes=rep(meanobs,T)
befs=matrix(.5,6,T)

for (i in 1:6){
  ni=size[i]
  xobs=robs[1:ni]
  meanobs=mean(xobs)
  dnormxobs=dnormrobs[1:ni]
  sdobs=1/sqrt(ni)

  if (i>1){ befs[i,1]=befs[i-1,T];thes[1]=thes[T]}
  curlik=sum(log((1-befs[i,1])*dnormxobs+befs[i,1]*dnorm(xobs,mean=thes[1])))

  for (t in 1:(T-1)){

  #Hasting independent proposal:
  prop=rnorm(1,meanobs,sdobs)
  proplik=sum(log((1-befs[i,t])*dnormxobs+befs[i,t]*dnorm(xobs,mean=prop)))
  thes[t+1]=thes[t]
  if (runif(1)<exp(proplik-curlik)*dnorm(prop)*dnorm(thes[t],meanobs,sdobs)/
	(dnorm(thes[t])*dnorm(prop,meanobs,sdobs))){
     thes[t+1]=prop
     curlik=proplik}

  z=1+(runif(ni)<1/(1+befs[i,t]*dnorm(xobs,mean=thes[t])/((1-befs[i,t])*dnormxobs)))
  befs[i,t+1]=1/(1+rgamma(1,a0[i]+sum((z==2)))/rgamma(1,a0[i]+sum((z==1))))
  curlik=sum(log((1-befs[i,t+1])*dnormxobs+befs[i,t+1]*dnorm(xobs,mean=thes[t+1])))
  prop=rbeta(1,a0[i],a0[i])
  proplik=sum(log((1-prop)*dnormxobs+prop*dnorm(xobs,mean=thes[t+1])))
  if (runif(1)<exp(proplik-curlik)){
    befs[i,t+1]=prop
    curlik=proplik}
  }}

befs=befs[,-1]

boxplot(t(befs),names=c("5","10","50","100","500","1000"),col="wheat2",outline=FALSE)
X11(width=12,height=5)
par(mfrow=c(6,1),mar=c(0,0,0,0))
for (i in 1:6) plot(befs[i,],type="l",col="sienna",axes=FALSE)
