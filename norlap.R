#Laplace versus Normal
n=10^3
sqrtt=sqrt(2)
xobs=(1-2*(runif(n)<.5))*rexp(n,rat=sqrtt) #rnorm(n)
#xobs=rnorm(n,sd=1)
#zobs=(1-2*(runif(n)<.5))*rexp(n,rat=sqrtt)
#xobs=(runif(n)<.7)*(yobs*zobs)-zobs
meanobs=mean(xobs)
sdobs=1/sqrt(n)

a0=c(0.1,.2,.3,.4,.5,1.0)
#Gibbs sampler
T=10^4
alfs=matrix(.3,6,T)
thes=rep(meanobs,T)

for (i in 1:6){
  for (t in 1:(T-1)){

  z=1+(runif(n)<1/(1+alfs[i,t]*dnorm(xobs,thes[t])/((1-alfs[i,t])*.5*dexp(abs(xobs-thes[t]),rat=sqrtt))))
  #Hasting independent proposal:
  curlik=sum(dnorm(xobs[z==1],thes[t]),log=TRUE)+sum(dexp(abs(xobs[z==2]-thes[t]),rat=sqrtt,log=TRUE))
  prop=rnorm(1,meanobs,sdobs)
  proplik=sum(dnorm(xobs[z==1],prop),log=TRUE)+sum(dexp(abs(xobs[z==2]-prop),rat=sqrtt,log=TRUE))
  thes[t+1]=thes[t]
  if (runif(1)<exp(proplik-curlik)*dnorm(thes[t],meanobs,sdobs)/dnorm(prop,meanobs,sdobs)) thes[t+1]=prop
  alfs[i,t+1]=1/(1+rgamma(1,a0[i]+sum((z==2)))/rgamma(1,a0[i]+sum((z==1))))
  }}

alfs=alfs[,-1]

boxplot(t(alfs),names=c(".1",".2",".3",".4",".5","1"),col="wheat2",outline=FALSE)

T=10^5
n=10^3
sqrtt=sqrt(2)
robs=rnorm(n)#(1-2*(runif(n)<.5))*rexp(n,rat=sqrtt) #rnorm(n)
size=c(5,10,50,100,500,1000)
a0=rep(.5,6)
befs=matrix(.3,6,T)
meanobs=mean(robs)
thes=rep(meanobs,T)

for (i in 1:6){
  xobs=robs[1:size[i]]
  meanobs=mean(xobs)
  sdobs=1/sqrt(size[i])

  for (t in 1:(T-1)){

  z=1+(runif(size[i])<1/(1+befs[i,t]*dnorm(xobs,thes[t])/((1-befs[i,t])*.5*dexp(abs(xobs-thes[t]),rat=sqrtt))))
  #Hasting independent proposal:
  curlik=sum(dnorm(xobs[z==1],thes[t]),log=TRUE)+sum(dexp(abs(xobs[z==2]-thes[t]),rat=sqrtt,log=TRUE))
  prop=rnorm(1,meanobs,sdobs)
  proplik=sum(dnorm(xobs[z==1],prop),log=TRUE)+sum(dexp(abs(xobs[z==2]-prop),rat=sqrtt,log=TRUE))
  thes[t+1]=thes[t]
  if (runif(1)<exp(proplik-curlik)*dnorm(thes[t],meanobs,sdobs)/dnorm(prop,meanobs,sdobs)) thes[t+1]=prop
  befs[i,t+1]=1/(1+rgamma(1,a0[i]+sum((z==2)))/rgamma(1,a0[i]+sum((z==1))))
  }}

befs=befs[,-1]

boxplot(t(befs),names=c("5","10","50","100","500","1000"),col="wheat2",outline=FALSE)

