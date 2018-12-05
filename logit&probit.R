# Logistic & probit with different parameters :
install.packages("mvtnorm") # For the multivariate normal
require(mvtnorm)
install.packages("corpcor") # For the positive definite matrix
require(corpcor)

# Analyse :
install.packages("MCMCpack")
 require(MCMCpack)

install.packages("coda")
 require(coda)

# Data package :
install.packages("MASS")
 require(MASS)

data(Pima.tr)
Nsim=15000
type=Pima.tr$type
 bmi=Pima.tr$bmi
X=cbind(1, bmi)
y=type
y=as.numeric(y)
y[y==1]=0
y[y==2]=1

c1= dim(X)[1]
c2=c1

flogis=function(theta1, y=y, X=X){
    prob=plogis(X%*%matrix(theta1,ncol=1))
    plgs=(prob^y) * (1-prob)^(1-y)
return(plgs)
}

fprobit=function(theta2, y=y, X=X){
  prob=pnorm(X%*%matrix(theta2,ncol=1))
  prbit=(prob^y) * (1-prob)^(1-y)
return(prbit)
}

# log function :
lflogis=function(theta1, y=y, X=X){
    prob=plogis(X%*%matrix(theta1,ncol=1))
prob[prob==1,]=0.9999999999
logprob= sum((y*log(prob)) + ((1-y) *log((1-prob))))
    return(logprob)
}

lfprobit=function(theta2, y=y, X=X){
  prob=pnorm(X%*%matrix(theta2,ncol=1))
prob[prob==1,]=0.9999999999
logprob= sum((y*log(prob)) + ((1-y) *log((1-prob))))
    return(logprob)
}
lpost_theta1=function(theta1, y, X, c1){
#c1= dim(X)[1]
if(is.vector(X)!= TRUE){
q=dim(X)[2]}else{
q=length(X)}

mean=rep(0,q)
txx=t(X)%*%X
stxx=solve(txx)
sigma=c1*stxx 
dnor= dmvnorm(theta1, mean, sigma, log=TRUE)
return(dnor+ lflogis(theta1, y=y, X=X))
}

lpost_theta2=function(theta2, y, X,c2){
#c2= dim(X)[1]
if(is.vector(X)!= TRUE){
q=dim(X)[2]}else{
q=length(X)}

mean=rep(0,q)
txx=t(X)%*%X
stxx=solve(txx)
sigma=c2*stxx 
dnor= dmvnorm(theta2, mean, sigma, log=TRUE)
return(dnor+ lfprobit(theta2, y=y, X=X))
}
# Gibbs parameters :Correct
GSmix_theta5=function(Nsim, e, tau, y, X){
k=length(e)
n=length(y)
if(is.vector(X)!= TRUE){
q=dim(X)[2]}else{
q=length(X)}

theta1=matrix(0, nrow=Nsim, ncol=q)
theta2=matrix(0, nrow=Nsim, ncol=q)
p=matrix(0, nrow=Nsim, ncol=k)
p[1,]= p0
mle1 = summary(glm(y~X, family=binomial(link="logit")))
  beta1= mle1$coefficients[,1]
  sigma.asymp1= mle1$cov.unscaled
mle2 = summary(glm(y~X, family=binomial(link="probit")))
  beta2 = mle2$coefficients[,1]
  sigma.asymp2= mle2$cov.unscaled
  theta1[1,]= beta1
  theta2[1,]= beta2

for(t in 2:Nsim){

z=vector("numeric", length=n)
for(l in 1:n){
z[l]=1+(runif(1)<1/(1+p[t-1,1]* flogis(theta1[t-1,], y[l], X[l,])/((1-p[t-1,1])*fprobit(theta2[t-1,], y[l], X[l,]))))
}

  #pz=p[t-1,]
  #       z=vector("numeric", length=n)
     #    probs=matrix(0, nrow=k, ncol=n)      
#probs[1,]=pz[1]*flogis(theta1[t-1,], y=y, X=X)
#probs[2,]=pz[2]*fprobit(theta2[t-1,], y=y, X=X)
#for(i in 1:n){
#                   z[i]=which.max(probs[,i])
#}

L=vector("numeric", length=k)
        for(I in 1:k){
L[I]=length(z[z==I])
}

alphap=L+e
p[t, ]=rdirichlet(1, alphap)

           xbar1=X[z==1,]
xbar1=matrix(xbar1, ncol=q)
           xx1=t(xbar1)%*%xbar1
           xbar2=X[z==2,]
xbar2=matrix(xbar2, ncol=q)
           xx2=t(xbar2)%*%xbar2
           yz1=y[z==1]
           yz2=y[z==2]
if(is.positive.definite(xx1)!=FALSE){
                            #proposal1=rmvnorm(1, theta1[t-1,], tau*sigma.asymp1)
proposal1=rmvnorm(1, beta1, tau*sigma.asymp1)

    logalpha1=min(0, lpost_theta1(proposal1, yz1, xbar1, c1) - lpost_theta1(theta1[t-1,], yz1, xbar1, c1))
    if(log(runif(1))<logalpha1){ theta1[t, ]= proposal1} else{theta1[t,]=theta1[t-1,]}
}else{theta1[t,]=theta1[t-1,]}

if(is.positive.definite(xx2)!=FALSE){
#proposal2=rmvnorm(1, theta2[t-1,], tau*sigma.asymp2)
proposal2=rmvnorm(1, beta2, tau*sigma.asymp2)
    logalpha2=min(0, lpost_theta2(proposal2, yz2, xbar2, c2) - lpost_theta2(theta2[t-1,], yz2, xbar2, c2))
    if(log(runif(1))<logalpha2){ theta2[t, ]= proposal2} else{theta2[t,]=theta2[t-1,]}
}else{theta2[t,]=theta2[t-1,]}
}
            TP=cbind(theta1, theta2, p)

return(TP)   
}

# Convergent MH :

GSmix_theta5=function(Nsim, e, tau, y, X){
k=length(e)
n=length(y)
if(is.vector(X)!= TRUE){
q=dim(X)[2]}else{
q=length(X)}
c=n
a0=e[1]

theta1=matrix(0, nrow=Nsim, ncol=q)
theta2=matrix(0, nrow=Nsim, ncol=q)
p=matrix(0, nrow=Nsim, ncol=k)
p[1,]= p0
mle1 = summary(glm(y~X, family=binomial(link="logit")))
  beta1= mle1$coefficients[,1]
  sigma.asymp1= mle1$cov.unscaled
mle2 = summary(glm(y~X, family=binomial(link="probit")))
  beta2 = mle2$coefficients[,1]
  sigma.asymp2= mle2$cov.unscaled

  theta1[1,]= beta1
  theta2[1,]= beta2

tx=t(X)%*%X
v1=solve(tx)
v2=c*v1
mean=rep(0, q)

for(t in 2 :Nsim){
#Hasting independent proposal theta1:
cur=(p[t-1,1]*flogis(theta1[t-1,], y=y, X=X))+(p[t-1,2]*fprobit(theta2[t-1,], y=y, X=X))
lcur=log(cur)

curlik=sum(lcur)+ log(dmvnorm(theta1[t-1,], mean, v2))
#proposal1=rmvnorm(1, theta1[t-1,], tau*sigma.asymp1)
proposal1=rmvnorm(1, beta1, tau*sigma.asymp1)

pr=(p[t-1,1]*flogis(proposal1, y=y, X=X))+(p[t-1,2]*fprobit(theta2[t-1,], y=y, X=X))
lpr=log(pr)
proplik= sum(lpr) + log(dmvnorm(c(proposal1), mean, v2))
theta1[t,]=theta1[t-1,]

alp= proplik-curlik

alpha1=min(0, alp)

if (log(runif(1))<alpha1){
     theta1[t,]=proposal1}
#Hasting independent proposal theta2:
curt=(p[t-1,1]*flogis(theta1[t,], y=y, X=X))+(p[t-1,2]*fprobit(theta2[t-1,], y=y, X=X))
lcurt=log(curt) 
curlikt=sum(lcurt)+ log(dmvnorm(theta2[t-1,], mean, v2))
#proposal2=rmvnorm(1, theta2[t-1,], tau*sigma.asymp2)
proposal2=rmvnorm(1, beta2, tau*sigma.asymp2)

prt=(p[t-1,1]*flogis(theta1[t,], y=y, X=X))+(p[t-1,2]*fprobit(proposal2, y=y, X=X))
lprt=log(prt)
proplikt= sum(lprt) + log(dmvnorm(c(proposal2), mean, v2))

theta2[t,]=theta2[t-1,]

alpt= proplikt-curlikt 

alpha2=min(0, alpt)

if (log(runif(1))<alpha2){
     theta2[t,]=proposal2}

z=1+(runif(n)<1/(1+p[t-1,1]*flogis(theta1[t-1,], y=y, X=X)/((1-p[t-1,1])*fprobit(theta2[t-1,], y=y, X=X))))
L=vector("numeric", length=k)
        for(I in 1:k){
L[I]=length(z[z==I])
}
alphap=L+e
p[t, ]=rdirichlet(1, alphap)

  #p[t,1]=1/(1+rgamma(1,a0+sum((z==2)))/rgamma(1,a0+sum((z==1))))
  #p[t,2]=1-p[t,1]

cur1= (p[t,1]*flogis(theta1[t,], y=y, X=X))+(p[t,2]*fprobit(theta2[t,], y=y, X=X))
lcur1=log(cur1)
curlik1= sum(lcur1) 

propos=rbeta(1,a0,a0)
pr1= (propos*flogis(theta1[t,], y=y, X=X))+((1-propos)*fprobit(theta2[t,], y=y, X=X))
lpr1=log(pr1)
proplik1= sum(lpr1) 
  if (runif(1)<exp(proplik1-curlik1)){
    p[t,1]=propos
p[t,2]=1- p[t,1]
 }
}
            TP=cbind(theta1, theta2, p)

return(TP)   
}

# Estimations of alpha for different values of tau :

tau=1
e=c(.1, .1)
estimats.e.1=GSmix_theta5(Nsim, e, tau, y, X)
e=c(.2, .2)
estimats.e.2=GSmix_theta5(Nsim, e, tau, y, X)
e=c(.3, .3)
estimats.e.3=GSmix_theta5(Nsim, e, tau, y, X)
e=c(.4, .4)
estimats.e.4=GSmix_theta5(Nsim, e, tau, y, X)
e=c(.5, .5)
estimats.e.5=GSmix_theta5(Nsim, e, tau, y, X)

# Another example with juste one explanatory variable like Pima data :
n=1000, 10000
I=rep(1, n)
X1=runif(n)
X=cbind(I, X1)

# GÈnÈration des donnÈes
beta0 = c(2, .5)
# ModËle Logit
y = rbinom(n, size = 1, prob = plogis(X %*% beta0))
 
n=2000, 20000
I=rep(1, n)
X1=runif(n)
X=cbind(I, X1)

beta0 = c(.4, -.3)
# ModËle Probit
y = rbinom(n, size = 1, prob = pnorm(X %*% beta0))
# plots :
par(mfrow=c(2, 3))
i=
plot(estimats.e.1[,i], type="l", ylab="", main="a0=.1")
plot(estimats.e.2[,i], type="l", ylab="", main="a0=.2")
plot(estimats.e.3[,i], type="l", ylab="", main="a0=.3")
plot(estimats.e.4[,i], type="l", ylab="", main="a0=.4")
plot(estimats.e.5[,i], type="l", ylab="", main="a0=.5")

# Hist all2.hist
par(mfrow=c(2, 5))
hist(estimats.e.1[,9], xlab="", main="a0=.1", col="gray")
hist(estimats.e.2[,9], xlab="", main="a0=.2", col="gray")
hist(estimats.e.3[,9], xlab="", main="a0=.3", col="gray")
hist(estimats.e.4[,9], xlab="", main="a0=.4", col="gray")
hist(estimats.e.5[,9], xlab="", main="a0=.5", col="gray")
hist(estimats.e.1[,9], xlab="", main="", col="black")
hist(estimats.e.2[,9], xlab="", main="", col="black")
hist(estimats.e.3[,9], xlab="", main="", col="black")
hist(estimats.e.4[,9], xlab="", main="", col="black")
hist(estimats.e.5[,9], xlab="", main="", col="black")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Both models have a common parameter theta :
 
flogis=function(theta1, y=y, X=X){
    prob=plogis(X%*%matrix(theta1,ncol=1))
    plgs=(prob^y) * (1-prob)^(1-y)
return(plgs)
}

fprobit=function(theta2, y=y, X=X){
  prob=pnorm(X%*%matrix(theta2,ncol=1))
  prbit=(prob^y) * (1-prob)^(1-y)
return(prbit)
}

# log function :
logfpost=function(theta, X, y1, y2, x1, x2, c){
problgs=plogis(x1%*%matrix(theta,ncol=1))
    plgs=(problgs^y1) * (1-problgs)^(1-y1)

probpbit=pnorm(x2%*%matrix(theta,ncol=1))
  prbit=(probpbit^y2) * (1-probpbit)^(1-y2)

pmelang=c(plgs, prbit)
logpmelang=log(pmelang)
slogpmelang=sum(logpmelang)

q=length(theta)

mean=rep(0,q)
txx=t(X)%*%X
stxx=solve(txx)
sigma=c*stxx 
dnor= dmvnorm(theta, mean, sigma, log=TRUE)
return(dnor+ slogpmelang)
}

# Gibbs parameters :Correct
GSmix_theta5=function(Nsim, e, tau, y, X){
k=length(e)
n=length(y)
if(is.vector(X)!= TRUE){
q=dim(X)[2]}else{
q=length(X)}

theta=matrix(0, nrow=Nsim, ncol=q)
p=matrix(0, nrow=Nsim, ncol=k)
p[1,]= p0
mle1 = summary(glm(y~X, family=binomial))
  beta1= mle1$coefficients[,1]
  sigma.asymp1= mle1$cov.unscaled
theta[1,]= beta1

for(t in 2:Nsim){

z=vector("numeric", length=n)
for(l in 1:n){
z[l]=1+(runif(1)<1/(1+p[t-1,1]* flogis(theta[t-1,], y[l], X[l,])/((1-p[t-1,1])*fprobit(theta[t-1,], y[l], X[l,]))))
}

  #pz=p[t-1,]
  #       z=vector("numeric", length=n)
     #    probs=matrix(0, nrow=k, ncol=n)      
#probs[1,]=pz[1]*flogis(theta1[t-1,], y=y, X=X)
#probs[2,]=pz[2]*fprobit(theta2[t-1,], y=y, X=X)
#for(i in 1:n){
#                   z[i]=which.max(probs[,i])
#}

L=vector("numeric", length=k)
        for(I in 1:k){
L[I]=length(z[z==I])
}

alphap=L+e
p[t, ]=rdirichlet(1, alphap)

           xbar1=X[z==1,]
xbar1=matrix(xbar1, ncol=q)
           xx1=t(xbar1)%*%xbar1
           xbar2=X[z==2,]
xbar2=matrix(xbar2, ncol=q)
           xx2=t(xbar2)%*%xbar2
           yz1=y[z==1]
           yz2=y[z==2]

proposal1=rmvnorm(1, theta[t-1,], tau*sigma.asymp1)
#proposal1=rmvnorm(1, beta1, tau*sigma.asymp1)

    logalpha1=min(0, logfpost(proposal1, X, yz1, yz2, xbar1, xbar2, c) - logfpost(theta[t-1,], X, yz1, yz2, xbar1, xbar2, c))
    if(log(runif(1))<logalpha1){ theta[t, ]= proposal1} else{theta[t,]=theta[t-1,]}
 
}
            TP=cbind(theta, p)

return(TP)   
}

# Convergent MH :

GSmix_theta5=function(Nsim, e, tau, y, X){
k=length(e)
n=length(y)
if(is.vector(X)!= TRUE){
q=dim(X)[2]}else{
q=length(X)}
c=n
a0=e[1]

theta=matrix(0, nrow=Nsim, ncol=q)
p=matrix(0, nrow=Nsim, ncol=k)
p[1,]= p0
mle1 = summary(glm(y~X, family=binomial))
  beta1= mle1$coefficients[,1]
  sigma.asymp1= mle1$cov.unscaled

  theta[1,]= beta1

tx=t(X)%*%X
v1=solve(tx)
v2=c*v1
mean=rep(0, q)

for(t in 2 :Nsim){
#Hasting independent proposal theta:
cur=(p[t-1,1]*flogis(theta[t-1,], y=y, X=X))+(p[t-1,2]*fprobit(theta[t-1,], y=y, X=X))
lcur=log(cur)

curlik=sum(lcur)+ log(dmvnorm(theta[t-1,], mean, v2))
proposal1=rmvnorm(1, theta[t-1,], tau*sigma.asymp1)
#proposal1=rmvnorm(1, beta1, tau*sigma.asymp1)

pr=(p[t-1,1]*flogis(proposal1, y=y, X=X))+(p[t-1,2]*fprobit(proposal1, y=y, X=X))
lpr=log(pr)
proplik= sum(lpr) + log(dmvnorm(c(proposal1), mean, v2))
theta[t,]=theta[t-1,]

alp= proplik-curlik

alpha1=min(0, alp)

if (log(runif(1))<alpha1){
     theta[t,]=proposal1}

z=1+(runif(n)<1/(1+p[t-1,1]*flogis(theta[t-1,], y=y, X=X)/((1-p[t-1,1])*fprobit(theta[t-1,], y=y, X=X))))
L=vector("numeric", length=k)
        for(I in 1:k){
L[I]=length(z[z==I])
}
alphap=L+e
p[t, ]=rdirichlet(1, alphap)

  #p[t,1]=1/(1+rgamma(1,a0+sum((z==2)))/rgamma(1,a0+sum((z==1))))
  #p[t,2]=1-p[t,1]

cur1= (p[t,1]*flogis(theta[t,], y=y, X=X))+(p[t,2]*fprobit(theta[t,], y=y, X=X))
lcur1=log(cur1)
curlik1= sum(lcur1) 

propos=rbeta(1,a0,a0)
pr1= (propos*flogis(theta[t,], y=y, X=X))+((1-propos)*fprobit(theta[t,], y=y, X=X))
lpr1=log(pr1)
proplik1= sum(lpr1) 
  if (runif(1)<exp(proplik1-curlik1)){
    p[t,1]=propos
p[t,2]=1- p[t,1]
 }
}
            TP=cbind(theta, p)

return(TP)   
}
c= dim(X)[1]
# Estimations of alpha for different values of tau :

tau=1
e=c(.1, .1)
estimats.e.1=GSmix_theta5(Nsim, e, tau, y, X)
e=c(.2, .2)
estimats.e.2=GSmix_theta5(Nsim, e, tau, y, X)
e=c(.3, .3)
estimats.e.3=GSmix_theta5(Nsim, e, tau, y, X)
e=c(.4, .4)
estimats.e.4=GSmix_theta5(Nsim, e, tau, y, X)
e=c(.5, .5)
estimats.e.5=GSmix_theta5(Nsim, e, tau, y, X)

# Another example with juste one explanatory variable like Pima data :
n=1000, 10000
I=rep(1, n)
X1=runif(n)
X=cbind(I, X1)

# GÈnÈration des donnÈes
beta0 = c(2, .5)
# ModËle Logit
y = rbinom(n, size = 1, prob = plogis(X %*% beta0))
 
n=2000, 20000
I=rep(1, n)
X1=runif(n)
X=cbind(I, X1)

beta0 = c(.4, -.3)
# ModËle Probit
y = rbinom(n, size = 1, prob = pnorm(X %*% beta0))
# plots :
par(mfrow=c(2, 3))
i=
plot(estimats.e.1[,i], type="l", ylab="", main="a0=.1")
plot(estimats.e.2[,i], type="l", ylab="", main="a0=.2")
plot(estimats.e.3[,i], type="l", ylab="", main="a0=.3")
plot(estimats.e.4[,i], type="l", ylab="", main="a0=.4")
plot(estimats.e.5[,i], type="l", ylab="", main="a0=.5")

# Hist all2.hist
par(mfrow=c(2, 5))
hist(estimats.e.1[,9], xlab="", main="a0=.1", col="gray")
hist(estimats.e.2[,9], xlab="", main="a0=.2", col="gray")
hist(estimats.e.3[,9], xlab="", main="a0=.3", col="gray")
hist(estimats.e.4[,9], xlab="", main="a0=.4", col="gray")
hist(estimats.e.5[,9], xlab="", main="a0=.5", col="gray")
hist(estimats.e.1[,9], xlab="", main="", col="black")
hist(estimats.e.2[,9], xlab="", main="", col="black")
hist(estimats.e.3[,9], xlab="", main="", col="black")
hist(estimats.e.4[,9], xlab="", main="", col="black")
hist(estimats.e.5[,9], xlab="", main="", col="black")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Marginal likelihood Probit:
mrg_probit=function(theta, thetahat, cv, y, X){
v=2*cv
q=dim(theta)[1]
h=dim(theta)[2]
n=length(y)
sm=rep(1, q)
tx=t(X)%*%X
v1=solve(tx)
v2=c*v1
mean=rep(0, h)

for(m in 1:q){
dnrm=dmvnorm(theta[m,], thetahat, v)
dnrm1= dmvnorm(theta[m,], mean, v2)
probpbit=pnorm(X%*%matrix(theta[m,],ncol=1))
  prbit=(probpbit^y) * (1-probpbit)^(1-y)
pprbit=prod(prbit)
sm[m]=pprbit*dnrm1*dnrm
}
marginal=sum(sm)/q
return(marginal)
}

thetahat=c(mean(estimats.e.1[,1]), mean(estimats.e.1[,2]))
theta=cbind(estimats.e.1[,1], estimats.e.1[,2])
cv=cov(theta)

n=length(y)
mrg_p=rep(0, n)
for(i in 1:n){
mrg_p[i]= mrg_probit(theta, thetahat, cv, y[i], X[i,])
}

# Marginal likelihood Logit:
mrg_logit=function(theta, thetahat, cv, y, X){
v=2*cv
q=dim(theta)[1]
n=length(y)
sm=rep(1, q)
tx=t(X)%*%X
v1=solve(tx)
v2=c*v1
mean=rep(0, h)

for(m in 1:q){
dnrm=dmvnorm(theta[m,], thetahat, v)
dnrm1= dmvnorm(theta[m,], mean, v2)
probpbit=plogis(X%*%matrix(theta[m,],ncol=1))
  prbit=(probpbit^y) * (1-probpbit)^(1-y)
pprbit=prod(prbit)
sm[m]=pprbit*dnrm1*dnrm
}
marginal=sum(sm)/q
return(marginal)
}

mrg_p=rep(0, n)
for(i in 1:n){
mrg_p[i]= mrg_logit(cc, thetahat, cv, y[i], X[i,])
}


Simple approximation to normal distribution
29 April 2010 John Math, Statistics
Here’s a simple approximation to the normal distribution I just ran across. The density function is
f(x) = (1 + cos(x))/2π
over the interval (-π, π). The plot below graphs this density with a solid blue line. For comparison, the density of a normal distribution with the same variance is plotted with a dashed orange line.

The approximation is good enough to use for teaching. Students may benefit from doing an exercise twice, once with this approximation and then again with the normal distribution. Having an approximation they can integrate in closed form may help take some of the mystery out of the normal distribution.
The approximation may have practical uses. The agreement between the PDFs isn’t great. However, the agreement between the CDFs (which is more important) is surprisingly good. The maximum difference between the two CDFs is only 0.018. (The differences between the PDFs oscillate, and so their integrals, the CDFs, are closer together.)
I ran across this approximation here. It goes back to the 1961 paper “A cosine approximation to the normal distribution” by D. H. Raab and E. H. Green, Psychometrika, Volume 26, pages 447-450.
Update 1: See the paper referenced in the first comment. It gives a much more accurate approximation using a logistic function. The cosine approximation is a little simpler and may be better for teaching. However, the logistic approximation has infinite support. That could be an advantage since students might be distracted by the finite support of the cosine approximation.
The logistic approximation for the standard normal CDF is
F(x) = 1/(1 + exp(-0.07056 x^3 – 1.5976 x))
and has a maximum error of 0.00014 at x = ± 3.16.

# using the approximation above we have :


So the algorithms change to the following :
# Both models have parameter theta with the same sens in both model:
 
flogis=function(theta1, y=y, X=X){
    prob=plogis(X%*%matrix(theta1,ncol=1))
    plgs=(prob^y) * (1-prob)^(1-y)
return(plgs)
}

fprobit=function(theta2, y=y, X=X){
  prob=pnorm(X%*%matrix(theta2,ncol=1))
  prbit=(prob^y) * (1-prob)^(1-y)
return(prbit)
}

# log function :
logfpost=function(theta, X, y, y1, y2, x1, x2, c){
problgs=plogis(x1%*%matrix(theta,ncol=1))
    plgs=(problgs^y1) * (1-problgs)^(1-y1)

k=sqrt(2/pi)
#rat= summary(glm(y~X, family=binomial(logit)))$coefficients[,1]/summary(glm(y~X, family=binomial(probit)))$coefficients[,1]
#thetap=theta/rat
thetap=theta/(2*k)
probpbit=pnorm(x2%*%matrix(thetap,ncol=1))
  prbit=(probpbit^y2) * (1-probpbit)^(1-y2)

pmelang=c(plgs, prbit)
logpmelang=log(pmelang)
slogpmelang=sum(logpmelang)

q=length(theta)

mean=rep(0,q)
txx=t(X)%*%X
stxx=solve(txx)
sigma=c*stxx 
dnor= dmvnorm(theta, mean, sigma, log=TRUE)
return(dnor+ slogpmelang)
}

# Gibbs parameters :Correct
GSmix_theta5=function(Nsim, e, tau, y, X){
k=length(e)
n=length(y)
if(is.vector(X)!= TRUE){
q=dim(X)[2]}else{
q=length(X)}
c=n

theta=matrix(0, nrow=Nsim, ncol=q)
p=matrix(0, nrow=Nsim, ncol=k)
p[1,]= p0
mle1 = summary(glm(y~X, family=binomial))
  beta1= mle1$coefficients[,1]
  sigma.asymp1= mle1$cov.unscaled
theta[1,]= beta1
kk=sqrt(2/pi)
thetap=theta/(2*kk)
#rat= summary(glm(y~X, family=binomial(logit)))$coefficients[,1]/summary(glm(y~X, family=binomial(probit)))$coefficients[,1]
#thetap=theta/rat

for(t in 2:Nsim){

z=vector("numeric", length=n)

for(l in 1:n){
z[l]=1+(runif(1)<1/(1+p[t-1,1]* flogis(theta[t-1,], y[l], X[l,])/((1-p[t-1,1])*fprobit(thetap[t-1,], y[l], X[l,]))))
}

  #pz=p[t-1,]
  #       z=vector("numeric", length=n)
     #    probs=matrix(0, nrow=k, ncol=n)      
#probs[1,]=pz[1]*flogis(theta1[t-1,], y=y, X=X)
#probs[2,]=pz[2]*fprobit(theta2[t-1,], y=y, X=X)
#for(i in 1:n){
#                   z[i]=which.max(probs[,i])
#}

L=vector("numeric", length=k)
        for(I in 1:k){
L[I]=length(z[z==I])
}

alphap=L+e
p[t, ]=rdirichlet(1, alphap)

           xbar1=X[z==1,]
xbar1=matrix(xbar1, ncol=q)
           xx1=t(xbar1)%*%xbar1
           xbar2=X[z==2,]
xbar2=matrix(xbar2, ncol=q)
           xx2=t(xbar2)%*%xbar2
           yz1=y[z==1]
           yz2=y[z==2]

proposal1=rmvnorm(1, theta[t-1,], tau*sigma.asymp1)
proposal1p= proposal1/(2*kk)
# proposal1p= proposal1/rat

#proposal1=rmvnorm(1, beta1, tau*sigma.asymp1)

    logalpha1=min(0, logfpost(proposal1, X, y, yz1, yz2, xbar1, xbar2, c) - logfpost(theta[t-1,], X, y, yz1, yz2, xbar1, xbar2, c))
    if(log(runif(1))<logalpha1){ theta[t, ]= proposal1 ; thetap[t, ]= proposal1p } else{theta[t,]=theta[t-1,] ; thetap[t, ]= thetap[t-1,]}
 
}
            TP=cbind(theta, thetap, p)

return(TP)   
}

# Convergent MH :

GSmix_theta5=function(Nsim, e, tau, y, X){
k=length(e)
n=length(y)
if(is.vector(X)!= TRUE){
q=dim(X)[2]}else{
q=length(X)}
c=n
a0=e[1]

theta=matrix(0, nrow=Nsim, ncol=q)
p=matrix(0, nrow=Nsim, ncol=k)
p[1,]= p0
mle1 = summary(glm(y~X, family=binomial))
  beta1= mle1$coefficients[,1]
  sigma.asymp1= mle1$cov.unscaled

  theta[1,]= beta1
kk=sqrt(2/pi)
thetap=theta/(2*kk)


tx=t(X)%*%X
v1=solve(tx)
v2=c*v1
mean=rep(0, q)

for(t in 2 :Nsim){
#Hasting independent proposal theta:
cur=(p[t-1,1]*flogis(theta[t-1,], y=y, X=X))+(p[t-1,2]*fprobit(thetap[t-1,], y=y, X=X))
lcur=log(cur)

curlik=sum(lcur)+ log(dmvnorm(theta[t-1,], mean, v2))
proposal1=rmvnorm(1, theta[t-1,], tau*sigma.asymp1)
proposal1p=proposal1/(2*kk)

#proposal1=rmvnorm(1, beta1, tau*sigma.asymp1)

pr=(p[t-1,1]*flogis(proposal1, y=y, X=X))+(p[t-1,2]*fprobit(proposal1p, y=y, X=X))
lpr=log(pr)
proplik= sum(lpr) + log(dmvnorm(c(proposal1), mean, v2))
theta[t,]=theta[t-1,]
thetap[t,]=thetap[t-1,]

alp= proplik-curlik

alpha1=min(0, alp)

if (log(runif(1))<alpha1){
     theta[t,]=proposal1
thetap[t,]=proposal1p}

z=1+(runif(n)<1/(1+p[t-1,1]*flogis(theta[t-1,], y=y, X=X)/((1-p[t-1,1])*fprobit(thetap[t-1,], y=y, X=X))))
L=vector("numeric", length=k)
        for(I in 1:k){
L[I]=length(z[z==I])
}
alphap=L+e
p[t, ]=rdirichlet(1, alphap)

  #p[t,1]=1/(1+rgamma(1,a0+sum((z==2)))/rgamma(1,a0+sum((z==1))))
  #p[t,2]=1-p[t,1]

cur1= (p[t,1]*flogis(theta[t,], y=y, X=X))+(p[t,2]*fprobit(thetap[t,], y=y, X=X))
lcur1=log(cur1)
curlik1= sum(lcur1) 

propos=rbeta(1,a0,a0)
pr1= (propos*flogis(theta[t,], y=y, X=X))+((1-propos)*fprobit(thetap[t,], y=y, X=X))
lpr1=log(pr1)
proplik1= sum(lpr1) 
  if (runif(1)<exp(proplik1-curlik1)){
    p[t,1]=propos
p[t,2]=1- p[t,1]
 }
}
            TP=cbind(theta, thetap, p)

return(TP)   
}
# Estimations of alpha for different values of tau :
c= dim(X)[1]
tau=1
e=c(.1, .1)
estimats.e.1=GSmix_theta5(Nsim, e, tau, y, X)
e=c(.2, .2)
estimats.e.2=GSmix_theta5(Nsim, e, tau, y, X)
e=c(.3, .3)
estimats.e.3=GSmix_theta5(Nsim, e, tau, y, X)
e=c(.4, .4)
estimats.e.4=GSmix_theta5(Nsim, e, tau, y, X)
e=c(.5, .5)
estimats.e.5=GSmix_theta5(Nsim, e, tau, y, X)

# Another example with juste one explanatory variable like Pima data :
n=100000
I=rep(1, n)
X1=runif(n)
X=cbind(I, X1)

# GÈnÈration des donnÈes
beta0 = c(.6, 1.5)
# ModËle Logit
y = rbinom(n, size = 1, prob = plogis(X %*% beta0))
 

beta0 = c(-1.2, 4)
# ModËle Probit
y = rbinom(n, size = 1, prob = pnorm(X %*% beta0))

summary(glm(y~X, family=binomial(logit)))$coefficients[,1]/summary(glm(y~X, family=binomial(probit)))$coefficients[,1]
(Intercept)         XX1 
   1.730236    1.741938
# In the case where the ratio of the max likelihood of logit over probit is 1.73, our algorithm is correct and alpha support the true model from which our data is simulated. 
# plots :
par(mfrow=c(2, 3))
i=
plot(estimats.e.1[,i], type="l", ylab="", main="a0=.1")
plot(estimats.e.2[,i], type="l", ylab="", main="a0=.2")
plot(estimats.e.3[,i], type="l", ylab="", main="a0=.3")
plot(estimats.e.4[,i], type="l", ylab="", main="a0=.4")
plot(estimats.e.5[,i], type="l", ylab="", main="a0=.5")
df=data.frame(estimats.e.1[,i], estimats.e.2[,i], estimats.e.3[,i], estimats.e.4[,i], estimats.e.5[,i])
 names(df)=xl
 boxplot(df, col="gray")
# Hist all2.hist
par(mfrow=c(2, 5))
hist(estimats.e.1[,9], xlab="", main="a0=.1", col="gray")
hist(estimats.e.2[,9], xlab="", main="a0=.2", col="gray")
hist(estimats.e.3[,9], xlab="", main="a0=.3", col="gray")
hist(estimats.e.4[,9], xlab="", main="a0=.4", col="gray")
hist(estimats.e.5[,9], xlab="", main="a0=.5", col="gray")
hist(estimats.e.1[,9], xlab="", main="", col="black")
hist(estimats.e.2[,9], xlab="", main="", col="black")
hist(estimats.e.3[,9], xlab="", main="", col="black")
hist(estimats.e.4[,9], xlab="", main="", col="black")
hist(estimats.e.5[,9], xlab="", main="", col="black")

# Some other experience : MH that always gives us the true model.
# Convergent MH :

GSmix_theta5=function(Nsim, e, tau, y, X){
k=length(e)
n=length(y)
if(is.vector(X)!= TRUE){
q=dim(X)[2]}else{
q=length(X)}
c=n
a0=e[1]

theta=matrix(0, nrow=Nsim, ncol=q)
p=matrix(0, nrow=Nsim, ncol=k)
p[1,]= p0
mle1 = summary(glm(y~X, family=binomial))
  beta1= mle1$coefficients[,1]
  sigma.asymp1= mle1$cov.unscaled
rat=summary(glm(y~X, family=binomial(logit)))$coefficients[,1]/summary(glm(y~X, family=binomial(probit)))$coefficients[,1]
  theta[1,]= beta1
 
thetap=theta/rat


tx=t(X)%*%X
v1=solve(tx)
v2=c*v1
mean=rep(0, q)

for(t in 2 :Nsim){
#Hasting independent proposal theta:
cur=(p[t-1,1]*flogis(theta[t-1,], y=y, X=X))+(p[t-1,2]*fprobit(thetap[t-1,], y=y, X=X))
lcur=log(cur)

curlik=sum(lcur)+ log(dmvnorm(theta[t-1,], mean, v2))
proposal1=rmvnorm(1, theta[t-1,], tau*sigma.asymp1)
proposal1p=proposal1/rat

#proposal1=rmvnorm(1, beta1, tau*sigma.asymp1)

pr=(p[t-1,1]*flogis(proposal1, y=y, X=X))+(p[t-1,2]*fprobit(proposal1p, y=y, X=X))
lpr=log(pr)
proplik= sum(lpr) + log(dmvnorm(c(proposal1), mean, v2))
theta[t,]=theta[t-1,]
thetap[t,]=thetap[t-1,]

alp= proplik-curlik

alpha1=min(0, alp)

if (log(runif(1))<alpha1){
     theta[t,]=proposal1
thetap[t,]=proposal1p}

z=1+(runif(n)<1/(1+p[t-1,1]*flogis(theta[t-1,], y=y, X=X)/((1-p[t-1,1])*fprobit(thetap[t-1,], y=y, X=X))))
L=vector("numeric", length=k)
        for(I in 1:k){
L[I]=length(z[z==I])
}
alphap=L+e
p[t, ]=rdirichlet(1, alphap)

  #p[t,1]=1/(1+rgamma(1,a0+sum((z==2)))/rgamma(1,a0+sum((z==1))))
  #p[t,2]=1-p[t,1]

cur1= (p[t,1]*flogis(theta[t,], y=y, X=X))+(p[t,2]*fprobit(thetap[t,], y=y, X=X))
lcur1=log(cur1)
curlik1= sum(lcur1) 

propos=rbeta(1,a0,a0)
pr1= (propos*flogis(theta[t,], y=y, X=X))+((1-propos)*fprobit(thetap[t,], y=y, X=X))
lpr1=log(pr1)
proplik1= sum(lpr1) 
  if (runif(1)<exp(proplik1-curlik1)){
    p[t,1]=propos
p[t,2]=1- p[t,1]
 }
}
            TP=cbind(theta, thetap, p)

return(TP)   
}

1)
n=10000
 I=rep(1, n)
 X1=rnorm(n)
 X=cbind(I, X1)
beta0=c(5, 1.5)
 y = rbinom(n, size = 1, prob = plogis(X %*% beta0))
 summary(glm(y~X, family=binomial(logit)))$coefficients[,1]/summary(glm(y~X, family=binomial(probit)))$coefficients[,1]
(Intercept)         XX1 
   2.006829    2.296949

beta0=c(2.3, .5)
  y = rbinom(n, size = 1, prob = pnorm(X %*% beta0))
  summary(glm(y~X, family=binomial(logit)))$coefficients[,1]/summary(glm(y~X, family=binomial(probit)))$coefficients[,1]
(Intercept)         XX1 
   2.147218    2.260361

2)
n=10000
I=rep(1, n)
X1=runif(n)
X=cbind(I, X1)

# GÈnÈration des donnÈes
beta0 = c(2, .5)
# ModËle Logit
y = rbinom(n, size = 1, prob = plogis(X %*% beta0))

 summary(glm(y~X, family=binomial(logit)))$coefficients[,1]/summary(glm(y~X, family=binomial(probit)))$coefficients[,1]
(Intercept)         XX1 
   1.699867    1.969545

beta0 = c(.4, -.3)
# ModËle Probit
y = rbinom(n, size = 1, prob = pnorm(X %*% beta0))

summary(glm(y~X, family=binomial(logit)))$coefficients[,1]/summary(glm(y~X, family=binomial(probit)))$coefficients[,1]
(Intercept)         XX1 
   1.606054    1.612486

3)
n=10000
 I=rep(1, n)
 X1=rnorm(n)
X2=runif(n)
X3=runif(n, -1, 1)
 X=cbind(I, X1, X2, X3)
beta0=c(5, 1.5, -.5, 2)
 y = rbinom(n, size = 1, prob = plogis(X %*% beta0))
 summary(glm(y~X, family=binomial(logit)))$coefficients[,1]/summary(glm(y~X, family=binomial(probit)))$coefficients[,1]
(Intercept)         XX1         XX2         XX3 
   1.961665    2.097110    2.183348    2.172798

beta0=c(-2.3, .5, 1, -.7)
  y = rbinom(n, size = 1, prob = pnorm(X %*% beta0))
  summary(glm(y~X, family=binomial(logit)))$coefficients[,1]/summary(glm(y~X, family=binomial(probit)))$coefficients[,1]
(Intercept)         XX1         XX2         XX3 
   1.840712    1.920629    1.919713    1.983582

4)
n=20000
 I=rep(1, n)
 X1=rnorm(n)
X2=runif(n)
X3=runif(n, -1, 1)
X4=rnorm(n, 1, 4)
X5=runif(n, 2, 3)
X6=runif(n, -1, 1)

 X=cbind(I, X1, X2, X3, X4, X5, X6)
beta0=c(3, 1.5, -.5, 2, -.3, 1.1, -.8)
 y = rbinom(n, size = 1, prob = plogis(X %*% beta0))
 summary(glm(y~X, family=binomial(logit)))$coefficients[,1]/summary(glm(y~X, family=binomial(probit)))$coefficients[,1]
(Intercept)         XX1         XX2         XX3         XX4         XX5         XX6 
   1.865496    2.035094    2.385370    1.999091    2.013937    2.032902    2.040194

beta0=c(-2.3, .5, 1, -.7, .4, 1.4, -3)
  y = rbinom(n, size = 1, prob = pnorm(X %*% beta0))
  summary(glm(y~X, family=binomial(logit)))$coefficients[,1]/summary(glm(y~X, family=binomial(probit)))$coefficients[,1]
(Intercept)         XX1         XX2         XX3         XX4         XX5         XX6 
   1.785770    1.804756    1.793429    1.796141    1.807808    1.796413    1.808213


5)
n=100000
 I=rep(1, n)
 X1=rnorm(n)
X2=runif(n)
X3=runif(n, -1, 1)
 X=cbind(I, X1, X2, X3)
beta0=c(5, 1.5, -.5, 2)
 y = rbinom(n, size = 1, prob = plogis(X %*% beta0))
 summary(glm(y~X, family=binomial(logit)))$coefficients[,1]/summary(glm(y~X, family=binomial(probit)))$coefficients[,1]
(Intercept)         XX1         XX2         XX3 
   1.959995    2.088305    2.153621    2.153811

beta0=c(-2.3, .5, 1, -.7)
  y = rbinom(n, size = 1, prob = pnorm(X %*% beta0))
  summary(glm(y~X, family=binomial(logit)))$coefficients[,1]/summary(glm(y~X, family=binomial(probit)))$coefficients[,1]
(Intercept)         XX1         XX2         XX3 
   1.863254    1.951561    1.965914    1.986160


6)
Prima dataset
