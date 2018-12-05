###program mixture for GoF tests with DPM 
## model : alpha N(mu,sigma^2) + (1-alpha)int N(mu',sigma^2)dP(mu')
## prior : alpha : beta(a1,a1)
## mu1 : N(0,sigma^2*tau), sigma^2 = IG(b1,b2)
## P = DP(A, N(0, sigma^2*tau))

### this program is a marginal version : where we also integrate out alpha to avoid the local 
#modes near alpha = 0 or 1

posteriorDP3 = function(T, data1, a1, a2, tau, b1, b2, A){
	
	 Y = data1 ## data
	 n = length(Y)
	 print(n)
	 
	 
	 #true parameters
	# Z = data1$Z
	 #sigma = data1$sigma
	 #mu1 = data1$mu1
	 
	 ## we initialize considering everyone in gourp 1 (normal) except 1
	 K = 1
	 
	 alpha = 0.5
	 n1 = n/2 ## be careful if n not even
	 n2 = n/2
	 mu1 = 0
	
	 
	 ## results stocked into
	 BZ = matrix(0, nrow = T,ncol = n)
	 vecalpha = rep(0,T)
	 vecmu = rep(0,T)
	 vecsigma1 = rep(0,T)
	 vecsigma2 = rep(0,T)
	 vecK = rep(0,T)
	 Bn1 = rep(0,T)
	 bt = 1
	 
	 
	 U = 0
	
	 vecalpha[1]= alpha
	 vecmu[1]= mu1
	 vecsigma1[1] = 1
	 vecsigma2[1] = 1
	 
	 Z = c(rep(0,n/2),rep(1,n/2))
	 Z1 = Z
	 Y1 = Y[Z==0]
	 n1 = sum(Z==0)
	 mY1 = sum(Y1)/n1
	 Y2 = Y[Z==1]
	 n2 = sum(Z==1)
	 mY2 = sum(Y2)/n2
	 nmY1 = mY1
	 nmY1i = mY1
	 nmY2= mY2
	 nmY2i= nmY2
	 SY1 = sum((Y1- mY1)^2)
	 SY2 = sum((Y2- mY2)^2)
	 nSY1 = SY1
	 nSY2= SY2
	  nSY1i = nSY1
	
	 nSY2i=nSY2
	 Yi = rep(0, n-1)
	 Zi = rep(0,n-1)
	 BZ[1, ] = Z
	 vecp = c(1,0)
	 vecn = c(n/2,n/2)
	 Vn2=0
	  K = 1 
	  nzi = 0
	  Knew= K
	  nnew= vecn
	  Z1=Z
    p = 0
	 
	 
	 
	## iterations
	for (t in (2:T)){
	 # print(t)

	  
	  n1 = vecn[1]
	  ## we now reactualise the means and "variance " in each group since the groups have changed
	 
	  ## update sigma1 = sigma2 now
	  bt = b2 + SY1/2+n1*mY1^2/(2*(tau^2*n1+1)) 
	
	 # vecsigma1[t] = 1/sqrt(rgamma(1,b1+n1/2,bt)) 
	
	    ## update sigma2
	  if (K>0){
	 bt = bt +sum(SY2)/2+ sum(vecn[2:(K+1)]*mY2^2/(2*(tau^2*vecn[2:(K+1)]+1)))  
	 }
	  vecsigma1[t] = 1/sqrt(rgamma(1,b1+ n/2,bt)) 
	  vecsigma2[t] = vecsigma1[t]
	  
		vecmu[t] = rnorm(1,n1*mY1/(n1+ 1/tau^2), vecsigma1[t]/sqrt(n1+1/tau^2))
	
		vecalpha[t] = rbeta(1,a1+n1,a2+n-n1)
		
		#allocations of both Z=0 / Z>0 and Z=j , j>0
		# computation of proba for Z and update of Z
		for (i in (1:n))
		  {
		  ### we define a new vector of nj' without i : called nnew
		  nzi = Z[i] ## old allocation for i
		  nnew = vecn #[] - ((0:K) == nzi)
		  
		  nnew[nzi+1] = nnew[nzi+1]-1
		  n1 = nnew[1]
		  
		  ## construction of Z2 : new vector of Z if one less group of size n-1;  same with Y
		  Zi = Z[-i]
		  Yi  = Y[-i]
		### empirical   means and variances for group 0
		  if (Z[i]==0){
		    if (n1>0){
		    nmY1 = (n1+1)*mY1/n1 - Y[i]/n1
	      nSY1 = sum((Yi[Zi==0]-nmY1)^2)
		    }
		    else{nmY1 = 0 
		         nSY1 = 0
		        }
	
		    nmY1i = mY1
		    nSY1i = SY1
		    nmY2 = mY2
		    nSY2 = SY2
		  }
		  
		  if (Z[i] >0){
		    nmY1 = mY1
		    nSY1 = SY1
		    nmY1i = mY1*n1/(n1+1) + Y[i]/(n1+1)
		    nSY1i = sum((Yi[Zi==0]-nmY1i)^2)+ (Y[i]-nmY1i)^2
		  }
		  
		 p= exp(lgamma(n1+1+a1)+lgamma(n-n1-1+a2)-lgamma(n+a1+a2) -(nSY1i-nSY1 + (n1+1)*nmY1i^2/((n1+1)*tau^2+1) - n1*nmY1^2/(n1*tau^2+1))/(2*vecsigma1[t]^2))*sqrt(n1*tau^2+1)/(vecsigma1[t]*sqrt((n1+1)*tau^2+1)
		 )		  
		  
		  ### determination of new means and variances (empirical)
		 
		  
		  
		   if (K>0){
		    Vn2 = nnew[2:(K+1)]
		    Knew = sum((Vn2>0))
		    
		  
		    if (Knew >0){
		      Z1=Z
		      nnew = c(n1,Vn2[(Vn2>0)])
		      
		      if (Knew < K){
		        Z1[Z>nzi] = Z1[Z>nzi]-1
		        Z = Z1
	#	        Zi= Z[-i]
		    
		        nmY2 = mY2[(Vn2>0)]
		#        print(nmY2)
		        nSY2 = SY2[(Vn2>0)]
		        nmY2i = nmY2*nnew[2:(Knew+1)]/(nnew[2:(Knew+1)]+1) + Y[i]/(nnew[2:(Knew+1)]+1)
		        nSY2i = nSY2 + nnew[2:(Knew+1)]*nmY2^2+ Y[i]^2 - (nnew[2:(Knew+1)]+1)*nmY2i^2
		      }
		      if (Knew==K){
		          if (Z[i]>0){
		            nmY2 = mY2*vecn[2:(Knew+1)]/nnew[2:(Knew+1)] - ((1:K)==Z[i])*Y[i]/nnew[2:(Knew+1)]
		            nSY2 = SY2 + mY2^2*vecn[2:(Knew+1)] - ((1:K)==Z[i])*Y[i]^2 - nnew[2:(Knew+1)]*nmY2^2
		           # print("vecn")
		          #  print(vecn)
		           # print(nnew)
		          #  print(Knew)
		          }
		        nmY2i = nmY2*nnew[2:(Knew+1)]/(nnew[2:(Knew+1)]+1) + Y[i]/(nnew[2:(Knew+1)]+1)
		        nSY2i = nSY2 + nnew[2:(Knew+1)]*nmY2^2+ Y[i]^2 - (nnew[2:(Knew+1)]+1)*nmY2i^2
		     #   if (nSY2i< 0){print("nSY2i if negative")
		    #    print(nSY2i)}
		      }
		      
		      
		  
		      vecp = c(p,exp(lgamma(n1+a1)+lgamma(n-n1+a2)-lgamma(n+a1+a2)-(nSY2i-nSY2 + (nnew[2:(Knew+1)]+1)*nmY2i^2/((nnew[2:(Knew+1)]+1)*tau^2+1) - nnew[2:(Knew+1)]*nmY2^2/(nnew[2:(Knew+1)]*tau^2+1))/(2*vecsigma2[t]^2))*sqrt(nnew[2:(Knew+1)]*tau^2+1)/(vecsigma2[t]*sqrt((nnew[2:(Knew+1)]+1)*tau^2+1))*nnew[2:(Knew+1)]/(n-n1+A-1))
		      vecp = c(vecp, exp(lgamma(n1+a1)+lgamma(n-n1+a2)-lgamma(n+a1+a2))*dnorm(Y[i],0,sd= vecsigma2[t]*sqrt(tau^2+1))*A/(A+n-n1-1))
		    }
		    if (Knew==0){
		       
		       vecp = c(p, exp(lgamma(n1+a1)+lgamma(n-n1+a2)-lgamma(n+a1+a2))*dnorm(Y[i],0,sd= vecsigma2[t]*sqrt(tau^2+1)))
		       nnew = c(n-1)
		       nmY2 =0
		       nSY2=0
		    }
		  }
		    if (K==0){
		    vecp = c(p, exp(lgamma(n1+a1)+lgamma(n-n1+a2)-lgamma(n+a1+a2))*dnorm(Y[i],0,sd= vecsigma2[t]*sqrt(tau^2+1))) 
		    Knew= 0
		    nnew= c(n-1)
		    Vn2 = as.vector(0)
		   
		    nmY2 = as.vector(0)
		    nSY2= as.vector(0)
		    }
		    
    
		 vecp = vecp/(sum(vecp))
		 
	
      
			Z[i] = sample((0:(Knew+1)),1,prob= vecp)
			
			if (Z[i] < (Knew+1)){
			vecn = nnew
			vecn[Z[i]+1] =vecn[Z[i]+1]+1
		      K = Knew
		      if (Z[i]==0){
		        mY1 = nmY1i
		        SY1 = nSY1i
		    
		     #     if (Knew>0){
		          mY2 = nmY2
		          SY2  = nSY2
		     #      }
		      }  
		      if (0<Z[i]){
		        mY1 = nmY1
		        SY1 = nSY1
		        mY2 = nmY2
		        SY2= nSY2
		        mY2[Z[i]]= nmY2i[Z[i]]
		        SY2[Z[i]]= nSY2i[Z[i]]
		      }
			}
			if (Z[i]== (Knew+1)){
			  vecn = c(nnew,1)
			  mY1 = nmY1
			  SY1 = nSY1
			  if (Knew>0){
			  mY2 = c(nmY2,Y[i])
			  SY2 = c(nSY2,0)}
			  else{
			    mY2= c(Y[i])
			    SY2 = c(0)}
			  K = Knew+1
			}
			
			
		
			
		}
		
		
		n1 = vecn[1]
		BZ[t,]= Z
    vecK[t] = K
    Bn1[t] = n1
   
		
	}
    
    
 #rbind(vecalpha=vecalpha, vecmu = vecmu, vecsigma1=vecsigma1, vecsigma2 =vecsigma2, BZ=BZ, a1 = a1, a2= a2, b1 = b1, b2 = b2, A =A,data1=data1, vecK = vecK, Bn1= Bn1)
 
  #list( kantalpha=quantile(vecalpha,probs=c(.25,.5,.75)),
 	#kantmu=quantile(vecmu,probs=c(.25,.5,.75)),	
	#kantsigma1=quantile(vecsigma1,probs=c(.25,.5,.75)),
	#kantsigma2=quantile(vecsigma2,probs=c(.25,.5,.75)),
	#kantkay=quantile(vecK,probs=c(.25,.5,.75)))
  c(quantile(vecalpha,probs=c(.25,.5,.75)),
 	quantile(vecmu,probs=c(.25,.5,.75)),	
	quantile(vecsigma1,probs=c(.25,.5,.75)),
	quantile(vecsigma2,probs=c(.25,.5,.75)),
	quantile(vecK,probs=c(.25,.5,.75)))


}


##
M=1e2
T=2e4
a1=a2=1/2
A=tau=b1=b2=1

X = posteriorDP3(T, rt(1e3,df=3), a1, a2, tau, b1, b2, A)
write.table(t(X),file = "resultatsH1t7n1e3", col.names=FALSE,row.names=FALSE)

for (i in 2:M){
 X = posteriorDP3(T, rt(1e3,df=3), a1, a2, tau, b1, b2, A)
 write.table(t(X),file = "resultatsH1t7n1e3",append=TRUE, col.names=FALSE,row.names=FALSE)
}
