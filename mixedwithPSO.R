##############################################################
###### Library/Packages needed ###########
##############################################################

library(pracma)

##############################################################
############ additional functions #################################
##############################################################



############## Function to delete zero entries in an array -- needed for reading the indices for independet variables #######

removezero = function(a)
{
 indzero = 0
 ind = 1
 while (ind <= length(a)) 
 {
	if (a[ind] != 0) 
	{
		indzero = ind
		ind = length(a)+5
	}
 ind = ind + 1
 }
 
 for (i in (indzero+1):length(a)) 
 {
	if (a[i] != 0)
	{
		indzero = cbind(indzero,i)
	}
 }
 
 b = matrix(0,nrow=1,ncol=length(indzero))
 
 for (j in 1:length(b))
 {
	b[j] = a[indzero[j]]
 }
 
 return(b)
}

###################################################################################################


###########################################
### (1) Particle swarm optimization (PSO) ####
###########################################

###############################################
#### some useful functions #########
###############################################


### Innertia weight function ####

omega <- function(t,tmaks){
	bobinersia = 0.9-(((0.9-0.4)/tmaks)*t);
	bobinersia
}

### A function to generate R1 and R2 ##########

generater12 <- function(nvar){
 entridiag = c(matrix(runif(1*nvar), ncol = nvar));
 if (nvar == 1)
 {
	matriksr = entridiag;
 }
 else
 {
	matriksr = diag(entridiag);
 }
 return(matriksr)
}



#########################################
#### An example for PSO input initialization ############
#########################################


#N = 5; ### The number of initial solutions
#itermaks = 100; ### number of iterations
#nvar = 5; ### number of variables (dependent + independent) 


################################################
############ PSO main function ########################
################################################

 
pso <- function(N,itermax,nvar,f,y,x,batasvar,Dalpha,degree,nknot) #### minimize f
{
 npartisidomainvar = 100
 X = matrix(0, nrow=N, ncol = nvar*nknot); ### generate initial solutions
 for (i in 1:N)
 {
   for (j in 1:nvar)
   {
	rn = sample.int(100,nknot)
	for (k in 1:nknot)
	{
		X[i,((j-1)*nknot)+k] = batasvar[j,1]+((rn[k]-1)*(batasvar[j,2]-batasvar[j,1])/npartisidomainvar)	
	}
     	
     	
   }
 }

 #print(X)

 p = X; ### initial personal best
 v = matrix(0,N,nvar*nknot); ### initial velocity


 fitness = matrix(0,1,N);

 
 for (i in 1:N){
	Xtemp = X[i,]
	fitness[i] = f(y,x,Xtemp,Dalpha,degree);
 } 

 sortfitness = sort(fitness, index.return = TRUE);

 g = X[sortfitness$ix[1],]; ### best solution so far

 #print(fitness)
 #print(X)

 for (iter in 2:itermaks){


 ### particle velocity update  ###

 for (i in 1:N){
 	v[i,] = (omega(iter,itermaks)*v[i,])+(2*(p[i,]-X[i,])%*%generater12(nvar*nknot))+(2*(g-X[i,])%*%generater12(nvar*nknot));
 }


 ##### solution update #####

 for (i in 1:N){
	X[i,] = X[i,] + v[i,];
	Xtemp2 = X[i,]
	if (f(y,x,Xtemp2,Dalpha,degree)< fitness[i]){
		p[i,] = X[i,];
	}
	
	if (f(y,x,X[i,],Dalpha,degree)< f(y,x,g,Dalpha,degree)){
		g = X[i,];
	}
 }

 for (i in 1:N){
	Xtemp3 = X[i,]
	fitness[i] = f(y,x,Xtemp3,Dalpha,degree);
 }

 }

 c(g,f(y,x,g,Dalpha,degree))

} ######################################### the end of pso pso function

#output = pso(N,itermaks,nvar,kuadratik,Dalpha,degree);  ### how to use PSO main function

#########################################################
#########################################################




#####################################################
############ (2) Non-parametric kernel regression ##################
#####################################################

#########################################
##### kernel function and Generalized Cross Validation (GCV) function#######
#########################################

kernel = function(u) #### bi-square kernel 
{
 if (abs(u)<=1) {ker = 0.9375*(1-u^2)^2}
 else {ker = 0}
 return(ker)
}

euclidnorm = function(b) #### Euclidean distance 
{
 norma = 0
 for (i in 1:length(b))
 {
	norma = norma + b[i]^2
 }
 return(norma)
}#######################################


nadwatson = function(x0,x1,x,h) ###### Function W in Nadaraya-Watson estimator
{
 m = length(x)
 K = 0

 for (s in 1:m)
 {
	if (h == 0){Ktemp = 0}
	else {Ktemp = kernel(((x0-x[s])/h))}
	K = K + Ktemp
 }
 
 if (h == 0){Ktemp2 = 0}
 else {Ktemp2 = kernel(((x0-x1)/h))}
 
 W = Ktemp2/K
 return(W)
}################################################################################

#########################################################################
##################### GCV main function* #######################################
#########################################################################
########## *Note : the variables Dalpha and degree fixed. added to make 'appropriate' input for PSO main function ########
#########################################################################

gcv = function(y,x,h,Dalpha,degree)############### GCV for univariate case
{
 n = length(y)
 D = matrix(0, ncol = n, nrow = n) ### Inisialisasi matriks D(h)
 
 for (k in 1:n) ##### D(h) matrix entries with function W of Nadaraya-Watson estimator 
 {
   for (l in 1:n)
   {
	D[k,l] = 1/n*(nadwatson(x[k],x[l],x,h))
   }
 }
 
 I = diag(n)
 HA = I-D
 HB = HA%*%y
 gcvvalue = (1/n*euclidnorm(HB)^2)/((1/n*sum(diag(HA)))^2)
 
 return(gcvvalue)
} ##################################################### the end of GCV main function 
########################################################################################################


#####################################################################################
############# (3) Spline Truncated ###############################################
#####################################################################################

#### truncated function ############

truncate = function(t,k,d)
{
  if (t >= k)
   {
   	trun = (t-k)^d 
   }
  else
   {
	trun = 0
   }
 return(trun)
}

################################################################################
################### GCV function for Mixed non-parametric regression* ######################################
################################################################################


gcvmix = function(y,t,knot,Dalpha,degree)
{
 n = nrow(t)
 nmat = ncol(t)
 nknot = length(knot)/nmat
 m = degree + 1 + nknot
 Gk = matrix(0,nrow = n, ncol = m) 
 
 for (indcol in 1:m) 
 {
	if (indcol <= degree+1)
	{
	  	for (indrow in 1:n)
	  	{
			Gk[indrow,indcol] = t[indrow,1]^(indcol-1)
	  	}
	}
	else
	{
		for (indrow in 1:n)
		{
			Gk[indrow,indcol] = truncate(t[indrow,1],knot[indcol-(degree+1)],degree)
		}	
	}
 }	
 
 for (indvar in 2:nmat) 
 {
	Gtemp = matrix(0, nrow = n, ncol = m)
	for (indcol in 1:m) 
 	{
		if (indcol <= (degree+1))
		{
	  		for (indrow in 1:n)
	  		{
				
				Gtemp[indrow,indcol] = t[indrow,indvar]^(indcol-1)
	  		}
		}
		else
		{
			for (indrow in 1:n)
			{
				Gtemp[indrow,indcol] = truncate(t[indrow,indvar],knot[((indvar-1)*nknot)+(indcol-(degree+1))],degree)
			}	
		}
 	}	

 	Gk = cbind(Gk,Gtemp)	
 }
 #print(Gk)
 #print(t(Gk)%*%Gk)
 I = diag(nrow(Dalpha))
 K = Gk%*%pinv((t(Gk)%*%Gk))%*%t(Gk)%*%(I-Dalpha)
 M = K + Dalpha
 H = I-M
 HA = 1/n*(euclidnorm(H%*%y))^2
 HB = (1/n*sum(diag(H)))^2
 gcvm = HA/HB
 
 return(gcvm)
}


################################################################
####### (4) Computation #############
################################################################

##### Read the data and indices for kernel and spline truncated independent variables #############

mydata = read.csv('datacoba.csv',header=TRUE)
dataindeks = read.csv('indeks.csv',header=TRUE)


####### dependent variable ###########

mydata = data.frame(mydata)
y = as.matrix(mydata[,3])

### Independent variables ####

dataindeks = data.frame(dataindeks)
indeks = as.matrix(dataindeks)



indeksx = removezero(indeks[,1]) ##### indices of independent variables (for kernel regression) start from the first column in indeks.csv 

x = matrix(ncol=length(indeksx),nrow=nrow(y))

for (i in 1:length(indeksx))
{
 for (j in 1:nrow(y))
 {
   x[j,i] = mydata[j,indeksx[i]+3] #### the indices of independent variables (for kernel and spline regressions) start from the 4-th column in datacoba.csv
 }
}

print(dim(y))
print(dim(x))
y[1]
x[1,1]
x[1,2]

###########################################################################################
########### (4.1) Optimal bandwidth search for Kernel Regression ########################
###########################################################################################

#### PSO main function input initialization for Kernel regression##############


N = 50; ### the number of initial solutions
itermaks = 1; ### number of iterations
nvar = 1; ### number of variables (dependent + independent)

#### fixed parameter #######
D = matrix(0,nrow(x),nrow(x))
degree = 1
nknot = 1 ### nknot = nvar for kernel regression 
############################

band = matrix(0, ncol = nrow(indeks), nrow = 1) ######### bandwidths matrix
nilaigcv = matrix(0, ncol = nrow(indeks), nrow = 1) ###### GCV values matrix

 for (indvar in 1:nrow(indeks)) 
 {

 ###### upper and lower bounds for each bandwidth ########

 batasvar = matrix(0, nrow = nvar, ncol = 2)

 for (i in 1:nvar)
 {
	batasvar[i,1] = 0.1
	batasvar[i,2] = 2
 }
 ###############################################################

 xin = x[,indvar]
 bandplusgcv = pso(N,itermaks,nvar,gcv,y,xin,batasvar,D,degree,nknot); ######## PSO search (note: nknot = nvar = 1 for kernel regression)

 band[indvar] = bandplusgcv[1]
 nilaigcv[indvar] = bandplusgcv[2]
 
} ###### the end of bandwidth search

print(band)

######### D(alpha) matrix ###############

n = length(y)
Dalpha = matrix(0, nrow = n, ncol = n)

for (indvar in 1:nrow(indeks))
{
 D = matrix(0, ncol = n, nrow = n) 
 
 for (k in 1:n) 
 {
   for (l in 1:n)
   {
	D[k,l] = 1/n*(nadwatson(x[k,indvar],x[l,indvar],x[,indvar],band[indvar]))
   }
 }
 Dalpha = Dalpha + D
}


######## print GCV values ############

ndata = nrow(Dalpha)
I = diag(ndata)
H = I-Dalpha
HA = (1/ndata)*euclidnorm((H%*%y))^2
HB = ((1/ndata)*sum(diag(H)))^2
gcvkernel = HA/HB
print(gcvkernel)



#################################################

########################################################################
############ (4.2) optimal Knots for Spline truncated regression #############
########################################################################

###### fixed parameters for 'appropriate' input in PSO main function ########

h = 3 ### h = number of independent variables in the case of spline truncated regression
###############################

###### Indices for independent variables ################

indekst = removezero(indeks[,2]) ##### indices of independent variables (for spline truncated regression) start from the 2nd column in indeks.csv

t = matrix(ncol=length(indekst),nrow=nrow(y))

for (i in 1:length(indekst))
{
 for (j in 1:nrow(y))
 {
   t[j,i] = mydata[j,indekst[i]+3] #### the indices of independent variables (for kernel and spline regressions) start from the 4-th column datacoba.csv
 }
}

#print(dim(t))
#print(t)

############ Collects lower and upper bounds for each independent variable ########

nvar = length(indekst) 
batasvar = matrix(0,nrow = nvar, ncol = 2)

for (i in 1:nvar)
{
	batasvar[i,1] = min(t[,i])
	batasvar[i,2] = max(t[,i])
}

#### PSO main function initial input for Spline ##############


N = 50; ### the number of initial solutions
itermaks = 10; ### number of iterations

degree = 3 ### the degree of polynomial in spline
nknot = 4 ## number of knots

knotplusmgcv = pso(N,itermaks,nvar,gcvmix,y,t,batasvar,Dalpha,degree,nknot); ######## PSO search 

knotplusmgcv


################# print predicted data #############

knotakhir = knotplusmgcv[1:(length(knotplusmgcv)-1)]

n = nrow(t)
nmat = ncol(t)
nknot = length(knotakhir)/nmat
m = degree + 1 + nknot
Gk = matrix(0,nrow = n, ncol = m) 
 
 for (indcol in 1:m) 
 {
	if (indcol <= degree+1)
	{
	  	for (indrow in 1:n)
	  	{
			Gk[indrow,indcol] = t[indrow,1]^(indcol-1)
	  	}
	}
	else
	{
		for (indrow in 1:n)
		{
			Gk[indrow,indcol] = truncate(t[indrow,1],knotakhir[indcol-(degree+1)],degree)
		}	
	}
 }	
 
 for (indvar in 2:nmat) 
 {
	Gtemp = matrix(0, nrow = n, ncol = m)
	for (indcol in 1:m) 
 	{
		if (indcol <= (degree+1))
		{
	  		for (indrow in 1:n)
	  		{
				
				Gtemp[indrow,indcol] = t[indrow,indvar]^(indcol-1)
	  		}
		}
		else
		{
			for (indrow in 1:n)
			{
				Gtemp[indrow,indcol] = truncate(t[indrow,indvar],knotakhir[((indvar-1)*nknot)+(indcol-(degree+1))],degree)
			}	
		}
 	}	

 	Gk = cbind(Gk,Gtemp)	
 }
 
 I2 = diag(nrow(Dalpha))
 K = Gk%*%pinv((t(Gk)%*%Gk))%*%t(Gk)%*%(I2-Dalpha)
 M = K + Dalpha

ypred = M%*%y

print(ypred)

#####################################
#### plot predicted vs original #####
####################################

plot(y,type="b", col='red')
plot(ypred, col = 'blue', add = TRUE)
