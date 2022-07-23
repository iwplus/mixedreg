##############################################################
###### Library/Packages yang diperlukan untuk run program ###########
##############################################################

library(pracma)

##############################################################
############ Fungsi tambahan #################################
##############################################################



############## Fungsi menghapus unsur nol di vektor #######

removezero = function(a)
{
 indzero = 0
 ind = 1
 while (ind <= length(a)) #### cari indeks entri tak-nol pertama
 {
	if (a[ind] != 0) 
	{
		indzero = ind
		ind = length(a)+5
	}
 ind = ind + 1
 }
 
 for (i in (indzero+1):length(a)) ### cari semua indeks entri tak-nol
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
### (1) Bagian particle swarm optimization ####
###########################################

###############################################
####Daftar fungsi-fungsi yang berguna #########
###############################################


### fungsi bobot inersia ####

omega <- function(t,tmaks){
	bobinersia = 0.9-(((0.9-0.4)/tmaks)*t);
	bobinersia
}

### fungsi pembangkit matriks R1 dan R2 ##########

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
#### Inisialisasi input PSO* ############
#########################################
### *sebagai contoh saja (di masing-masing bagian sudah ditaruh inisialisasi juga) ################

#N = 5; ### Banyaknya calon solusi awal
#itermaks = 100; ### banyaknya iterasi maksimum
#nvar = 5; ### diisi dengan banyaknya variabel 


################################################
############ Fungsi PSO ########################
################################################

 
pso <- function(N,itermax,nvar,f,y,x,batasvar,Dalpha,degree,nknot) #### bagian awal fungsi pso yang meminimumkan f
{
 npartisidomainvar = 100
 X = matrix(0, nrow=N, ncol = nvar*nknot); ### bangkitkan solusi awal
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

 p = X; ### personal best awal
 v = matrix(0,N,nvar*nknot); ### kecepatan awal


 fitness = matrix(0,1,N);

 
 for (i in 1:N){
	Xtemp = X[i,]
	fitness[i] = f(y,x,Xtemp,Dalpha,degree);
 } 

 sortfitness = sort(fitness, index.return = TRUE);

 g = X[sortfitness$ix[1],]; ### solusi minimum

 #print(fitness)
 #print(X)

 for (iter in 2:itermaks){


 ### update kecepatan particle ###

 for (i in 1:N){
 	v[i,] = (omega(iter,itermaks)*v[i,])+(2*(p[i,]-X[i,])%*%generater12(nvar*nknot))+(2*(g-X[i,])%*%generater12(nvar*nknot));
 }


 ##### update solusi #####

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

} ######################################### akhir dari fungsi pso

#output = pso(N,itermaks,nvar,kuadratik,Dalpha,degree);

#########################################################
#########################################################




#####################################################
############ (2) Bagian Regresi Kernel ##################
#####################################################

#########################################
##### fungsi kernel dan fungsi GCV#######
#########################################

kernel = function(u) #### kernel bi-square
{
 if (abs(u)<=1) {ker = 0.9375*(1-u^2)^2}
 else {ker = 0}
 return(ker)
}

euclidnorm = function(b) #### Norma Euclid 
{
 norma = 0
 for (i in 1:length(b))
 {
	norma = norma + b[i]^2
 }
 return(norma)
}#######################################


nadwatson = function(x0,x1,x,h) ###### Fungsi W pada estimator Nadaraya-Watson
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
##################### Fungsi GCV* #######################################
#########################################################################
########## * keterangan : parameter Dalpha dan degree diset fiks. Parameter tsb Ditambahkan hanya untuk keseragaman input fungsi di PSO ########
#########################################################################

gcv = function(y,x,h,Dalpha,degree)############### Fungsi GCV yang univariat
{
 n = length(y)
 D = matrix(0, ncol = n, nrow = n) ### Inisialisasi matriks D(h)
 
 for (k in 1:n) ##### pengisian entri matriks D(h) dengan fungsi W dari estimator Nadaraya-Watson
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
} ##################################################### Akhir fungsi GCV 
########################################################################################################


#####################################################################################
############# (3) Bagian Spline Truncated ###############################################
#####################################################################################

#### Fungsi truncated ############

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
################### Fungsi GCV untuk Mix* ######################################
################################################################################


gcvmix = function(y,t,knot,Dalpha,degree)
{
 n = nrow(t)
 nmat = ncol(t)
 nknot = length(knot)/nmat
 m = degree + 1 + nknot
 Gk = matrix(0,nrow = n, ncol = m) ###### inisialisasi matriks G(k)
 
 for (indcol in 1:m) ###### menghitung G(k_1)
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
 
 for (indvar in 2:nmat) #### menghitung G(k_2),...,G(k_p) dan digabung jadi G(k)
 {
	Gtemp = matrix(0, nrow = n, ncol = m)
	for (indcol in 1:m) ###### menghitung G(k_1)
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
####### (4) Coba Pengolahan data (eksekusi) #############
################################################################

##### Baca data dan parameter lain yang diperlukan #############

mydata = read.csv('e:2019/PSO plus Mixed Spline and Kernel/datacoba.csv',header=TRUE)
dataindeks = read.csv('e:2019/PSO plus Mixed Spline and Kernel/indeks.csv',header=TRUE)


####### Data (variabel) respon/dependent ###########

mydata = data.frame(mydata)
y = as.matrix(mydata[,3])

### Input indeks variabel yang sudah tercatat di indeks.csv####

dataindeks = data.frame(dataindeks)
indeks = as.matrix(dataindeks)

##### Data (variabel) independent ########################

indeksx = removezero(indeks[,1]) ##### indeks variabel x ada di kolom ke-1 di file indeks.csv 

x = matrix(ncol=length(indeksx),nrow=nrow(y))

for (i in 1:length(indeksx))
{
 for (j in 1:nrow(y))
 {
   x[j,i] = mydata[j,indeksx[i]+3] #### variabel respon dimulai dari kolom ke-4 di file datacoba.csv
 }
}

print(dim(y))
print(dim(x))
y[1]
x[1,1]
x[1,2]

###########################################################################################
########### (4.1) Pencarian bandwidth optimum untuk Regresi Kernel ########################
###########################################################################################

#### Inisialisasi input PSO untuk Regresi Kernel##############


N = 50; ### Banyaknya calon solusi awal
itermaks = 1; ### banyaknya iterasi maksimum
nvar = 1; ### diisi dengan banyaknya variabel yang akan diproses

#### parameter fiks #######
D = matrix(0,nrow(x),nrow(x))
degree = 1
nknot = 1 ### nknot = nvar untuk kasus regresi kernel 
############################

band = matrix(0, ncol = nrow(indeks), nrow = 1) ######### Matriks tempat menyimpan bandwidth
nilaigcv = matrix(0, ncol = nrow(indeks), nrow = 1) ###### Matriks tempat menyimpan GCV

 for (indvar in 1:nrow(indeks)) ##### awal iterasi utama pencarian bandwidth untuk masing-masing variabel 
 {

 ###### batas atas dan bawah masing-masing bandwidth ########

 batasvar = matrix(0, nrow = nvar, ncol = 2)

 for (i in 1:nvar)
 {
	batasvar[i,1] = 0.1
	batasvar[i,2] = 2
 }
 ###############################################################

 xin = x[,indvar]
 bandplusgcv = pso(N,itermaks,nvar,gcv,y,xin,batasvar,D,degree,nknot); ######## pencarian bandwidth optimum dengan PSO (Ingat: nknot = nvar = 1)

 band[indvar] = bandplusgcv[1]
 nilaigcv[indvar] = bandplusgcv[2]
 
 #bandplusgcv[1:length(bandplusgcv)-1]
 #print(gcv(y,xin,bandplusgcv[1:length(bandplusgcv)-1]))

 #print(gcv(y,xin,5))
 #print(kernel((xin[1]-xin[2])/5))
 #print(nadwatson(xin[1],xin[2],xin,5))
 #print(euclidnorm(xin))
} ###### akhir iterasi pencarian bandwidth

print(band)
#print(nilaigcv) #### yang diprint adalah GCV bersama di bagian bawah

######### Pembentukan matriks D(alpha) ###############

n = length(y)
Dalpha = matrix(0, nrow = n, ncol = n)

for (indvar in 1:nrow(indeks))
{
 D = matrix(0, ncol = n, nrow = n) ### Inisialisasi matriks D(h)
 
 for (k in 1:n) ##### pengisian entri matriks D(h) dengan fungsi W dari estimator Nadaraya-Watson
 {
   for (l in 1:n)
   {
	D[k,l] = 1/n*(nadwatson(x[k,indvar],x[l,indvar],x[,indvar],band[indvar]))
   }
 }
 Dalpha = Dalpha + D
}


###### Coba program jalan/tidak :D ############

print(Dalpha)

######## cetak GCV ############

ndata = nrow(Dalpha)
I = diag(ndata)
H = I-Dalpha
HA = (1/ndata)*euclidnorm((H%*%y))^2
HB = ((1/ndata)*sum(diag(H)))^2
gcvkernel = HA/HB
print(gcvkernel)


#knot1 = c(3,4.5,5.5)
#knot2 = c(310,390,425)
#knot = rbind(knot1,knot2)
#h = 3
#gcvmix(y,x,h,Dalpha,1,knot)

#################################################

########################################################################
############ (4.2) Pencarian Knot optimum untuk Spline truncated #############
########################################################################

###### parameter fiks ########

h = 3
###############################

###### Ambil data t ################

indekst = removezero(indeks[,2]) ##### indeks variabel x ada di kolom ke-2 di file indeks.csv

t = matrix(ncol=length(indekst),nrow=nrow(y))

for (i in 1:length(indekst))
{
 for (j in 1:nrow(y))
 {
   t[j,i] = mydata[j,indekst[i]+3] #### variabel respon dimulai dari kolom ke-4 di file datacoba.csv
 }
}

#print(dim(t))
#print(t)

############ Range nilai variabel ########

nvar = length(indekst) ##### banyaknya variabel
batasvar = matrix(0,nrow = nvar, ncol = 2)

for (i in 1:nvar)
{
	batasvar[i,1] = min(t[,i])
	batasvar[i,2] = max(t[,i])
}

#### Inisialisasi input PSO untuk Spline ##############


N = 50; ### Banyaknya calon solusi awal
itermaks = 10; ### banyaknya iterasi maksimum

degree = 3
nknot = 4

knotplusmgcv = pso(N,itermaks,nvar,gcvmix,y,t,batasvar,Dalpha,degree,nknot); ######## pencarian bandwidth optimum dengan PSO 

knotplusmgcv


################# Cetak data prediksi oleh model #############

knotakhir = knotplusmgcv[1:(length(knotplusmgcv)-1)]

n = nrow(t)
nmat = ncol(t)
nknot = length(knotakhir)/nmat
m = degree + 1 + nknot
Gk = matrix(0,nrow = n, ncol = m) ###### inisialisasi matriks G(k)
 
 for (indcol in 1:m) ###### menghitung G(k_1)
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
 
 for (indvar in 2:nmat) #### menghitung G(k_2),...,G(k_p) dan digabung jadi G(k)
 {
	Gtemp = matrix(0, nrow = n, ncol = m)
	for (indcol in 1:m) ###### menghitung G(k_1)
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
 #print(Gk)
 #print(t(Gk)%*%Gk)
 I2 = diag(nrow(Dalpha))
 K = Gk%*%pinv((t(Gk)%*%Gk))%*%t(Gk)%*%(I2-Dalpha)
 M = K + Dalpha

ypred = M%*%y

print(ypred)

#plot(y,type="b", col='red')
#plot(ypred, col = 'blue', add = TRUE)
