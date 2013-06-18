#####################################################
#--  This is where we would start from a matrix
#--  of observations and model output
#--  ens_matrix and obs_vector need to come in as ppm
#####################################################

optimize_betas = function(betas_file,Rdiag_vector,ens_matrix,obs_vector,method=2,add_error=FALSE,localize=FALSE,
                          estimate_land_ocean_bias=FALSE,diags=FALSE)
{

#-- Subtract one for "obs" in last column, no longer
nens = dim(ens_matrix)[2] #- 1

mat.sqrt = function(a)
{
 a.eig <- eigen(a)
 a.sqrt <- a.eig$vectors %*% diag(sqrt(a.eig$values)) %*% solve(a.eig$vectors)
 return(a.sqrt)
}
	
out = as.data.frame(ens_matrix)

#-- Convert data frame columns from factors back to numeric

#-- for stations data
#out[2:204] = apply(out[2:204],2,FUN=function(x){as.numeric(as.character(x))})
#out[,4:204] = out[,4:204] * 10^10

#-- For ASCENDS
out = apply(out,2,FUN=function(x){as.numeric(as.character(x))})

#-- ***********KEEP EYE ON THIS, MIGHT NEED FOR ASCENDS**************
#out = out * 10^6

#-- This is for pseudo experiment, adds the assumed 1 ppm sd error
if(add_error)
{
	rnormpert =  rnorm(dim(out)[1],mean=0,sd=sqrt(Rdiag_vector))
    #out[,nens+1] = out[,nens+1] + rnormpert
    obs_vector = obs_vector + rnormpert 
}

#-- Convert REAL NOAA data to numeric from factor
#out[,"NOAA"] = as.numeric(as.character(out[,"NOAA"]))

#-- Pull land/water/ice mask data (from GEOSCHEM) into matrix

OCEAN.MASK = !is.na(rotate270.matrix(read.fwf("/discover/nobackup/aschuh/run/ENSCODE/Regions_ocean.dat",widths=rep(1,144))))
LAND.MASK = !is.na(rotate270.matrix(read.fwf("/discover/nobackup/aschuh/run/ENSCODE/Regions_land.dat",widths=rep(1,144))))

#-- Not sure, the GEOS masks seem to have the land mask w/ two extra rows which
#-- appear to be at north pole, have to check this one out, probably bug in their routines

LAND.MASK = LAND.MASK[,3:93]

#-- Pull betas used for the ensemble runs (eventually needs to be tagged w/ cycle #)

require(ncdf4)

if(class(betas_file)=="character")
{
 fil = nc_open(betas_file)

 betas_gpp = ncvar_get(fil,"BETAGPP")
 betas_resp = ncvar_get(fil,"BETARESP")
 betas_ocean = ncvar_get(fil,"BETAOCEAN")

 if(estimate_land_ocean_bias)
 {
    landbias_obs = ncvar_get(fil,"OBSLANDBIAS")
    oceanbias_obs = ncvar_get(fil,"OBSOCEANBIAS")
 }

 nc_close(fil)

 betas_vect_gpp = apply(betas_gpp,c(3),FUN=function(x){x[LAND.MASK]})  #betas_gpp[LAND.MASK] 
 betas_vect_resp = apply(betas_resp,c(3),FUN=function(x){x[LAND.MASK]})
 betas_vect_ocean = apply(betas_ocean,c(3),FUN=function(x){x[OCEAN.MASK]})

 #-- BETA vectorized, GPP FIRST!!
 betas_vect = rbind(betas_vect_gpp,betas_vect_resp,betas_vect_ocean)

 if(estimate_land_ocean_bias)
 {
    betas_vect = rbind(betas_vect,landbias_obs,oceanbias_obs)
 }

 rm(betas_gpp)
 rm(betas_resp)
 rm(betas_ocean)

 rm(betas_vect_gpp)
 rm(betas_vect_resp)
 rm(betas_vect_ocean)
}

if(class(betas_file)=="matrix")
{
  betas_vect = betas_file	
}

###################################################################
#-- EnKF equations, Tippett et al 2003, "direct method"
###################################################################

if(method==2)
{
	
chi_sq_fun = function(alpha,obs_ind_sq)
{
	Y2 = 1/sqrt(sum(obs_ind_sq))*mat.sqrt(solve(HA[obs_ind_sq,] %*% t(HA[obs_ind_sq,]) + alpha*R[obs_ind_sq,obs_ind_sq])) %*% as.matrix(  obs_vector[obs_ind_sq]  - out[obs_ind_sq,1])
	return(abs(sum(diag(Y2 %*% t(Y2)))-1))
}

#require(SparseM)

obs_ind = rep(TRUE,dim(out)[1])
#obs_ind = rep(FALSE,dim(out)[1])
#obs_ind[seq(1,length(obs_ind),by=20)] = TRUE
#ind = !is.na(out$NOAA)

nobs = dim(out[obs_ind,])[1]

#-- Just diagonal matrices right now
R = diag(Rdiag_vector[obs_ind])
Rinv = diag(Rdiag_vector[obs_ind]^-1)

#R = as(Rdiag_vector[obs_ind],"matrix.diag.csr")
#R.inv.sq <- as(sqrt(Rdiag_vector[obs_ind])^-1,"matrix.diag.csr")  #vars for obs

#-- These are standard residuals, not normalized by ensemble number

A =  1/sqrt(nens-1)*(betas_vect[,2:nens] - betas_vect[,1])


if(estimate_land_ocean_bias)
 {
    #-- Not sure yet
 }

#-- HA is forecast observation residuals, normalized by ensemble size

HA = 1/sqrt(nens-1) * apply(out[obs_ind,2:nens],2,FUN=function(x){x - out[obs_ind,1]})

#HA[fulldat$landmask,] = HA[fulldat$landmask,] + matrix(rep(bb[-1,1],sum(fulldat$landmask)),ncol=length(bb[-1,1]),byrow=TRUE)
#HA[fulldat$landmask==0,] = HA[fulldat$landmask==0,] + matrix(rep(bb[-1,2],sum(fulldat$landmask==0)),ncol=length(bb[-1,2]),byrow=TRUE)

#-- Provide subset of obs to calc chi-sq stats on

obs_ind_sq = obs_ind

if(sum(obs_ind_sq)>1000)
	{
		trs = grep(TRUE,obs_ind_sq)
		#-- This turns most (except about 1000) TRUE's to FALSE's, 'slicing along the time dim'
		fls = trs[!trs %in% trs[seq(1,length(trs),by=floor(length(trs)/1000))]]
		obs_ind_sq[fls] = FALSE
	}
	
alpha = optimize(f=chi_sq_fun,interval=c(0,100),tol=0.01,obs_ind_sq=obs_ind_sq)$minimum

#-- We just adjust Rinv because nothing else is currently being used (R related), past this pt
Rinv = alpha^-1 * Rinv

#-- Create diagnostic infl matrix
S0 = diag(Rinv) * apply(HA,1,FUN=function(x){sum(x^2)})

#-- This is what S should be, we'll directly find S^-1 instead
#S = (1/(nens-2))*HA %*% t(HA) + R
cent = solve(diag(rep(1,dim(HA)[2])) + t(HA) %*% Rinv %*% HA )

Sinv = Rinv %*% HA %*% cent %*%
                   t(HA) %*% Rinv

Sinv2 = as.matrix(Rinv) - as.matrix(Sinv)

Z = Sinv2 %*% HA

TTmat = diag(rep(1,dim(HA)[2])) - t(HA) %*% Z

#-- Really have two choices here, going w/ SVD for now

Tmat = mat.sqrt(TTmat)

#Tmat2 = chol(TTmat)


#-- Calc posterior mean, ENS0001 is always the mean/prior 
#X_postmean = betas_vect[,1] + A %*% t(HA) %*% Sinv2 %*% as.matrix(  out[obs_ind,nens+1]  - out[obs_ind,1])
X_postmean = betas_vect[,1] + A %*% t(HA) %*% Sinv2 %*% as.matrix(  obs_vector[obs_ind]  - out[obs_ind,1])

#-- For NOAA in situ
#X_postmean = betas_vect[,1] + 1/198*A %*% t(HA) %*% solve(P) %*% as.matrix( ( out$PSEUDO[obs_ind] + rnormpert ) - out[obs_ind,4])

#-- Using SVD sqrt instead of chol sqrt, seems a bit more stable
#-- also, we need to be very careful about sqrt() scalar in front
#-- which removes "normalization" in definition of A
X_postsd = sqrt(nens-1) * A %*% Tmat

#-- Create posterior ensemble
X_post = cbind(X_postmean,matrix(rep(X_postmean,nens-1),ncol=nens-1,byrow=F) + X_postsd)

if(localize)
{
######################################################################################
#-- Now, we mask by low variance spots, localization based on priorvar/postvar ratio
#-- Then we rerun same equations, w/ KGAIN "tapered" to only high infl areas
######################################################################################

#-- Determine cutoff for LAND localization, based on RESP footprint (as opposed to GPP)
postsd_land = apply(X_postsd[(sum(LAND.MASK)+1):(sum(LAND.MASK)*2),],1,sd)
prsd_land = apply(A[(sum(LAND.MASK)+1):(sum(LAND.MASK)*2),],1,sd)
pr_post_land_ratio = prsd_land/postsd_land
quant.land_upper40 = quantile(pr_post_land_ratio,p=0.6)

#-- Determine cutoff for OCEAN localization
postsd_ocean = apply(X_postsd[(sum(LAND.MASK)*2+1):(dim(A)[1]),],1,sd)
prsd_ocean = apply(A[(sum(LAND.MASK)*2+1):(dim(A)[1]),],1,sd)
pr_post_ocean_ratio = prsd_ocean/postsd_ocean
quant.ocean_upper10 = quantile(pr_post_ocean_ratio,p=0.9)


land_state_ind = pr_post_land_ratio > quant.land_upper40
ocean_state_ind = pr_post_ocean_ratio > quant.ocean_upper10

#-- Mask state based on beta vector of GPP,RESP,OCEANNEE
state_ind = rep(0,sum(LAND.MASK)*2+sum(OCEAN.MASK))
state_ind[c(land_state_ind,land_state_ind,ocean_state_ind)] = 1


#-- Repeat MEAN calc w/ heaviside filtering based on covariance ratio
#-- Calc posterior mean, ENS0001 is always the mean/prior 

X_postmean_adj = A %*% t(HA) %*% Sinv2 
X_postmean_adj = apply(X_postmean_adj,2,FUN=function(x){x*state_ind})

X_postmean = betas_vect[,1] + X_postmean_adj %*% as.matrix(  obs_vector[obs_ind]  - out[obs_ind,1])

#-------------------------------------------------------------

#-- Essentially this updates the perturbations for localized portions and replaces the 
#-- the non-localized portions with the "old" X_postsd
X_postsd_new = sqrt(nens-1) * A %*% Tmat
X_postsd_new[!state_ind,] = sqrt(nens-1) * A[!state_ind,]
X_postsd = X_postsd_new

#-- Create posterior ensemble
X_post = cbind(X_postmean,matrix(rep(X_postmean,nens-1),ncol=nens-1,byrow=F) + X_postsd)

}


}
##################################################
#-- EnKF equations, (linear operators), METHOD 2
##################################################

if(method==1)
{
#require(SparseM)

obs_ind = rep(TRUE,dim(out)[1])
#obs_ind = rep(FALSE,dim(out)[1])
#obs_ind[seq(1,length(obs_ind),by=20)] = TRUE
#ind = !is.na(out$NOAA)

nobs = dim(out[obs_ind,])[1]

#-- Just diagonal matrices right now
R = diag(Rdiag_vector[obs_ind])
Rinv = diag(Rdiag_vector[obs_ind]^-1)

#R = as(Rdiag_vector[obs_ind],"matrix.diag.csr")
#R.inv <- as(Rdiag_vector[obs_ind]^-1,"matrix.diag.csr")  #vars for obs

A = betas_vect[,2:nens] - betas_vect[,1]

HA = apply(out[obs_ind,2:nens],2,FUN=function(x){x - out[obs_ind,1]})

P = (1/(nens-2))*HA %*% t(HA) + R

#-- Calc posterior mean, ENS0001 is always the mean/prior 
#X_postmean = betas_vect[,1] + 1/(nens-2)*A %*% t(HA) %*% solve(P) %*% as.matrix(  out[obs_ind,nens+1]  - out[obs_ind,1])
X_postmean  = betas_vect[,1] + 1/(nens-2)*A %*% t(HA) %*% solve(P) %*% as.matrix(  obs_vector[obs_ind]  - out[obs_ind,1])

#-- For NOAA in situ
#X_postmean = betas_vect[,1] + 1/198*A %*% t(HA) %*% solve(P) %*% as.matrix( ( out$PSEUDO[obs_ind] + rnormpert ) - out[obs_ind,4])

#-- Calc posterior sd, do we want 1/(nens-2) here?, can't find support
temp = 1/(nens-2) * t(HA) %*% Rinv %*% HA
X_postsd = A %*% mat.sqrt(solve(diag(rep(1,dim(temp)[1])) + temp))

#-- Create posterior ensemble
X_post = matrix(rep(X_postmean,nens-1),ncol=nens-1,byrow=F) + X_postsd

if(diags)
{
###########################################################################################
#-- CHI SQ DIAGNOSTIC  #################################################################### 
###########################################################################################
#zis = apply(out[,2:nens],2,FUN=function(x){x - out[,1]})
#zzt = zis %*% t(zis)
#zzt.inv = solve(diag(rep(1,dim(zzt)[1]))+zzt)
#chisq = 1/(dim(out)[1])*t(R.inv.sq %*% (out[,nens+1] - out[,1] )) %*% zzt.inv %*% 
#                             (R.inv.sq %*% (out[,nens+1] - out[,1] ) )
}

if(localize)
{
######################################################################################
#-- Now, we mask by low variance spots, localization based on priorvar/postvar ratio
#-- Then we rerun same equations, w/ KGAIN "tapered" to only high infl areas
######################################################################################

#-- Determine cutoff for LAND localization, based on RESP footprint (as opposed to GPP)
postsd_land = apply(X_postsd[(sum(LAND.MASK)+1):(sum(LAND.MASK)*2),],1,sd)
prsd_land = apply(A[(sum(LAND.MASK)+1):(sum(LAND.MASK)*2),],1,sd)
pr_post_land_ratio = prsd_land/postsd_land
quant.land_upper40 = quantile(pr_post_land_ratio,p=0.6)

#-- Determine cutoff for OCEAN localization
postsd_ocean = apply(X_postsd[(sum(LAND.MASK)*2+1):(dim(A)[1]),],1,sd)
prsd_ocean = apply(A[(sum(LAND.MASK)*2+1):(dim(A)[1]),],1,sd)
pr_post_ocean_ratio = prsd_ocean/postsd_ocean
quant.ocean_upper10 = quantile(pr_post_ocean_ratio,p=0.9)


land_state_ind = pr_post_land_ratio > quant.land_upper40
ocean_state_ind = pr_post_ocean_ratio > quant.ocean_upper10

#-- Mask state based on beta vector of GPP,RESP,OCEANNEE
state_ind = rep(0,sum(LAND.MASK)*2+sum(OCEAN.MASK))
state_ind[c(land_state_ind,land_state_ind,ocean_state_ind)] = 1


#-- Repeat MEAN calc w/ heaviside filtering based on covariance ratio
#-- Calc posterior mean, ENS0001 is always the mean/prior 

KGAIN = 1/198*A %*% t(HA) %*% solve(P)
splits = seq(1,dim(KGAIN)[2],length=10)
splits = round(splits)

for(i in 1:(length(splits)-1))
  {
	KGAIN[,(splits[i]:splits[i+1])] = apply(KGAIN[,(splits[i]:splits[i+1])],2,FUN=function(x){x*state_ind})
  }
  
X_postmean = betas_vect[,1] + KGAIN %*% as.matrix( ( out$PSEUDO[obs_ind] + rnormpert ) - out[obs_ind,4])


#-- Need to inflate betas and "adjust the X_postsd STILL for localization effect !!!!!!!!!!!!!!!!!!!

 }

}
########################
#-End EnKF eqs METHOD 1
########################

####################################
#-- Linear Operators, Zupanski style
####################################

if(method==3){
#-- Via Zupanski et al 2005 
require(SparseM)

#-- Print initial mismatch
print(paste("initial Model/Data mismatch, Mean:",mean(out[,nens+1] - out[,1])," SD:",sd(out[,nens+1] - out[,1])))
#-- Same as below but easier to see
#zis = apply(out[,2:500],2,FUN=function(x){(out[,501] - out[,1]) - (out[,501] - x)})
zis = apply(out[,2:nens],2,FUN=function(x){(x - out[,1])})

C = t(zis) %*% zis

#-- Cost function, equation (11) from Zupanski et al 2005
#--
#--  cost(eta) = A*eta - B*zis*R^0.5 X (y - H(xb) - H())
#--

A = solve(diag(rep(1,nens-1)) + C/(nens-2))
B = mat.sqrt(A)
#R = as(dim(zis)[1],"matrix.diag.csr")
R = as(Rdiag_vector,"matrix.diag.csr")
Rinv = as(Rdiag_vector^-1,"matrix.diag.csr")
Rinv.sq	 = as(sqrt(Rdiag_vector)^-1,"matrix.diag.csr")
resids = matrix(out[,nens+1] - out[,1],ncol=1)
part1 = t(zis) %*% R
adj = (-B %*% part1)  %*% resids

###################################################################
#-- MISSING NEGATIVE SOMEWHERE, CAN'T FIND, ADDING HERE BUT WE NEED 
#-- TO FIND THIS, AESCHUH

#adj = - adj

####################################################################

Pf = 1/sqrt(nens-2) * (betas_vect[,2:nens] - betas_vect[,1])
part2 = Pf %*% t(A)

xadj = part2 %*% as.matrix(adj)

#-- Calculate the standard deviation for the ensemble 

xperts = t(zis) %*% Rinv

xperts = xperts %*% zis

xperts = solve(as.matrix(xperts + diag(dim(xperts)[1])))

#xperts = as.matrix(xperts)

xperts = mat.sqrt(xperts)

betamat = betas_vect[,2:nens] - betas_vect[,1]

xperts  = betamat %*% xperts

#-- Create posterior ensemble
X_post = matrix(rep(as.vector(xadj),nens-1),ncol=nens-1,byrow=F) + xperts

if(diags)
{
###########################################################################################
#-- CHI SQ DIAGNOSTIC  #################################################################### 
###########################################################################################
#zis = apply(out[,2:nens],2,FUN=function(x){x - out[,1]})
#zzt = zis %*% t(zis)
#zzt.inv = solve(diag(rep(1,dim(zzt)[1]))+zzt)
#chisq = 1/(dim(out)[1])*t(R.inv.sq %*% (out[,nens+1] - out[,1] )) %*% zzt.inv %*% 
#                             (R.inv.sq %*% (out[,nens+1] - out[,1] ) )
 
Cdecomp = eigen(C)
V = Cdecomp$vectors
lambda = diag(Cdecomp$values)

Psi = - solve(diag(rep(1,dim(lambda)[1])) + lambda)
sigma = matrix(0,ncol=dim(lambda)[2],nrow=dim(lambda)[1])
gamma = matrix(0,ncol=dim(lambda)[2],nrow=dim(lambda)[1])
N = 3

for(k in 1:N)
{
  sigma = 0.5*(sigma + Psi - gamma - Psi %*% lambda %*% gamma)	
  gamma = sigma %*% solve( diag(rep(1,dim(lambda)[1])) + lambda %*% sigma)
} 

#-- Ginv.sqrt = (I + Z %*% V %*% sigma %*% t(V) %*% t(Z) )

part1 = 1/(dim(out)[1])*t(Rinv.sq %*% (out[,nens+1] - out[,1] )) %*% 
          (Rinv.sq %*% (out[,nens+1] - out[,1] ))
#part2 = 1/(dim(out)[1])*t(Rinv.sq %*% (out[,nens+1] - out[,1] )) %*% zis %*% V
#part3 = t(zis) %*% (Rinv.sq %*% (out[,nens+1] - out[,1] ))
#part3 = Psi %*% t(V) %*% part3
#chisq = part1 + part2 %*% part3


part2 = 1/(dim(out)[1])*t(Rinv.sq %*% (out[,nens+1] - out[,1] )) %*% zis %*% V %*% Psi
part3 = t(V) %*% t(zis) %*% (Rinv.sq %*% (out[,nens+1] - out[,1] ))
chisq = part1 + part2 %*% part3

#chisq = 1/(dim(out)[1])*t(R.inv.sq %*% (out[,nens+1] - out[,1] )) %*% Ginv.sqrt %*% 
#                             (R.inv.sq %*% (out[,nens+1] - out[,1] ) )
                             
print(paste("chisq: ",chisq) )                             
###########################################################################################

###########################################################################################
#-- Normality checks, 3 common tests, shapiro pretty conservative for rejection
###########################################################################################

#  1/(dim(out)[1])* (I + ZVsigmat(V)%*%t(Z))     (Rinv.sq %*% (out[,nens+1] - out[,1] ))

p1 = 1/sqrt(dim(out)[1])*(Rinv.sq %*% (out[,nens+1] - out[,1] )) 
p2 = 1/sqrt(dim(out)[1])*zis %*% V %*% sigma
p2 = p2 %*% t(V)
p3 = t(zis) %*% (Rinv.sq %*% (out[,nens+1] - out[,1] ))
innov = p2 %*% p3

require(nortest)

cvm.test(innov)
ad.test(innov)
#shapiro.test(innov)

print(paste("Innovations, Mean:",mean(out[,nens+1] - out[,1])," SD:",sd(out[,nens+1] - out[,1])))
###########################################################################################

rm(zzt);rm(zzt.inv);

}
rm(betamat);rm(part1);rm(adj);rm(Pf);rm(C)
rm(zis);
rm(A);rm(B);rm(R);rm(Rinv);rm(resids);
}

################################
#-  End of Zupanski, method 4
################################

################################
#-  This method unfinished
################################
if(method==4)
{
	
	obs_ind = rep(TRUE,dim(out)[1])
	
    #-- 
	A = betas_vect[,2:500] - betas_vect[,1]
	
	HA = as.matrix(out[obs_ind,2:500]) - (1/499)*(as.matrix(out[obs_ind,2:500]) %*% matrix(rep(1,499),ncol=1) ) %*% matrix(rep(1,499),nrow=1) 
	
	tHA = t(HA)
	
	#tHA = apply(tHA,1,FUN=function(x){x*Rdiag_vector*1/498})
	
	part1 = tHA %*% HA 
	part2 = part1 + diag(rep(1,dim(part1)[1]))
	part3 = solve(part2)
	part4 = 1/498 * HA %*% part3
	part5 = part4 %*% tHA
}

return(list(X_post=X_post,S0=S0))
}
