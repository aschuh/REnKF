##############################
#--  Variables
#-- you must set these
#-- dim/gridcellsize correspond
#-- to GEOSCHEM run to be used
#-- fileout is where bias ens
#-- is stored.
##############################
nugget.var.perc.prior = 0
var.mult.prior = 0.2^2

origdim.x_2x2.5 =  144
origdim.y_2x2.5 =  91
grid.x_2x2.5 = 2.5
grid.y_2x2.5 = 2

origdim.x_4x5 =  72
origdim.y_4x5 =  46
grid.x_4x5 = 5
grid.y_4x5 = 4

ens.size = 1000
infl.factor.mean = 0.15
infl.factor.var = 0.03^2

obslandbias_H_mean = 0
obslandbias_H_sd = 0.5

obslandbias_M_mean = 0
obslandbias_M_sd = 0.5

obsoceanbias_mean = 0
obsoceanbias_sd = 0.5

b1_offset_mean = 0.5
b1_offset_sd = 0.25

s32_mean = 34
s32_sd = 8

albedo_2_H_mean = 0
albedo_2_H_sd = 5

albedo_2_M_mean = 0
albedo_2_M_sd = 5

dp_cld_mean = 0
dp_cld_sd = 0.1



#fileout = "/user1/aschuh/temp/betas.061113.nc"
#fileout = "~/Desktop/betas.061113.nc"
fileout = "/user1/aschuh/temp/betas.112513.nc"
###############################

#-- Libraries needed

require(fields)
require(geoR)
require(mnormt)
#require(matrixcalc)

#-- COLD START, generation of initial ensemble

grid.pr_2x2.5 = expand.grid(x=seq(-180,180-grid.x_2x2.5,by=grid.x_2x2.5)+0.5*grid.x_2x2.5,
                      y=c(-90+0.5*(0.5*grid.y_2x2.5),seq(-90+0.5*grid.y_2x2.5,
                      90-1.5*grid.y_2x2.5,by=grid.y_2x2.5)+0.5*grid.y_2x2.5,90-0.5*(0.5*grid.y_2x2.5)))

dist.pr_2x2.5 = rdist.earth(grid.pr_2x2.5,miles=FALSE)

ltriang.dist.pr_2x2.5 = lower.tri(dist.pr_2x2.5)


grid.pr_4x5 = expand.grid(x=seq(-180,180-grid.x_4x5,by=grid.x_4x5)+0.5*grid.x_4x5,
                      y=c(-90+0.5*(0.5*grid.y_4x5),seq(-90+0.5*grid.y_4x5,
                      90-1.5*grid.y_4x5,by=grid.y_4x5)+0.5*grid.y_4x5,90-0.5*(0.5*grid.y_4x5)))

dist.pr_4x5 = rdist.earth(grid.pr_4x5,miles=FALSE)

ltriang.dist.pr_4x5 = lower.tri(dist.pr_4x5)


#-- FOR 2x2.5, LAND

VV.land_2x2.5 = varcov.spatial(dists.lowertri=dist.pr_2x2.5[lower.tri(dist.pr_2x2.5)],cov.model="exponential",
                            cov.pars=c(var.mult.prior,800))

ensembs.land_2x2.5 = rmnorm(ens.size*2,mean=rep(0,grid.x_2x2.5*grid.y_2x2.5),varcov=VV.land_2x2.5$varcov)

ensembs.land.gpp_2x2.5 = ensembs.land_2x2.5[1:ens.size,]

ensembs.land.resp_2x2.5 = ensembs.land_2x2.5[(ens.size+1):(ens.size*2),]

rm(VV.land_2x2.5)

#-- FOR 4x5, LAND

VV.land_4x5 = varcov.spatial(dists.lowertri=dist.pr_4x5[lower.tri(dist.pr_4x5)],cov.model="exponential",
                            cov.pars=c(var.mult.prior,800))

ensembs.land_4x5 = rmnorm(ens.size*2,mean=rep(0,grid.x_4x5*grid.y_4x5),varcov=VV.land_4x5$varcov)

ensembs.land.gpp_4x5 = ensembs.land_4x5[1:ens.size,]

ensembs.land.resp_4x5 = ensembs.land_4x5[(ens.size+1):(ens.size*2),]

rm(VV.land_4x5)


#--  Write out ensemble reals to netcdf

#-- FOR 2x2.5
VV.ocean_2x2.5 = varcov.spatial(dists.lowertri = dist.pr_2x2.5[lower.tri(dist.pr_2x2.5)],cov.model="exponential",cov.pars=c(0.1^2,1600))

ensembs.ocean_2x2.5 = rmnorm(ens.size,mean=rep(0,grid.x_2x2.5*grid.y_2x2.5),varcov=VV.ocean_2x2.5$varcov)
                                           
rm(VV.ocean_2x2.5)

#-- FOR 4x5
VV.ocean_4x5 = varcov.spatial(dists.lowertri = dist.pr_4x5[lower.tri(dist.pr_4x5)],cov.model="exponential",cov.pars=c(0.1^2,1600))

ensembs.ocean_4x5 = rmnorm(ens.size,mean=rep(0,grid.x_4x5*grid.y_4x5),varcov=VV.ocean_4x5$varcov)
                                           
rm(VV.ocean_4x5)




#--  Write out ensemble reals to netcdf, 2x2.5

ens = list(ensembs.land.gpp_2x2.5=ensembs.land.gpp_2x2.5,
           ensembs.land.resp_2x2.5=ensembs.land.resp_2x2.5,
           ensembs.ocean_2x2.5=ensembs.ocean_2x2.5)                                        

#--  This sets the first ens member to the "control", just 0 in this case because
#--  GEOSCHEM is expecting the add '1' to this, ie. newGPP = GPP (1 + beta)

ens$ensembs.land.gpp[1,] = 0
ens$ensembs.land.resp[1,] = 0
ens$ensembs.ocean[1,] = 0

#save(ens,file="/user1/aschuh/temp/ens.040912.rda")
#load("/user1/aschuh/temp/ens.050112.rda")

#-- This begins to build the netcdf file
require(ncdf4)

x  = ncdim_def( "lon", "degrees_east", as.double(seq(-180,180-grid.x_2x2.5,by=grid.x_2x2.5))) 
y  = ncdim_def( "lat", "degrees_north", as.double(c(-90+0.25*grid.y_2x2.5,seq((-90+0.25*grid.y_2x2.5)+0.5*grid.y_2x2.5,
                                                         90-0.75*grid.y_2x2.5,by=grid.y_2x2.5),90-0.25*grid.y_2x2.5)))
                                                         
nens =  ncdim_def( "ensemble", units="integer", as.double(seq(1:ens.size)))
                                                         
t  = ncdim_def( "time", "hours since 1900-01-01", 1, unlim=TRUE) 
    	
varlist = list()

varlist[[1]]  <- ncvar_def(name="BETAGPP",units="",
                                 dim=list(x,y,nens,t), missval=NA,prec="float")
                                 
varlist[[2]]  <- ncvar_def(name="BETARESP",units="",
                                 dim=list(x,y,nens,t), missval=NA,prec="float")
                                 
varlist[[3]]  <- ncvar_def(name="BETAOCEAN",units="",
                                 dim=list(x,y,nens,t), missval=NA,prec="float")
                                 
#varlist[[4]]  <- ncvar_def(name="BETAGPP_INFL_MEAN",units="",
#                                 dim=list(x,y,t), missval=NA,prec="float")

#varlist[[5]]  <- ncvar_def(name="BETAGPP_INFL_VAR",units="",
#                                 dim=list(x,y,t), missval=NA,prec="float")

#varlist[[6]]  <- ncvar_def(name="BETARESP_INFL_MEAN",units="",
#                                 dim=list(x,y,t), missval=NA,prec="float")

#varlist[[7]]  <- ncvar_def(name="BETARESP_INFL_VAR",units="",
#                                 dim=list(x,y,t), missval=NA,prec="float")

#varlist[[8]]  <- ncvar_def(name="BETAOCEAN_INFL_MEAN",units="",
#                                 dim=list(x,y,t), missval=NA,prec="float")

#varlist[[9]]  <- ncvar_def(name="BETAOCEAN_INFL_VAR",units="",
#                                 dim=list(x,y,t), missval=NA,prec="float")
                                 
#varlist[[1]] <- ncvar_def(name="OBSB1_OFFSETCOEF",units="",
#                                 dim=list(nens,t), missval=NA,prec="float")                                 

#varlist[[2]] <- ncvar_def(name="OBSS32_COEF",units="",
#                                 dim=list(nens,t), missval=NA,prec="float") 
                                 
#varlist[[3]] <- ncvar_def(name="OBSALBEDO_2_M_COEF",units="",
#                                 dim=list(nens,t), missval=NA,prec="float")                                 
                                              
#varlist[[4]] <- ncvar_def(name="OBSALBEDO_2_H_COEF",units="",
#                                 dim=list(nens,t), missval=NA,prec="float")                                 

#varlist[[5]] <- ncvar_def(name="OBSDP_CLD_COEF",units="",
#                                 dim=list(nens,t), missval=NA,prec="float") 
                                 
#varlist[[6]] <- ncvar_def(name="OBSLANDHBIAS",units="",
#                                 dim=list(nens,t), missval=NA,prec="float")                                 

#varlist[[7]] <- ncvar_def(name="OBSLANDMBIAS",units="",
#                                 dim=list(nens,t), missval=NA,prec="float") 
                                                                  
#varlist[[8]] <- ncvar_def(name="OBSOCEANBIAS",units="",
#                                 dim=list(nens,t), missval=NA,prec="float") 
                                 
#fileout = "/discover/nobackup/aschuh/data/betas/betas.041514_2x2.5.nc"

ncnew <- nc_create(fileout,varlist)

ncvar_put(ncnew, varlist[[1]],t(ens$ensembs.land.gpp))
ncvar_put(ncnew, varlist[[2]],t(ens$ensembs.land.resp))
ncvar_put(ncnew, varlist[[3]],t(ens$ensembs.ocean))
#ncvar_put(ncnew, varlist[[4]],rep(infl.factor.mean,length(x$vals)*length(y$vals)))
#ncvar_put(ncnew, varlist[[5]],rep(infl.factor.var,length(x$vals)*length(y$vals)))
#ncvar_put(ncnew, varlist[[6]],rep(infl.factor.mean,length(x$vals)*length(y$vals)))
#ncvar_put(ncnew, varlist[[7]],rep(infl.factor.var,length(x$vals)*length(y$vals)))
#ncvar_put(ncnew, varlist[[8]],rep(infl.factor.mean,length(x$vals)*length(y$vals)))
#ncvar_put(ncnew, varlist[[9]],rep(infl.factor.var,length(x$vals)*length(y$vals)))

#ncvar_put(ncnew, varlist[[1]],rnorm(ens.size,b1_offset_mean,b1_offset_sd))
#ncvar_put(ncnew, varlist[[2]],rnorm(ens.size,s32_mean,s32_sd))
#ncvar_put(ncnew, varlist[[3]],rnorm(ens.size,albedo_2_M_mean,albedo_2_M_sd))
#ncvar_put(ncnew, varlist[[4]],rnorm(ens.size,albedo_2_H_mean,albedo_2_H_sd))
#ncvar_put(ncnew, varlist[[5]],rnorm(ens.size,dp_cld_mean,dp_cld_sd))

#ncvar_put(ncnew, varlist[[6]],rnorm(ens.size,obslandbias_H_mean,obslandbias_H_sd))
#ncvar_put(ncnew, varlist[[7]],rnorm(ens.size,obslandbias_M_mean,obslandbias_M_sd))
#ncvar_put(ncnew, varlist[[8]],rnorm(ens.size,obsoceanbias_mean,obsoceanbias_sd))

nc_close(ncnew) 





#--  Write out ensemble reals to netcdf, 4x5

ens = list(ensembs.land.gpp_4x5=ensembs.land.gpp_4x5,
           ensembs.land.resp_4x5=ensembs.land.resp_4x5,
           ensembs.ocean_4x5=ensembs.ocean_4x5)            

ens$ensembs.land.gpp[1,] = 0
ens$ensembs.land.resp[1,] = 0
ens$ensembs.ocean[1,] = 0

#-- This begins to build the netcdf file
require(ncdf4)

x  = ncdim_def( "lon", "degrees_east", as.double(seq(-180,180-grid.x_4x5,by=grid.x_4x5))) 
y  = ncdim_def( "lat", "degrees_north", as.double(c(-90+0.25*grid.y_4x5,seq((-90+0.25*grid.y_4x5)+0.5*grid.y_4x5,
                                                         90-0.75*grid.y_4x5,by=grid.y_4x5),90-0.25*grid.y_4x5)))
                                                         
nens =  ncdim_def( "ensemble", units="integer", as.double(seq(1:ens.size)))
                                                         
t  = ncdim_def( "time", "hours since 1900-01-01", 1, unlim=TRUE) 
    	
varlist = list()

varlist[[1]]  <- ncvar_def(name="BETAGPP",units="",
                                 dim=list(x,y,nens,t), missval=NA,prec="float")
                                 
varlist[[2]]  <- ncvar_def(name="BETARESP",units="",
                                 dim=list(x,y,nens,t), missval=NA,prec="float")
                                 
varlist[[3]]  <- ncvar_def(name="BETAOCEAN",units="",
                                 dim=list(x,y,nens,t), missval=NA,prec="float")
                                                 
#fileout = "/discover/nobackup/aschuh/data/betas/betas.041514_4x5.nc"

ncnew <- nc_create(fileout,varlist)

ncvar_put(ncnew, varlist[[1]],t(ens$ensembs.land.gpp))
ncvar_put(ncnew, varlist[[2]],t(ens$ensembs.land.resp))
ncvar_put(ncnew, varlist[[3]],t(ens$ensembs.ocean))

nc_close(ncnew) 