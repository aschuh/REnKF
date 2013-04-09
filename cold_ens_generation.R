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
origdim.x =  144
origdim.y =  91
grid.x = 2.5
grid.y = 2
ens.size = 1000
infl.factor.mean = 0.15
infl.factor.var = 0.03^2
fileout = "/home/aschuh/betas.040913.nc"
###############################

#-- Libraries needed

require(fields)
require(geoR)
require(mnormt)
#require(matrixcalc)

#-- COLD START, generation of initial ensemble

grid.pr = expand.grid(x=seq(-180,180-grid.x,by=grid.x)+0.5*grid.x,
                      y=c(-90+0.5*(0.5*grid.y),seq(-90+0.5*grid.y,
                      90-1.5*grid.y,by=grid.y)+0.5*grid.y,90-0.5*(0.5*grid.y)))

dist.pr = rdist.earth(grid.pr,miles=FALSE)

ltriang.dist.pr = lower.tri(dist.pr)

VV.land = varcov.spatial(dists.lowertri=dist.pr[lower.tri(dist.pr)],cov.model="exponential",
                            cov.pars=c(var.mult.prior,800))

ensembs.land = rmnorm(ens.size*2,mean=rep(0,grid.x*grid.y),varcov=VV.land$varcov)

ensembs.land.gpp = ensembs.land[1:ens.size,]

ensembs.land.resp = ensembs.land[(ens.size+1):(ens.size*2),]

rm(VV.land)

#--  Write out ensemble reals to netcdf

VV.ocean = varcov.spatial(dists.lowertri = dist.pr[lower.tri(dist.pr)],cov.model="exponential",cov.pars=c(0.1^2,1600))

ensembs.ocean = rmnorm(ens.size,mean=rep(0,grid.x*grid.y),varcov=VV.ocean$varcov)
                                           
rm(VV.ocean)

#--  Write out ensemble reals to netcdf

ens = list(ensembs.land.gpp=ensembs.land.gpp,ensembs.land.resp=ensembs.land.resp,ensembs.ocean=ensembs.ocean)                                        

#--  This sets the first ens member to the "control", just 0 in this case because
#--  GEOSCHEM is expecting the add '1' to this, ie. newGPP = GPP (1 + beta)

ens$ensembs.land.gpp[1,] = 0
ens$ensembs.land.resp[1,] = 0
ens$ensembs.ocean[1,] = 0

#save(ens,file="/user1/aschuh/temp/ens.040912.rda")
#load("/user1/aschuh/temp/ens.050112.rda")

#-- This begins to build the netcdf file
require(ncdf4)

x  = ncdim_def( "lon", "degrees_east", as.double(seq(-180,180-grid.x,by=grid.x))) 
y  = ncdim_def( "lat", "degrees_north", as.double(c(-90+0.25*grid.y,seq((-90+0.25*grid.y)+0.5*grid.y,
                                                         90-0.75*grid.y,by=grid.y),90-0.25*grid.y)))
                                                         
nens =  ncdim_def( "ensemble", units="integer", as.double(seq(1:ens.size)))
                                                         
t  = ncdim_def( "time", "hours since 1900-01-01", 1, unlim=TRUE) 
    	
varlist = list()

varlist[[1]]  <- ncvar_def(name="BETAGPP",units="",
                                 dim=list(x,y,nens,t), missval=NA,prec="float")
                                 
varlist[[2]]  <- ncvar_def(name="BETARESP",units="",
                                 dim=list(x,y,nens,t), missval=NA,prec="float")
                                 
varlist[[3]]  <- ncvar_def(name="BETAOCEAN",units="",
                                 dim=list(x,y,nens,t), missval=NA,prec="float")
                                 
varlist[[4]]  <- ncvar_def(name="BETAGPP_INFL_MEAN",units="",
                                 dim=list(x,y,t), missval=NA,prec="float")

varlist[[5]]  <- ncvar_def(name="BETAGPP_INFL_VAR",units="",
                                 dim=list(x,y,t), missval=NA,prec="float")

varlist[[6]]  <- ncvar_def(name="BETARESP_INFL_MEAN",units="",
                                 dim=list(x,y,t), missval=NA,prec="float")

varlist[[7]]  <- ncvar_def(name="BETARESP_INFL_VAR",units="",
                                 dim=list(x,y,t), missval=NA,prec="float")

varlist[[8]]  <- ncvar_def(name="BETAOCEAN_INFL_MEAN",units="",
                                 dim=list(x,y,t), missval=NA,prec="float")

varlist[[9]]  <- ncvar_def(name="BETAOCEAN_INFL_VAR",units="",
                                 dim=list(x,y,t), missval=NA,prec="float")

#fileout = "/user1/aschuh/temp/betas.112712.nc"

ncnew <- nc_create(fileout,varlist)

ncvar_put(ncnew, varlist[[1]],t(ens$ensembs.land.gpp))
ncvar_put(ncnew, varlist[[2]],t(ens$ensembs.land.resp))
ncvar_put(ncnew, varlist[[3]],t(ens$ensembs.ocean))
ncvar_put(ncnew, varlist[[4]],rep(infl.factor.mean,length(x$vals)*length(y$vals)))
ncvar_put(ncnew, varlist[[5]],rep(infl.factor.var,length(x$vals)*length(y$vals)))
ncvar_put(ncnew, varlist[[6]],rep(infl.factor.mean,length(x$vals)*length(y$vals)))
ncvar_put(ncnew, varlist[[7]],rep(infl.factor.var,length(x$vals)*length(y$vals)))
ncvar_put(ncnew, varlist[[8]],rep(infl.factor.mean,length(x$vals)*length(y$vals)))
ncvar_put(ncnew, varlist[[9]],rep(infl.factor.var,length(x$vals)*length(y$vals)))

nc_close(ncnew) 


