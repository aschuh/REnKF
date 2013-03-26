##################################################
#-- This is outline of code needed to write the beta
#-- out file, which will be fed into GEOSCHEM
#-- at the next time step.
##################################################

#-- Land=true for land, otherwise for Ocean
remap2full = function(v,mask)
{
	#-- mask is fed in internally now
    #if(land){wind <- LAND.MASK == 0}
    #else{wind <- OCEAN.MASK == 0}    
    wind = mask == 0
    tmp = vector(length=144*91)
    tmp[!wind] = v
    tmpm = t(matrix(tmp,ncol=144,byrow=T))
    return(tmpm)
}
   
output2ncdf = function(betas = X_post,fileout)
{

#-- Need to review this stuff a bit, masks aren't great, land.mask is a bit bigger than ocean.mask, why?
OCEAN.MASK = !is.na(rotate270.matrix(read.fwf("/discover/nobackup/aschuh/run/ENSCODE/Regions_ocean.dat",widths=rep(1,144))))
LAND.MASK = !is.na(rotate270.matrix(read.fwf("/discover/nobackup/aschuh/run/ENSCODE/Regions_land.dat",widths=rep(1,144))))

#-- Not sure, the GEOS masks seem to have the land mask w/ two extra rows which
#-- appear to be at north pole, have to check this one out, probably bug in their routines

LAND.MASK = LAND.MASK[,3:93]


map.GPP = apply(betas[1:3557,],2,FUN=remap2full,mask=LAND.MASK)
map.RESP = apply(betas[3558:(3557*2),],2,FUN=remap2full,mask=LAND.MASK)
map.OCEAN = apply(betas[(3557*2+1):(15829),],2,FUN=remap2full,mask=OCEAN.MASK)

#-- Add mean control in as first vector
map.GPP = cbind(apply(map.GPP,1,mean),map.GPP)
map.RESP = cbind(apply(map.RESP,1,mean),map.RESP)
map.OCEAN = cbind(apply(map.OCEAN,1,mean),map.OCEAN)


####################################################################
#--  Variables you must set these dim/gridcellsize correspond
#-- to GEOSCHEM run to be used fileout is where bias ens
#-- is stored.  Betas is nens - 1 for mean.  So, ens.size is this
#-- plus one. We MAY WANT TO WRITE THIS TO SOME KIND OF LOG FILE
####################################################################
origdim.x =  144
origdim.y =  91
grid.x = 2.5
grid.y = 2
ens.size = dim(betas)[2] 
####################################################################

gpp_bias_ncfile_feed = array(map.GPP,dim=c(144,91,ens.size,1))
resp_bias_ncfile_feed = array(map.RESP,dim=c(144,91,ens.size,1))
ocean_bias_ncfile_feed = array(map.OCEAN,dim=c(144,91,ens.size,1))

require(ncdf4)

x  = ncdim_def( "lon", "degrees_east", as.double(seq(-180,180-grid.x,by=grid.x))) 
y  = ncdim_def( "lat", "degrees_north", as.double(c(-90+0.25*grid.y,seq((-90+0.25*grid.y)+0.5*grid.y,
                                                         90-0.75*grid.y,by=grid.y),90-0.25*grid.y)))
                                                        
nens =  ncdim_def( "ensemble", units="integer", as.double(seq(1:(ens.size))))
                                                         
#z  = ncdim_def("Levels","hPa",as.double(1:47))
t  = ncdim_def( "time", "hours since 1900-01-01", 1, unlim=TRUE) 

varlist = list()

varlist[[1]]  <- ncvar_def(name="BETAGPP",units="",
                                 dim=list(x,y,nens,t), missval=NA,prec="float")
                                 
varlist[[2]]  <- ncvar_def(name="BETARESP",units="",
                                 dim=list(x,y,nens,t), missval=NA,prec="float")
                                 
varlist[[3]]  <- ncvar_def(name="BETAOCEAN",units="",
                                 dim=list(x,y,nens,t), missval=NA,prec="float")
                                 
###########fileout = "/Volumes/user1/aschuh/temp/betasSCOTT.nc"
#-- fileout = "~/Desktop/tobedeleted/tester.nc"
#-- Pull data

ncnew <- nc_create(fileout,varlist)

ncvar_put(nc=ncnew, varid=varlist[[1]], vals=gpp_bias_ncfile_feed, start=c(1,1,1,1), count= c(-1,-1,-1,1))
ncvar_put(ncnew, varlist[[2]],resp_bias_ncfile_feed)
ncvar_put(ncnew, varlist[[3]],ocean_bias_ncfile_feed)

nc_close(ncnew) 

}



write_new_priors_nc = function(BETAOCEAN,BETAGPP,BETARESP,fileout,grid.x,grid.y)
{
require(ncdf4)

x  = ncdim_def( "lon", "degrees_east", as.double(seq(-180,180-grid.x,by=grid.x))) 
y  = ncdim_def( "lat", "degrees_north", as.double(c(-90+0.25*grid.y,seq((-90+0.25*grid.y)+0.5*grid.y,
                                                         90-0.75*grid.y,by=grid.y),90-0.25*grid.y)))
ens.size = dim(BETAOCEAN)[3]                                                        
nens =  ncdim_def( "ensemble", units="integer", as.double(seq(1:(ens.size))))
                                                         
#z  = ncdim_def("Levels","hPa",as.double(1:47))
t  = ncdim_def( "time", "hours since 1900-01-01", 1, unlim=TRUE) 

varlist = list()

varlist[[1]]  <- ncvar_def(name="BETAGPP",units="",
                                 dim=list(x,y,nens,t), missval=NA,prec="float")
                                 
varlist[[2]]  <- ncvar_def(name="BETARESP",units="",
                                 dim=list(x,y,nens,t), missval=NA,prec="float")
                                 
varlist[[3]]  <- ncvar_def(name="BETAOCEAN",units="",
                                 dim=list(x,y,nens,t), missval=NA,prec="float")
                                 
###########fileout = "/Volumes/user1/aschuh/temp/betasSCOTT.nc"
#-- fileout = "~/Desktop/tobedeleted/tester.nc"
#-- Pull data

ncnew <- nc_create(fileout,varlist)

ncvar_put(nc=ncnew, varid=varlist[[1]], vals=BETAGPP, start=c(1,1,1,1), count= c(-1,-1,-1,1))
ncvar_put(ncnew, varlist[[2]],vals=BETARESP)
ncvar_put(ncnew, varlist[[3]],vals=BETAOCEAN)

nc_close(ncnew) 

}
