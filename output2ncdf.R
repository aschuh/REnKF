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
   
output2ncdf = function(betas = X_post,fileout,output_biases) #OBSLANDBIAS=FALSE,OBSOCEANBIAS=FALSE)
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

#-- If we are going for land/ocean bias
#if(OBSOCEANBIAS){
#     bias.OCEAN = betas[15831,]
#     bias.OCEAN = c(mean(bias.OCEAN),bias.OCEAN)
#}

#if(OBSLANDBIAS){
#     bias.LAND = betas[15830,]
#     bias.LAND = c(mean(bias.LAND),bias.LAND)
#}

if(output_biases)
{
   bias.list = list()
   bias.list[[1]] = c(mean(betas[15830,-1]),betas[15830,-1])
   bias.list[[2]] = c(mean(betas[15831,-1]),betas[15831,-1])
   bias.list[[3]] = c(mean(betas[15832,-1]),betas[15832,-1])
   bias.list[[4]] = c(mean(betas[15833,-1]),betas[15833,-1])
   bias.list[[5]] = c(mean(betas[15834,-1]),betas[15834,-1])
   bias.list[[6]] = c(mean(betas[15835,-1]),betas[15835,-1])
   bias.list[[7]] = c(mean(betas[15836,-1]),betas[15836,-1])
   bias.list[[8]] = c(mean(betas[15837,-1]),betas[15837,-1])
}

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

#if(OBSLANDBIAS){
#land_offset_ncfile_feed = array(bias.LAND,dim=c(ens.size,1))
#}

#if(OBSOCEANBIAS){
#ocean_offset_ncfile_feed = array(bias.OCEAN,dim=c(ens.size,1))
#}

if(output_biases)
{
        ncfile_feed = list()
	for(k in 1:8)
        {
           ncfile_feed[[k]] = array(bias.list[[k]],dim=c(ens.size,1))
        }
}

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

#if(OBSLANDBIAS)
#{
#  varlist[[4]] <- ncvar_def(name="OBSLANDBIAS",units="",
#                                 dim=list(nens,t), missval=NA,prec="float")     
#}

#if(OBSOCEANBIAS)
#{
#   varlist[[5]] <- ncvar_def(name="OBSOCEANBIAS",units="",
#                                 dim=list(nens,t), missval=NA,prec="float")     
#}

output_bias_names = c("OBSLANDHBIAS","OBSLANDMBIAS",
                               "OBSOCEANBIAS","OBSS32_COEF","OBSB1_OFFSETCOEF",
                               "OBSALBEDO_2_H_COEF","OBSALBEDO_2_M_COEF","OBSDP_CLD_COEF")
 
if(output_biases)
{
  for(k in 4:(3+length(bias.list)))
  {
    varlist[[k]] <- ncvar_def(name=output_bias_names[k-3],units="",dim=list(nens,t), missval=NA,prec="float")
  }
}
                                
###########fileout = "/Volumes/user1/aschuh/temp/betasSCOTT.nc"
#-- fileout = "~/Desktop/tobedeleted/tester.nc"
#-- Pull data

ncnew <- nc_create(fileout,varlist)

ncvar_put(nc=ncnew, varid=varlist[[1]], vals=gpp_bias_ncfile_feed, start=c(1,1,1,1), count= c(-1,-1,-1,1))
ncvar_put(ncnew, varlist[[2]],resp_bias_ncfile_feed)
ncvar_put(ncnew, varlist[[3]],ocean_bias_ncfile_feed)

#if(OBSLANDBIAS)
#{
#  ncvar_put(ncnew, varlist[[4]],land_offset_ncfile_feed)
#}

#if(OBSOCEANBIAS)
#{
#  ncvar_put(ncnew, varlist[[5]],ocean_offset_ncfile_feed)
#}


if(output_biases)
{
 for(k in 1:length(bias.list))
  {
    ncvar_put(ncnew, varlist[[k+3]],ncfile_feed[[k]])
  }
}

nc_close(ncnew) 

}


output2ncdf_old = function(betas = X_post,fileout,OBSLANDBIAS=FALSE,OBSOCEANBIAS=FALSE)
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

#-- If we are going for land/ocean bias
if(OBSOCEANBIAS){
     bias.OCEAN = betas[15830,]
     bias.OCEAN = c(mean(bias.OCEAN),bias.OCEAN)
}

if(OBSLANDBIAS){
     bias.LAND = betas[15831,]
     bias.LAND = c(mean(bias.LAND),bias.LAND)
}


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

if(OBSLANDBIAS){
land_offset_ncfile_feed = array(bias.LAND,dim=c(ens.size,1))
}

if(OBSOCEANBIAS){
ocean_offset_ncfile_feed = array(bias.OCEAN,dim=c(ens.size,1))
}

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

if(OBSLANDBIAS)
{
  varlist[[4]] <- ncvar_def(name="OBSLANDBIAS",units="",
                                 dim=list(nens,t), missval=NA,prec="float")
}

if(OBSOCEANBIAS)
{
   varlist[[5]] <- ncvar_def(name="OBSOCEANBIAS",units="",
                                 dim=list(nens,t), missval=NA,prec="float")
}


ncnew <- nc_create(fileout,varlist)

ncvar_put(nc=ncnew, varid=varlist[[1]], vals=gpp_bias_ncfile_feed, start=c(1,1,1,1), count= c(-1,-1,-1,1))
ncvar_put(ncnew, varlist[[2]],resp_bias_ncfile_feed)
ncvar_put(ncnew, varlist[[3]],ocean_bias_ncfile_feed)

if(OBSLANDBIAS)
{
  ncvar_put(ncnew, varlist[[4]],land_offset_ncfile_feed)
}

if(OBSOCEANBIAS)
{
  ncvar_put(ncnew, varlist[[5]],ocean_offset_ncfile_feed)
}

nc_close(ncnew)

}


write_new_priors_nc = function(BETAOCEAN,BETAGPP,BETARESP,OBSLANDHBIAS=NULL,OBSLANDMBIAS=NULL,
                               OBSOCEANBIAS=NULL,OBSS32COEF=NULL,OBSB1_OFFSETCOEF=NULL,
                               OBSALBEDO_2_H_COEF=NULL,OBSALBEDO_2_M_COEF=NULL,OBSDP_CLD_COEF=NULL,fileout,grid.x,grid.y)
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

if(!is.null(OBSLANDHBIAS)){

varlist[[4]]  <- ncvar_def(name="OBSLANDHBIAS",units="",
                                 dim=list(nens,t), missval=NA,prec="float")
}

if(!is.null(OBSLANDMBIAS)){

varlist[[5]]  <- ncvar_def(name="OBSLANDMBIAS",units="",
                                 dim=list(nens,t), missval=NA,prec="float")
}

if(!is.null(OBSOCEANBIAS)){                                 
varlist[[6]]  <- ncvar_def(name="OBSOCEANBIAS",units="",
                                 dim=list(nens,t), missval=NA,prec="float")                                                                                                  
}                                 

if(!is.null(OBSS32COEF)){
varlist[[7]]  <- ncvar_def(name="OBSS32_COEF",units="",
                                 dim=list(nens,t), missval=NA,prec="float")
}

if(!is.null(OBSB1_OFFSETCOEF)){
varlist[[8]]  <- ncvar_def(name="OBSB1_OFFSETCOEF",units="",
                                 dim=list(nens,t), missval=NA,prec="float")
}

if(!is.null(OBSALBEDO_2_H_COEF)){
varlist[[9]]  <- ncvar_def(name="OBSALBEDO_2_H_COEF",units="",
                                 dim=list(nens,t), missval=NA,prec="float")
}

if(!is.null(OBSALBEDO_2_M_COEF)){
varlist[[10]]  <- ncvar_def(name="OBSALBEDO_2_M_COEF",units="",
                                 dim=list(nens,t), missval=NA,prec="float")
}

if(!is.null(OBSDP_CLD_COEF)){
varlist[[11]]  <- ncvar_def(name="OBSDP_CLD_COEF",units="",
                                 dim=list(nens,t), missval=NA,prec="float")
}

###########fileout = "/Volumes/user1/aschuh/temp/betasSCOTT.nc"
#-- fileout = "~/Desktop/tobedeleted/tester.nc"
#-- Pull data

ncnew <- nc_create(fileout,varlist)

ncvar_put(nc=ncnew, varid=varlist[[1]], vals=BETAGPP, start=c(1,1,1,1), count= c(-1,-1,-1,1))
ncvar_put(ncnew, varlist[[2]],vals=BETARESP)
ncvar_put(ncnew, varlist[[3]],vals=BETAOCEAN)

if(!is.null(OBSLANDHBIAS))
{
	ncvar_put(ncnew, varlist[[4]],vals=OBSLANDHBIAS)
}

if(!is.null(OBSLANDMBIAS))
{
        ncvar_put(ncnew, varlist[[5]],vals=OBSLANDMBIAS)
}

if(!is.null(OBSOCEANBIAS))
{
	ncvar_put(ncnew, varlist[[6]],vals=OBSOCEANBIAS)
}

if(!is.null(OBSS32COEF))
{
        ncvar_put(ncnew, varlist[[7]],vals=OBSS32COEF)
}

if(!is.null(OBSB1_OFFSETCOEF))
{
        ncvar_put(ncnew, varlist[[8]],vals=OBSB1_OFFSETCOEF)
}

if(!is.null(OBSALBEDO_2_H_COEF))
{
        ncvar_put(ncnew, varlist[[9]],vals=OBSALBEDO_2_H_COEF)
}

if(!is.null(OBSALBEDO_2_M_COEF))
{
        ncvar_put(ncnew, varlist[[10]],vals=OBSALBEDO_2_M_COEF)
}

if(!is.null(OBSDP_CLD_COEF))
{
        ncvar_put(ncnew, varlist[[11]],vals=OBSDP_CLD_COEF)
}
nc_close(ncnew) 

}

write_new_priors_nc2 = function(BETAOCEAN,BETAGPP,BETARESP,OBSOCEANBIAS,OBSLANDBIAS,fileout,grid.x,grid.y)
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

if(!is.null(OBSOCEANBIAS)){

varlist[[4]]  <- ncvar_def(name="OBSOCEANBIAS",units="",
                                 dim=list(nens,t), missval=NA,prec="float")
}

if(!is.null(OBSLANDBIAS)){

varlist[[5]]  <- ncvar_def(name="OBSLANDBIAS",units="",
                                 dim=list(nens,t), missval=NA,prec="float")
}


ncnew <- nc_create(fileout,varlist)

ncvar_put(nc=ncnew, varid=varlist[[1]], vals=BETAGPP, start=c(1,1,1,1), count= c(-1,-1,-1,1))
ncvar_put(ncnew, varlist[[2]],vals=BETARESP)
ncvar_put(ncnew, varlist[[3]],vals=BETAOCEAN)

if(!is.null(OBSOCEANBIAS))
{
        ncvar_put(ncnew, varlist[[4]],vals=OBSOCEANBIAS)
}

if(!is.null(OBSLANDBIAS))
{
        ncvar_put(ncnew, varlist[[5]],vals=OBSLANDBIAS)
}


nc_close(ncnew)

}
