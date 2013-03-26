#############################################################################
#--   This converts Becky's current bias files, based on Scott prescribed 
#--   biases to the "ensemble" format above, with ens_size=1.  This is the 
#--   the case that provides the pseudo-CO2.
#############################################################################

require(ncdf4)

x  = ncdim_def( "lon", "degrees_east", as.double(seq(-180,180-grid.x,by=grid.x))) 
y  = ncdim_def( "lat", "degrees_north", as.double(c(-90+0.25*grid.y,seq((-90+0.25*grid.y)+0.5*grid.y,
                                                         90-0.75*grid.y,by=grid.y),90-0.25*grid.y)))

ens_size = 1                                                         

nens =  ncdim_def( "ensemble", units="integer", as.double(seq(1:ens_size)))
                                                         
t  = ncdim_def( "time", "hours since 1900-01-01", 1, unlim=TRUE) 
	
varlist = list()

varlist[[1]]  <- ncvar_def(name="BETAGPP",units="",
                                 dim=list(x,y,nens,t), missval=NA,prec="float")
                                 
varlist[[2]]  <- ncvar_def(name="BETARESP",units="",
                                 dim=list(x,y,nens,t), missval=NA,prec="float")
                                 
varlist[[3]]  <- ncvar_def(name="BETAOCEAN",units="",
                                 dim=list(x,y,nens,t), missval=NA,prec="float")
                                 
fileout = "/Volumes/user1/aschuh/temp/betasSCOTT.nc"

#-- Pull data

fil.gpp = nc_open("/Volumes/scratch3/bmckeown/data/betas/betaGPP_2x25.nc")
fil.resp = nc_open("/Volumes/scratch3/bmckeown/data/betas/betaResp_2x25.nc")
fil.ocean = nc_open("/Volumes/scratch3/bmckeown/data/betas/betaOcean_2x25.nc")

gpp.bias = ncvar_get(fil.gpp,"betaGPP")
resp.bias = ncvar_get(fil.resp,"betaResp")
ocean.bias = ncvar_get(fil.ocean,"betaOcean")

nc_close(fil.gpp);nc_close(fil.resp);nc_close(fil.ocean)

ncnew <- nc_create(fileout,varlist)

ncvar_put(ncnew, varlist[[1]],as.vector(gpp.bias))
ncvar_put(ncnew, varlist[[2]],as.vector(resp.bias))
ncvar_put(ncnew, varlist[[3]],as.vector(ocean.bias))

nc_close(ncnew) 