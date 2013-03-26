#--  This function creates a list of observations
#--  from available data (NOAA CT right now).  Right now
#--  we will assume that SITE by TIME will uniquely identify them
#--  to be matched against the ensemble output.  We may need
#--  x,y,z later.

#-- Test args
#-- noaa_data_dir = "/projects/School/geoschem/data/carbontracker.obs"
#-- min.time=2009;max.time=2009+60/365;
create_noaa_data = function(noaa_data_dir=NULL,min.time=2009,max.time=2009+60/365)
 {
  #-- check for NULLS
  if(is.null(noaa_data_dir))
  {
    stop("Need NOAA DATA DIR specified")	
  }
  
  require(ncdf4)
  
  #--  Pull all the NOAA carbontracker preprocessed data
                                                  
  fls = list.files(noaa_data_dir,full.names=TRUE)

  fls.short = list.files(noaa_data_dir,full.names=FALSE)

  #-- Basing on CT useage right now

  ind_no_use = fls.short %in% c("BAO_01C3_02LST_obs.nc",
                 "EFS_03C0_02LST_obs.nc",
                 "FEF_03C0_02LST_obs.nc",
                 "LEF_01C3_02LST_obs.nc",
                 "NWR_03C0_14LST_obs.nc",
                 "SCT_01C3_02LST_obs.nc",
                 "SNP_01C3_02LST_obs.nc",
                 "WBI_01C3_02LST_obs.nc",
                 "WKT_01C3_02LST_obs.nc",
                 "HDP_03C0_02LST_obs.nc",
                 "HDP_03C0_14LST_obs.nc",
                 "SPL_03C0_14LST_obs.nc")
 
  ind_no_use_plane = ((1:length(fls.short)) %in% grep("P2",fls.short))
 
  fls = fls[!ind_no_use & !ind_no_use_plane]
  fls.short = fls.short[!ind_no_use & !ind_no_use_plane]
    
  #-- Output matrix
  output = matrix(0,nrow=0,ncol=3)

  for(i in 1:length(fls))
  {
	print(paste("working on ",i," of ",length(fls)))
	fil = nc_open(fls[i])
	
	#-- Just need site name, time and meas
	site = substring(fls.short[i],1,8)
	#tim = ncvar_get(fil,"decimal_date")
	tim = ncvar_get(fil,"date")
        meas = ncvar_get(fil,"measured_value")
	
	ind = tim > min.time & tim < max.time
    
    datmat = matrix(NA,ncol=3,nrow=sum(ind))
    datmat[,1] = rep(site,sum(ind))
    datmat[,2] = tim[ind]
    datmat[,3] = meas[ind]
                
    output = rbind(output,datmat)	
	nc_close(fil)	
    }
   return(output)
  }
