#--  This function merges the ensemble realizations with
#--  the observation matrix from "create_noaa_data" function

#-- Test args
#-- ensemble.dir = "~/Desktop/tobedeleted/cycle1/"
#-- observation.matrix=output
merge_ens_data = function(observation.matrix,ensemble.dir=NULL,evaluate.final=FALSE)
{

		#-- This is where I need to pull in extra info on tower locations, grid-wise,
		#-- adjust coastal sites to ocean and mountain sites to higher levels based
		#-- on David's adjustments
		
		#statdat = read.table("/projects/School/geoschem/data/NOAAsites_for_gchem.txt",header=TRUE)
		#statdat = read.delim("/discover/nobackup/aschuh/run/ENSCODE/NOAAsites_for_gchem.122012.txt",header=TRUE)
   statdat = read.delim("/discover/nobackup/aschuh/run/ENSCODE/NOAAsites_for_gchem_nonames_053013_werrors.txt",header=TRUE)		

   #-- Restrict to stations that have ROWS 1,2,4 and 5 available
		
   #ind_na = apply(statdat[,c(1,2,4,5)],1,FUN=function(x){any(is.na(x))})
   ind_na = apply(statdat[,c(1,2,4,5,6)],1,FUN=function(x){any(is.na(x))})
   statdat = statdat[!ind_na,]
				
   statdat = statdat[statdat[,1] %in% unique(observation.matrix[,1]),]
		
   #-- Now statdat has only stations where I have vertical level info and error info
   #-- from David's info
   
   require(akima)  #for interp
   
   #fls = sort(list.files(ensemble.dir,full.names=TRUE,pattern="ens"))
   #fls.short = sort(list.files(ensemble.dir,full.names=FALSE,pattern="ens"))

   fls.model           = sort(list.files(ensemble.dir,full.names=TRUE))
   fls.finalrun.model  = list.files(ensemble.dir,full.names=TRUE,pattern="FINALRUN")
   fls.priors.model  = list.files(ensemble.dir,full.names=TRUE,pattern="PRIOR")
   if(evaluate.final){fls.exc = fls.priors.model}else{fls.exc = c(fls.finalrun.model,fls.priors.model)}
   fls.model = fls.model[!(fls.model %in% fls.exc)]
  
   fls.short.model     = sort(list.files(ensemble.dir,full.names=FALSE))
   fls.short.finalrun.model  = list.files(ensemble.dir,full.names=FALSE,pattern="FINALRUN") 
   fls.short.priors.model  = list.files(ensemble.dir,full.names=FALSE,pattern="PRIOR")
   if(evaluate.final){fls.short.exc = fls.short.priors.model}else{
           fls.short.exc = c(fls.short.finalrun.model,fls.short.priors.model)}    
   fls.short.exc = c(fls.short.finalrun.model,fls.short.priors.model)
   fls.short.model = fls.short.model[!(fls.short.model %in% fls.short.exc)]

   fls = fls.model
   fls.short = fls.short.model

   #-- Drop in the PSEUDO test case
   #fls = c(fls,"/Users/andrewschuh/Desktop/tobedeleted/scottbetas/stations.20090101.ens.0001.nc")
   #fls.short = c(fls.short,"/Users/andrewschuh/Desktop/tobedeleted/scottbetas/stations.20090101.ens.0001.nc")

   ensnum = sapply(fls.short,
                    FUN=function(x){
                    	  xx = strsplit(x,"\\.")[[1]]
                    	  return(xx[length(xx)-1])
                    	  })
                    	  
	ensnum = as.numeric(as.character(ensnum))
	
	#-- THis is based on Becky's specifid origin
        #tim.offset = 6*24*366 + 18*24*365
        tim.offset = 3*24*366 + 12*24*365	

	stations.in.data = unique(observation.matrix[,1])
	
	#-- determine the "common" stations to data and GEOSCHEM output
	fil = nc_open(fls[1])
	#nc_stations = ncvar_get(fil,"station")
        nc_stations = fil$dim$station$vals	

	#-- Subset to stations in GEOSCHEM output and those w/ LEVEL/ERR info from David
	fil_statnames = nc_stations[nc_stations %in% statdat[,1]]
	
	#-- Subset these stations to those that are ALSO in CT NOAA data
	common.stations = fil_statnames[fil_statnames %in% stations.in.data]
	nc_close(fil)

        #-- need to convert dates in observation matrix, Obspack data in decimal format
        #-- Becky's stations output netcdf file in hours since 1985-01-01
        observation.matrix[,2] = julian(date_decimal(as.numeric(as.character(observation.matrix[,2]))),
                                        origin = as.Date("1985-01-01"))*24
	
	obs.adj.matrix = observation.matrix[observation.matrix[,1] %in% common.stations,]
	
	for(i in 1:length(fls))
	{
		print(paste("working on ...",fls[i]))
		fil = nc_open(fls[i])
		fil_statnames = ncvar_get(fil,"station")
		co2 = ncvar_get(fil,"CO2")
		
		#-- only pull stations from ncdf that are in the obs data
		#-- and check that the stations are in same order
		
		co2 = co2[fil_statnames %in% common.stations,]
		
		if(any(fil_statnames[fil_statnames %in% common.stations] != stations.in.data[stations.in.data %in% common.stations])==TRUE)
		{
			stop("station problem in merge_ens_data.R")
		}
		
		co2 = as.vector(t(co2))
		#tim = 2009 + (fil$dim$time$vals - tim.offset)/(365*24)

                # CT is in days since 1/1/2000 and GEOSCHEM output is hours since 1/1/1985
                #  This sets them both in hours and starts them at same time
                #tim = (fil$dim$time$vals - tim.offset) # CT is in days since 1/1/2000		

                tim = fil$dim$time$vals

                nc_close(fil)
		
		statvec = as.vector(t(matrix(rep(statdat[,1],length(tim)),
		                     nrow=length(statdat[,1]),byrow=F)))
		
		errvec = as.vector(t(matrix(rep(statdat[,6],length(tim)),
		                     nrow=length(statdat[,6]),byrow=F)))
		
		timvec = rep(tim,dim(statdat)[1])
		dat = cbind(statvec,timvec,co2,errvec)     
		for(j in 1:length(statdat[,1]))
		{
	            tempstation  = statdat[j,1]		    
		    tempobs      = matrix(obs.adj.matrix[obs.adj.matrix[,1]==tempstation,],ncol=3)
		    tempmodelout = dat[dat[,1]==tempstation,]
		    	if(i==1 & j==1){
		    	co2_out = approx(x=as.numeric(as.character(tempmodelout[,"timvec"])),
		                         y=as.numeric(as.character(tempmodelout[,"co2"])),
		                         xout = as.numeric(as.character(tempobs[,2])))$y
		        err_out = rep(unique(dat[dat[,1]==tempstation,4]),dim(tempobs)[1])
		                   }
		    else{
		    	co2_out = c(co2_out,approx(x=as.numeric(as.character(tempmodelout[,"timvec"])),
		                                   y=as.numeric(as.character(tempmodelout[,"co2"])),
		                                   xout = as.numeric(as.character(tempobs[,2])))$y)	
		        if(i==1){err_out = c(err_out,rep(unique(dat[dat[,1]==tempstation,4]),dim(tempobs)[1]))}                         	    
                        }
		    		}                
	}
	
	co2_out_mat = matrix(as.numeric(as.character(co2_out)),nrow=dim(obs.adj.matrix)[1],byrow=FALSE)
	
	#-- For now we have some NAs because we have some observations which occur in first hour of 14 day chunk
	#-- but it looks like we might only have model observations starting at the end of the first hour
	
	ind_na = is.na(co2_out_mat[,1])
	
	#finalout = cbind(obs.adj.matrix,co2_out_mat)
	
	#-- For PSEUDO CASE
	#dimnames(finalout) = list(NULL,c("STATION","TIME","NOAA",paste("ENS",sapply(1:200,FUN=pad,width=4,fill="0"),sep=""),"PSEUDO"))
	#dimnames(finalout) = list(NULL,c("STATION","TIME","NOAA",paste("ENS",sapply(1:200,FUN=pad,width=4,fill="0"),sep="")))
	
	return.data = list(fullensdat=co2_out_mat[!ind_na,],err = as.numeric(as.character(err_out[!ind_na])),
	                     obs=as.numeric(as.character(obs.adj.matrix[!ind_na,3])))
    return(return.data)
}
