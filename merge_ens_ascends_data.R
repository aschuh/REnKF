#--  This function merges the ensemble realizations with
#--  the observation matrix from "create_noaa_data" function

#-- Test args
#-- ensemble.dir = "/scratch4/aschuh/GEOS-CHEM_output/scottbetas2/ascends/"
#-- observation.matrix=output

merge_ens_ascends_data = function(ensemble.dir=NULL,model.res_x=2.5,model.res_y=2)
{
   #require(akima)  #for interp
   
   fls.model           = sort(list.files(ensemble.dir,full.names=TRUE))
   fls.finalrun.model  = list.files(ensemble.dir,full.names=TRUE,pattern="FINALRUN")
   fls.priors.model  = list.files(ensemble.dir,full.names=TRUE,pattern="PRIOR")
   fls.exc = c(fls.finalrun.model,fls.priors.model)
   fls.model = fls.model[!(fls.model %in% fls.exc)]
  
   fls.short.model     = sort(list.files(ensemble.dir,full.names=FALSE))
   fls.short.finalrun.model  = list.files(ensemble.dir,full.names=FALSE,pattern="FINALRUN") 
   fls.short.priors.model  = list.files(ensemble.dir,full.names=FALSE,pattern="PRIOR")    
   fls.short.exc = c(fls.short.finalrun.model,fls.short.priors.model)
   fls.short.model = fls.short.model[!(fls.short.model %in% fls.short.exc)]


   #-- Drop in the PSEUDO test case
   #pseudo.fls = sort(list.files("/scratch4/aschuh/GEOS-CHEM_output/control/ascends/",full.names=TRUE))
   #pseudo.fls.short = sort(list.files("/scratch4/aschuh/GEOS-CHEM_output/control/ascends/",full.names=FALSE))
   #fls = c(fls.model,pseudo.fls[pseudo.fls.short %in% fls.short.model])
   #fls.short = c(fls.short.model,pseudo.fls.short[pseudo.fls.short %in% fls.short.model])

   fls = fls.model
   fls.short = fls.short.model

   ensnum = sapply(fls.short,
                    FUN=function(x){
                    	  xx = strsplit(x,"\\.")[[1]]
                    	  return(xx[length(xx)])
                    	  })
    
    ensnum = as.numeric(as.character(ensnum))
    
    #-- For pseudo case
    #ensnum[(length(fls.model)+1):(length(fls))] = max(ensnum) + 1                	 
    #-------------------
                    	  
	

	tim.offset = 6*24*366 + 18*24*365
	
	for(i in 1:max(ensnum))
	{
		print(paste("working on ens ",i," of ",max(ensnum),sep=""))
		tmpfls = fls[ensnum == i]
		for(j in 1:length(tmpfls))
		 {
		  tmpdat = read.table(tmpfls[j],header=TRUE)
		  if(j==1){
		  	         #nobs = dim(tmpdat)[1]
		  	         enstmpdat = tmpdat
		  	      
		  }else{
		  	    enstmpdat = rbind(enstmpdat,tmpdat)
		  	  }
		 }
		 if(i==1)
		  {
		   nobs = dim(enstmpdat)[1]
		   fullensdat = array(NA,dim=c(nobs,max(ensnum)))
		   fullensdat[,1] = enstmpdat$TRA_001  
		   
		   #-- Get observations
		   obs = enstmpdat$XCO2
		      
		   #-- needed to adjust errors
		   lat.gcell = floor(enstmpdat$LAT + (90+0.5*model.res_y)/2)+1
		   lon.gcell = floor(enstmpdat$LON + (90+0.5*model.res_x)/2)+1
		   count = c(1)
		   for(k in 2:(dim(enstmpdat)[1]))
		    {
		      if(lat.gcell[k] == lat.gcell[k-1] & lon.gcell[k] == lon.gcell[k-1])
		   	   {
		   	   	count[k] = count[k-1] + 1
		   	   	}
		   	  else
		   	  {
		   	  	count[k] = 1
		   	  }
		    }
		    indx = dim(enstmpdat)[1]
		    while(indx > 0)
		    {
		    	if(count[indx] > 1)
		    	{
		    		count[indx:(indx-count[indx]+1)] = count[indx]
		    		indx = indx - count[indx]
		    	}else{indx = indx - 1}
		    }
		   err.adj = sqrt(count) * enstmpdat$ERROR
                 }
		  
		 else{
		 	   fullensdat[,i] = enstmpdat$TRA_001
		 	  }
    }
		
		
	#-- Add pseudo situation
	
	#for(j in 1:length(pseudo.fls))
	#	 {
	#	  tmpdat = read.table(pseudo.fls[j],header=TRUE)
	#	  if(j==1){enstmpdat = tmpdat}else{enstmpdat = rbind(enstmpdat,tmpdat)}
	#	 }
	#	 fullensdat= cbind(fullensdat,enstmpdat$TRA_001)
   
   return.data = list(fullensdat=fullensdat,err = err.adj,obs=obs)
   return(return.data)
   
   }
