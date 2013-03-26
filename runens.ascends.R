#--- ENSEMBLE SCRIPT FOR GEOS-CHEM

#-- Propagation Choices are 'none', 'pure' or 'ct'
prop_model = "ct"

#-- Required libs
require(Rsge)
require(ncdf4)

source("/user1/aschuh/temp/ENSCODE/merge_ens_ascends_data.R")
source("/user1/aschuh/temp/ENSCODE/optimize_betas.R")
source("/user1/aschuh/temp/ENSCODE/output2ncdf.R")
source("/user1/aschuh/temp/ENSCODE/utils.R")

#-- Use options
ensembles = 200
cycle_length = 14
cycles = 1:26
inflation.factor = 1.15
startdate = as.POSIXlt(strptime('2007-01-01 00:00:00', '%Y-%m-%d %H:%M:%S'),tz="GMT")
startdate$isdst = 0

#-- User directories
run_dir = "/nacp-scratch2/aschuh/GEOS-CHEM/run"
outdir = "/scratch4/aschuh/GEOS-CHEM_output/whitepaper"
input_geos_file = paste(run_dir,"/input.geos",sep="")

#-- Checking outdir against a few lines in input.geos which MUST match
geos_inputfile_check(input_geos_file,outdir,prop_model)

#-- Set working directory and sun grid eng. options
setwd(run_dir)
sge.options(sge.use.cluster=TRUE,sge.save.global=TRUE,sge.remove.ﬁles=TRUE)

orig_betas_file = "/user1/aschuh/temp/betas.071212.nc"

#-- This needs to be fixed for sequoia, only works for bamboo
#system(paste("/usr/local/nco4-gcc/bin/ncks -F -d ensemble,0,",ensembles-1,
#            " ",orig_betas_file," -o ",outdir,"/betas/betas_cycle_prior_000.nc",sep=""))
#
#system(paste("cp ",outdir,"/betas/betas_cycle_prior_000.nc ",outdir,"/betas/beta_cycle_post_000.nc",sep=""))


#-- *Need to check that output directories are there
#--  We should check for existence of output files now*

for(i in 1:length(cycles))
{
	#-- Adjust ensemble start date to cycle start date
	
        #rundate = as.POSIXlt(startdate + 3600*24*(cycles[i]-1)*cycle_length)
        rundate = as.POSIXlt(startdate + 3600*24*((cycles[i]-1)*cycle_length))
        rdate_arg = paste( rundate$year + 1900, pad(rundate$mon+1,width=2,fill="0"),  
	                          pad(rundate$mday, width=2, fill="0"),sep="")
	
        prior_betas_tminus1_file = paste(outdir,"/betas/betas_cycle_prior_",pad(cycles[i]-1,width=3,fill="0"),".nc",sep="")
        post_betas_tminus1_file =  paste(outdir,"/betas/betas_cycle_post_",pad(cycles[i]-1,width=3,fill="0"),".nc",sep="")
        post_betas_t_file =  paste(outdir,"/betas/betas_cycle_post_",pad(cycles[i],width=3,fill="0"),".nc",sep="")

        #-- HERE IS WHERE WE NEED CDO CALL TO CREATE betas_cycle_prior_00X.nc file (X=cycles[i]-1)
        #-- as function of (X-1) and (X-2) files, (betas_000 + betas_(X-1) + betas+(X-2) ) / 3
        #-- Try to create in same ../betas folder if possible.  Final output betas should still
        #-- write to 'regularly' named betas file, w/o 'prior'
        #-- call will look something like (haven't tested yet):

        if(prop_model == "ct")
        {

        #vv = sapply(c(0,max(0,cycles[i]-2):(cycles[i]-1)),pad,width=3,fill="0")
        #ifiles = paste(outdir,"/betas/betas_cycle_post_",vv,".nc",sep="",collapse=" ")
        #system(paste("ncea ",ifiles," ",prior_betas_tminus1_file,sep=""))
        
        vv = sapply(c(0,max(0,cycles[i]-2):(cycles[i]-1)),pad,width=3,fill="0")
        ifiles = paste(outdir,"/betas/betas_cycle_post_",vv,".nc",sep="")

        for(k in 1:length(ifiles))
        {
        	if(k==1)
        	{
        		ifiles.fil = nc_open(ifiles[k])
        	    BETAGPP = ncvar_get(ifiles.fil,"BETAGPP")/length(ifiles)	
        	    BETARESP = ncvar_get(ifiles.fil,"BETARESP")/length(ifiles)
        	    BETAOCEAN = ncvar_get(ifiles.fil,"BETAOCEAN")/length(ifiles)  
        	    nc_close(ifiles.fil)      		
        	}else
        	{
        		ifiles.fil = nc_open(ifiles[k])
        	    BETAGPP = BETAGPP + ncvar_get(ifiles.fil,"BETAGPP")/length(ifiles)	
        	    BETARESP = BETARESP + ncvar_get(ifiles.fil,"BETARESP")/length(ifiles)
        	    BETAOCEAN = BETAOCEAN + ncvar_get(ifiles.fil,"BETAOCEAN")/length(ifiles)  
        	    nc_close(ifiles.fil)         		
        	}
        }
        
        write_new_priors_nc(BETAOCEAN,BETAGPP,BETARESP,prior_betas_tminus1_file,grid.x=2.5,grid.y=2)

        }
        else{
        system(paste("cp ",post_betas_tminus1_file," ",prior_betas_tminus1_file,sep=""))

        }

        print(paste("Working on cycle",i))

        print(paste("running stuff like: ./geos ",1," ",cycles[i]," ",
                               rdate_arg," ",cycle_length," 0",sep=""))


        #-- Ens run (sge sends to Grid Engine)
        if(prop_model %in% c('pure','ct') )
        {
                print("Using mean propagation...")
                sge.parSapply(1:ensembles,
             function(x) {system(paste("./geos ",x," ",cycles[i]," ",
                               rdate_arg," ",cycle_length," 0",sep=""))},
                               njobs=ensembles)
        }
    if(prop_model == 'none')
    {
        print("Not using mean propagation...")
        #-- This option fixes cycles[i] = 1 so that we don't propagate correction factors
        #-- but start at same prior every time
                sge.parSapply(1:ensembles,
             function(x) {system(paste("./geos ",x," ",1," ",
                               rdate_arg," ",cycle_length," 0",sep=""))},
                               njobs=ensembles)
    }
	

    #-- Pull ensemble data from run
        print("merging data...")
        fulldat = merge_ens_ascends_data(ensemble.dir=paste(outdir,"/ascends/",sep=""),model.res_x=2.5,model.res_y=2)

         #-- Optimize the betas
         print("optimizing ....")

         #-- All the optimization is done here

         #-- Loop over all, chunk by chunk

         if(dim(fulldat$fullensdat)[1]>5000)
         {

            #-- Watch out for this if it is exactly length 5001
         	brks = c(seq(1,dim(fulldat$fullensdat)[1],by=5000),dim(fulldat$fullensdat)[1])
         	for(k in 1:(length(brks)-1))
         	{
         	   print(paste("working on obs: ",brks[k]," to ",brks[k+1],sep=""))
         	   betas_arg = prior_betas_tminus1_file
         	   X_post = optimize_betas(betas_file=betas_arg,
         	                          Rdiag_vector=fulldat$err[brks[k]:brks[k+1]]^2,
                                        ens_matrix=fulldat$fullensdat[brks[k]:brks[k+1],],
                                        obs_vector=fulldat$obs[brks[k]:brks[k+1]],method=2,diags=TRUE)
               betas_arg = X_post
               
         	}
         }
         else
         {
         	X_post = optimize_betas(betas_file=prior_betas_tminus1_file,Rdiag_vector=fulldat$err^2,
                         ens_matrix=fulldat$fullensdat,obs_vector=fulldat$obs,method=2,diags=TRUE)
         }
         

        #-- Inflate variance of ensemble for "mean propagation" case
        if(prop_model %in% c('pure','ct') )
        {
          X_post = (X_post - X_post[,1])*inflation.factor + X_post[,1]
        }

        #-- Output the betas to netcdf for next cycle
        print("outputting new betas ...")
        output2ncdf(betas=X_post,fileout=paste(outdir,"/betas/betas_cycle_post_",
                      pad(cycles[i],width=3,fill="0"),".nc",sep=""))

        #-- Create diagnostics
        #diags(X_post)

        #-- Remove objects
        rm(fulldat)
        rm(X_post)

      #-- Single rerun of mean to get restart CO2 for next cycle
	print("launching mean beta rerun to generate new CO2 field ....")
	sge.parSapply(1,
             function(x) {system(paste("./geos 1 ",cycles[i]," ",rdate_arg," ",cycle_length," 1",sep=""))},
             njobs=1)

	#-- Remove ascends files, but leave FINALRUN files, this is important 
	#-- because the model assumes that current cycle's 
	#-- data is in this directory, regardless of name.  It simply lists and sorts files.
      
	modelfiles = list.files(paste(outdir,"/ascends/",sep=""),full.names=TRUE)     
	
	#-- We want to keep "FINALRUN" files as well as prior CO2 guesses for each cyle
        first_files  = list.files(paste(outdir,"/ascends/",sep=""),full.names=TRUE,pattern="ens.0001")

        noremove_files  = list.files(paste(outdir,"/ascends/",sep=""),full.names=TRUE,pattern="ens.0001.FINALRUN")

        prior_files = list.files(paste(outdir,"/ascends/",sep=""),full.names=TRUE,pattern="ens.0001.PRIOR")

        move_files = first_files[ (!first_files %in% noremove_files) & (!first_files %in% prior_files)]

        if(length(move_files) > 0)
        {
                for(i in 1:length(move_files))
            {
                 newfilenm = paste(move_files[i],".PRIOR",sep="")
                 file.copy(move_files[i],newfilenm)
                 file.remove(move_files[i])
             }
        }

    remove_files = modelfiles[!modelfiles %in% first_files]

    file.remove(remove_files)
}	
