#--- ENSEMBLE SCRIPT FOR GEOS-CHEM

set.seed(27)

#--Queing software, nasa, torque or sge
que_soft = "nasa"

#-- Propagation Choices are 'none', 'pure' or 'ct'
prop_model = "ct"

#-- Required libs

if(que_soft=="sge"){require(Rsge)}

#-- Required libs
require(ncdf4)
require(plyr)

#-- Use options
ensembles = 200
cycle_length = 14
cycles = 1:78
startcycle = cycles[1]
endcycle = cycles[length(cycles)]
inflation.factor = 1.15
#startdate = as.POSIXlt(strptime('2009-06-01 00:00:00', '%Y-%m-%d %H:%M:%S'),tz="GMT")
startdate = as.POSIXlt(strptime('2010-01-01 00:00:00', '%Y-%m-%d %H:%M:%S'),tz="GMT")
startdate$isdst = 0
estimate_land_ocean_bias = TRUE


#-- Command line args
args = commandArgs(TRUE)

print(args)

 if(length(args)==0){
    print("No arguments supplied.")
 }else{
    for(i in 1:length(args)){
         eval(parse(text=args[[i]]))
    }
   print(paste("startcycle now:",startcycle))
   print(paste("endcycle now:",endcycle))
}

#-- User directories
run_dir = "/discover/nobackup/aschuh/run.ncdf4.test/"
outdir = "/discover/nobackup/aschuh/GEOS-CHEM_output/ncdf4.test"
#run_dir = "/discover/nobackup/aschuh/run.landoceanbias_andfluxbias_test2/"
#outdir = "/discover/nobackup/aschuh/GEOS-CHEM_output/obsbias_plus_fluxbias_TEST2"
input_geos_file = paste(run_dir,"/input.geos",sep="")
orig_betas_file = "/discover/nobackup/aschuh/data/betas/betas.112613_wbiases.nc"
#orig_betas_file = "/discover/nobackup/aschuh/data/betas/betas.061213.nc"

#-- Set working directory and sun grid eng. options
setwd(paste(run_dir,"../run/ENSCODE",sep=""))

#-- Necessary code
source("merge_ens_gosat_ncdf_data.R")
#source("merge_ens_gosat_data.R")
source("optimize_betas_gosat.R")
source("output2ncdf.R")
source("utils.R")
source("create_noaa_data.R")
source("create_prior.R")
#source("jobscript.NASA.pbs.R")

#-- Checking outdir against a few lines in input.geos which MUST match
geos_inputfile_check(input_geos_file,outdir,prop_model)

if(que_soft=="sge"){sge.options(sge.use.cluster=TRUE,sge.save.global=TRUE,sge.remove.ﬁles=FALSE)}

#system(paste("cp ",outdir,"/betas/betas_cycle_prior_000.nc ",outdir,"/betas/betas_cycle_post_000.nc",sep=""))


#-- *Need to check that output directories are there
#--  We should check for existence of output files now*

for(i in startcycle:endcycle)
{
	#-- Adjust ensemble start date to cycle start date
	#-- Next time we run, CHANGE BACK THE RUNDATE EQUATION CODE BELOW
    #-- NEED TO CLEAN UP TIME, RIGHT NOW ONLY USING DAY BUT HOURS ARE SORT 
    #-- OF SCREWED UP IN STARTDATE AND RUNDATE
    
	#rundate = as.POSIXlt(startdate + 3600*24*(cycles[i]-1)*cycle_length)
	rundate = as.POSIXlt(startdate + 3600*24*((cycles[i]-1)*cycle_length))
	rdate_arg = paste( rundate$year + 1900, pad(rundate$mon+1,width=2,fill="0"),  
	                          pad(rundate$mday, width=2, fill="0"),sep="")
	   
        #prior_betas_tminus1_file = paste(outdir,"/betas/betas_cycle_prior_",pad(cycles[i]-1,width=3,fill="0"),".nc",sep="")
        #post_betas_tminus1_file =  paste(outdir,"/betas/betas_cycle_post_",pad(cycles[i]-1,width=3,fill="0"),".nc",sep="")
        #post_betas_t_file =  paste(outdir,"/betas/betas_cycle_post_",pad(cycles[i],width=3,fill="0"),".nc",sep="")

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
         #ifiles = paste(outdir,"/betas/betas_cycle_post_",vv,".nc",sep="")

         vv = sapply(c(0,max(0,cycles[i]-2):(cycles[i]-1)),pad,width=3,fill="0")
         ifiles = paste(outdir,"/betas/betas_cycle_post_",vv,".nc",sep="")
         pr_ind = vv == "000"
         ifiles[pr_ind] = orig_betas_file

         # This select last prior file
         ifiles = ifiles[length(pr_ind)]
         pr_ind = pr_ind[length(pr_ind)]

         if(estimate_land_ocean_bias)
         {
        	ret2 = create_prior_landoceanbias(ifiles=ifiles[length(ifiles)],pr_ind=pr_ind[length(ifiles)],ensembles=ensembles) 
        	#filf = nc_open(ifiles[length(ifiles)])
                #ret2 = list(OBSLANDBIAS=ncvar_get(filf,"OBSLANDBIAS"),OBSOCEANBIAS=ncvar_get(filf,"OBSOCEANBIAS"))                
                #if(dim(ret2$OBSLANDBIAS)>200){ret2$OBSLANDBIAS=ret2$OBSLANDBIAS[sample(1:1000,200)];ret2$OBSOCEANBIAS=ret2$OBSOCEANBIAS[sample(1:1000,200)]}
                #nc_close(filf)
                ret = create_prior(ifiles=ifiles,pr_ind=pr_ind,ensembles=ensembles)
        	write_new_priors_nc(BETAOCEAN=ret$BETAOCEAN,BETAGPP=ret$BETAGPP,
        	                    BETARESP=ret$BETARESP,
        	                    OBSLANDHBIAS=ret2$OBSLANDHBIAS,
        	                    OBSLANDMBIAS = ret2$OBSLANDMBIAS,
                                    OBSOCEANBIAS=ret2$OBSOCEANBIAS,
                                    OBSS32COEF=ret2$OBSS32COEF,
                                    OBSB1_OFFSETCOEF=ret2$OBSB1_OFFSETCOEF,
                                    OBSALBEDO_2_H_COEF=ret2$OBSALBEDO_2_H_COEF,
                                    OBSALBEDO_2_M_COEF=ret2$OBSALBEDO_2_M_COEF,
                                    OBSDP_CLD_COEF=ret2$OBSDP_CLD_COEF,
        	                    fileout=prior_betas_tminus1_file,grid.x=2.5,grid.y=2)
                #write_new_priors_nc2(BETAOCEAN=ret$BETAOCEAN,BETAGPP=ret$BETAGPP,BETARESP=ret$BETARESP,
                #                     OBSOCEANBIAS=ret2$OBSOCEANBIAS,OBSLANDBIAS=ret2$OBSLANDBIAS,
                #                     fileout=prior_betas_tminus1_file,grid.x=2.5,grid.y=2)
          }else{       	
       	    ret = create_prior(ifiles=ifiles,pr_ind=pr_ind,ensembles=ensembles)
       	    write_new_priors_nc(BETAOCEAN=ret$BETAOCEAN,BETAGPP=ret$BETAGPP,
       	                        BETARESP=ret$BETARESP,
       	                        fileout=prior_betas_tminus1_file,grid.x=2.5,grid.y=2)
          }

        }else
        {
         system(paste("cp ",post_betas_tminus1_file," ",prior_betas_tminus1_file,sep=""))
        }


        print(paste("Working on cycle",i))
        
        print(paste("running stuff like: ./geos ",1," ",cycles[i]," ",
                               rdate_arg," ",cycle_length," 0",sep=""))


      #-- Ensemble  run
      if(prop_model %in% c('pure','ct') )
       {
        print("Using mean propagation...")

        ###########################################
        #-- SGE (Sun Grid Engine)
        ###########################################
        if(que_soft=="sge"){
                      sge.parSapply(X=1:ensembles,
                               FUN=function(x) {system(paste("./geos ",x," ",cycles[i]," ",
                                rdate_arg," ",cycle_length," 0",sep=""))},
                                njobs=100) }
        ###########################################
        #-- end SGE
        ###########################################

        ###########################################
        #-- NASA (PBS), experimental
        ###########################################
        if(que_soft=="nasa2"){

         #working_dir = "/discover/nobackup/aschuh/run"
         reg.folder = "/discover/nobackup/aschuh/reg_folders/my_job_dir38"

         if(file.exists(reg.folder)){
          system(paste("rm -rf ",reg.folder,sep=""))
         }

         system(paste("mkdir ",reg.folder,sep=""))

         system(paste("ln -s ",run_dir,"/geos ",reg.folder,"/geos",sep=""))

         system(paste("cp ",run_dir,"/input.geos ",reg.folder,"/input.geos",sep=""))

         clusterExport(cluster,list=c("i","reg.folder","cycles","rdate_arg",
                               "cycle_length"))

         parSapply(cluster,X=1:ensembles,FUN=function(x) {system(paste("cd ",reg.folder,";./geos ",x," ",cycles[i]," ",
                                rdate_arg," ",cycle_length," 0 > output",x,sep=""))})

        }

       ###########################################
       #-- end NASA (PBS)
       ###########################################

       ###########################################
       #-- NASA (pods.sh), default NASA
       ###########################################

       if(que_soft=="nasa"){

         #working_dir = "/discover/nobackup/aschuh/run"
         reg.folder = "/discover/nobackup/aschuh/reg_folders/my_job_dir38"

         if(file.exists(reg.folder)){
          system(paste("rm -rf ",reg.folder,sep=""))
         }

         system(paste("mkdir ",reg.folder,sep=""))

         system(paste("ln -s ",run_dir,"/geos ",reg.folder,"/geos",sep=""))

         system(paste("cp ",run_dir,"/input.geos ",reg.folder,"/input.geos",sep=""))

         con = file(description=paste(reg.folder,"/exec.script",sep=""),open="w")

         for(k in 1:ensembles)
          {
           writeLines(paste("./geos ",k," ",cycles[i]," ",
                               rdate_arg," ",cycle_length," 0 > ",reg.folder,"/outfile.",k,sep=""),con=con)
          }

         close(con)

         system("/discover/nobackup/aschuh/pods.sh /discover/nobackup/aschuh/reg_folders/my_job_dir38/exec.script 6")

         #stop("forced stop")
       }

       ###########################################
       #-- end NASA (pods.sh)
       ###########################################


       ##############################################
       #--  use BatchJobs package and "torque" (CIRA)
       ##############################################

       if(que_soft=="torque"){

       #default.resources = list(queue="batch1", walltime="96:00:00")
       library(BatchJobs)
       working_dir = "/discover/nobackup/aschuh/run"
       reg.folder = "/discover/nobackup/aschuh/reg_folders/my_job_dir38"

       if(file.exists(reg.folder)){
         system(paste("rm -rf ",reg.folder,sep=""))
       }
       reg <- makeRegistry(id="my_reg", seed=123, work.dir=working_dir,file.dir=reg.folder)

       run.geos = function(x,cycles,cycle_length,i,rdate_arg,rerun,run_dir)
       {
 	      system(paste("cd ",run_dir,sep=""))
 	      system(paste(run_dir,"geos ","ENS_NUMBER"," ",cycles[i]," ",
                                rdate_arg," ",cycle_length," ",rerun,sep=""))
       }

       xs <- 1:ensembles

       batchMap(reg, run.geos,xs,more.args=list(run_dir=run_dir,rdate_arg=rdate_arg,
                              i=i,cycle_length=cycle_length,cycles=cycles,rerun=0))

       ids <- getJobIds(reg)
       chunk = 4
       ids_list = lapply(1:(ceiling(max(xs/chunk))),
              FUN=function(x){return(seq((x-1)*chunk+1,x*chunk,1))})
       llen = length(ids_list)

       ids_list[llen] =   list(as.vector(unlist(ids_list[llen]))[as.vector(unlist(ids_list[llen])) <= max(xs) ] )

       submitJobs(reg,ids_list,job.delay=TRUE)

       success = waitForJobs(reg,sleep=60)

       if(success){

                   #log_files = list.files(working_dir,full.names=TRUE,pattern="master.log")

                   #file.remove(log_files)

                  }
               }
         }

       ###############################################
       ##  End BatchJobs/torque (CIRA)
       ###############################################

     if(prop_model == 'none')
     {
         print("Not using mean propagation...")
         #-- This option fixes cycles[i] = 1 so that we don't propagate correction factors
         #-- but start at same prior every time
                 sge.parSapply(cl,1:ensembles,
              function(x) {system(paste("./geos ",x," ",1," ",
                                rdate_arg," ",cycle_length," 0",sep=""))},
                                njobs=ensembles)
     }

        if(estimate_land_ocean_bias){
                                      return.landmask.arg = 1
                                    }else{
                                      return.landmask.arg = NULL
                                         }

        #-- Pull ensemble data from run         
	   print("merging data...")
	   fulldat = merge_ens_gosat_ncdf_data(ensemble.dir=paste(outdir,"/gosat/",sep=""),
                           return.landmask=estimate_land_ocean_bias)

           if(estimate_land_ocean_bias){landmask = fulldat$landmask}else{landmask=NULL}
	
	    #-- Optimize the betas
	    print("optimizing ....")

         #-- All the optimization is done here

          #-- Loop over all, chunk by chunk

          if(dim(fulldat$fullensdat)[1]>10000)
          {

             #-- Watch out for this if it is exactly length 5001
          	brks = c(seq(1,dim(fulldat$fullensdat)[1],by=5000),dim(fulldat$fullensdat)[1])
          	for(k in 1:(length(brks)-1))
          	{
          	   print(paste("working on obs: ",brks[k]," to ",brks[k+1],sep=""))
          	   if(k==1){betas_arg = prior_betas_tminus1_file}
                   err_vec = (2*(1.2*fulldat$err[brks[k]:brks[k+1]]+0.25))^2 
          	   ret  = optimize_betas(betas_file=betas_arg,
          	                          Rdiag_vector=err_vec,
                                         ens_matrix=fulldat$fullensdat[brks[k]:brks[k+1],],
                                         obs_vector=fulldat$obs[brks[k]:brks[k+1]],method=2,
                                         localize=FALSE,diags=TRUE,estimate_land_ocean_bias=estimate_land_ocean_bias,
                                          landmask=landmask[brks[k]:brks[k+1]])
                  betas_arg = ret$X_post
                 
                 if(k==1){
                   S0_diags = cbind(fulldat$coords[brks[k]:brks[k+1],],ret$S0)
                        }else{
                   S0_diags = rbind(S0_diags,cbind(fulldat$coords[brks[k]:brks[k+1],],ret$S0))
          	 }
                }
          }else
          {
                #err_vec = (2*(1.2*fulldat$err+0.25))^2
          	err_vec = (fulldat$err)^2
                ret = optimize_betas(betas_file=prior_betas_tminus1_file,Rdiag_vector=err_vec,
          	                        ens_matrix=fulldat$fullensdat,obs_vector=fulldat$obs,method=2,
          	                        localize=FALSE,diags=TRUE,estimate_land_ocean_bias=estimate_land_ocean_bias,
                                        landmask=landmask,test=i)
                S0_diags = cbind(fulldat$coords,ret$S0)
                S1_diags = cbind(fulldat$coords,t(ret$S1))
           }

        X_post = ret$X_post

        write.table(S0_diags,file=paste(outdir,"/diags/diags_S0_post_",i,sep=""),sep="\t",row.names=FALSE)
        write.table(S1_diags,file=paste(outdir,"/diags/diags_S1_post_",i,sep=""),sep="\t",row.names=FALSE)

        #-- Inflate variance of ensemble for "mean propagation" case
        if(prop_model %in% c('pure','ct') & estimate_land_ocean_bias==TRUE )
        {
          #-- We are inflating the constant obs biases by larger amounts than the other pieces
          bias_inds = (dim(X_post)[1]-7):(dim(X_post)[1])
          X_post[-bias_inds,] = (X_post[-bias_inds,] - X_post[-bias_inds,1])*inflation.factor + X_post[-bias_inds,1]

          #-- Here we trying to inflate the sd's back up to about 0.15 while maintaining correlation pattern
          #bias_inds = c(15830,15831)
          #inflation.factor_consts = 0.25/mean(sd(X_post[bias_inds[1],]),sd(X_post[bias_inds[2],]))

          inflation.factor_consts = c(0.15,0.15,0.15,3,0.1,1,1,0.02)/apply(X_post[bias_inds,],1,sd)
          X_post[bias_inds,] = (X_post[bias_inds,] - X_post[bias_inds,1])*inflation.factor_consts + X_post[bias_inds,1]
        }

        if(prop_model %in% c('pure','ct') & estimate_land_ocean_bias==FALSE )
        {
          X_post = (X_post - X_post[,1])*inflation.factor + X_post[,1]
        }

	#-- Output the betas to netcdf for next cycle
	print("outputting new betas ...")
        #--output2ncdf_old
	output2ncdf(betas=X_post,fileout=paste(outdir,"/betas/betas_cycle_post_",
	              pad(cycles[i],width=3,fill="0"),".nc",sep=""),
                      output_biases=estimate_land_ocean_bias)
        #             OBSLANDBIAS=estimate_land_ocean_bias,OBSOCEANBIAS=estimate_land_ocean_bias)
	              
	#-- Create diagnostics
	#diags(X_post)
	
	#-- Remove objects
	rm(fulldat)
	rm(X_post)	

	#-- Single rerun of mean to get restart CO2 for next cycle
	print("launching mean beta rerun to generate new CO2 field ....")

       ###########################################
       #-- NASA (pods.sh), default NASA
       ###########################################

       if(que_soft == "nasa")
       {

 #working_dir = "/discover/nobackup/aschuh/run"
         reg.folder = "/discover/nobackup/aschuh/reg_folders/my_job_dir38"

         if(file.exists(reg.folder)){
          system(paste("rm -rf ",reg.folder,sep=""))
         }

         system(paste("mkdir ",reg.folder,sep=""))

         system(paste("ln -s ",run_dir,"/geos ",reg.folder,"/geos",sep=""))

         system(paste("cp ",run_dir,"/input.geos ",reg.folder,"/input.geos",sep=""))

         con = file(description=paste(reg.folder,"/exec.script",sep=""),open="w")

         # for(k in ((j-1)*block.size+1):(min(j*block.size,ensembles) ) )
         for(k in 1:1)
          {
           writeLines(paste("./geos ",k," ",cycles[i]," ",
                               rdate_arg," ",cycle_length," 1 > ",reg.folder,"/outfile.",k,sep=""),con=con)
          }

         close(con)

         system("/discover/nobackup/aschuh/pods.sh /discover/nobackup/aschuh/reg_folders/my_job_dir38/exec.script 1")

        }

       ###########################################
       #-- End  NASA (pods.sh), default NASA
       ###########################################

       ##############################################
       #--  use BatchJobs package and "torque" (CIRA)
       ##############################################

       if(que_soft == "torque")
       {

       #default.resources = list(queue="batch1", walltime="96:00:00")
       library(BatchJobs)
       working_dir = "/home/aschuh/run.noaa"
       reg.folder = "/home/aschuh/reg_folders/my_job_dir38"

       if(file.exists(reg.folder)){
         system(paste("rm -rf ",reg.folder,sep=""))
       }
       reg <- makeRegistry(id="my_reg", seed=123, work.dir=working_dir,file.dir="/home/aschuh/reg_folders/my_job_dir38")

       xs <- 1
       #batchMap(reg, run.geos.reload, xs)
       batchMap(reg, run.geos, xs,more.args=list(rdate_arg=rdate_arg,
                                 i=i,cycle_length=cycle_length,cycles=cycles,rerun=1))
       ids <- getJobIds(reg)

       submitJobs(reg, resources=default.resources,ids)

       success = waitForJobs(reg,sleep=300)

       if(success){

                   log_files = list.files(working_dir,full.names=TRUE,pattern="master.log")

                   file.remove(log_files)

               }
        }

       ##############################################
       #--  End    BatchJobs package and "torque" (CIRA)
       ##############################################

       ##############################################
       #--  Start SGE (sun grid engine)
       ##############################################

       if(que_soft == "sge" )
        {
         sge.parSapply(1,
              function(x) {system(paste("./geos 1 ",cycles[i]," ",rdate_arg," ",cycle_length," 1",sep=""))},
              njobs=1)
        }

     ##############################################
     #--  End  SGE (sun grid engine)
     ##############################################

     ##############################################
     #-
     #-  Cleanup/removal of files
     #-
     ##############################################

     #-- Remove gosat files, but leave FINALRUN files, this is important
     #-- because the model assumes that current cycle's
     #-- data is in this directory, regardless of name.  It simply lists and sorts files.

     modelfiles = list.files(paste(outdir,"/gosat/",sep=""),full.names=TRUE)

     #-- We want to keep "FINALRUN" files as well as prior CO2 guesses for each cyle
     first_files  = list.files(paste(outdir,"/gosat/",sep=""),full.names=TRUE,pattern="ens.0001")

     noremove_files  = list.files(paste(outdir,"/gosat/",sep=""),full.names=TRUE,pattern="ens.0001.nc.FINALRUN")

     prior_files = list.files(paste(outdir,"/gosat/",sep=""),full.names=TRUE,pattern="ens.0001.PRIOR")

     move_files = first_files[ (!first_files %in% noremove_files) & (!first_files %in% prior_files)]

     if(length(move_files) > 0)
      {
                 for(k in 1:length(move_files))
             {
                  #newfilenm = paste(move_files[k],".PRIOR",sep="")
                  newfilenm = gsub("\\.nc",".PRIOR.nc",move_files[k])
                  file.copy(move_files[k],newfilenm)
                  file.remove(move_files[k])
              }
       }

     remove_files = modelfiles[!modelfiles %in% first_files]

     file.remove(remove_files)

     #############################################
     #-- Same for surface "STATIONS" files
     #############################################

     modelfiles = list.files(paste(outdir,"/stations/",sep=""),full.names=TRUE)

     #-- We want to keep "FINALRUN" files as well as prior CO2 guesses for each cyle
     first_files  = list.files(paste(outdir,"/stations/",sep=""),full.names=TRUE,pattern="ens.0001")

     noremove_files  = list.files(paste(outdir,"/stations/",sep=""),full.names=TRUE,pattern="FINALRUN.ens.0001")

     prior_files = list.files(paste(outdir,"/stations/",sep=""),full.names=TRUE,pattern="ens.0001.PRIOR")

     move_files = first_files[ (!first_files %in% noremove_files) & (!first_files %in% prior_files)]

     if(length(move_files) > 0)
      {
                 for(k in 1:length(move_files))
             {
                  newfilenm = gsub("//.nc",".PRIOR.nc",move_files[k])
                  file.copy(move_files[k],newfilenm)
                  file.remove(move_files[k])
              }
       }

     remove_files = modelfiles[!modelfiles %in% first_files]

     file.remove(remove_files)

}
