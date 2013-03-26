 #--- ENSEMBLE SCRIPT FOR GEOS-CHEM

 #--Queing software, nasa, torque or sge
 que_soft = "nasa"

 #-- Propagation Choices are 'none', 'pure' or 'ct'
 prop_model = "ct"

 #-- Required libs

 if(que_soft=="sge"){require(Rsge)}

 #-- If NASA, initialize cluster for duration of script run
 #if(que_soft=="nasa")
 #{
 #   library(Rmpi)
 #   np <- mpi.universe.size()
 #   library(snow)
 #   cluster <- makeMPIcluster(np)
 #}

 require(ncdf4)

 #-- Use options
 ensembles = 200
 cycle_length = 14
 cycles = 1:26
 inflation.factor = 1.15
 startdate = as.POSIXlt(strptime('2005-01-01 00:00:00', '%Y-%m-%d %H:%M:%S'),tz="GMT")
 startdate$isdst = 0

 #-- User directories
 run_dir = "/discover/nobackup/aschuh/run/"
 outdir = "/discover/nobackup/aschuh/GEOS-CHEM_output/longrun"
 input_geos_file = paste(run_dir,"/input.geos",sep="")

 #-- Set working directory and sun grid eng. options
 setwd(paste(run_dir,"ENSCODE",sep=""))

 #-- Necessary code
 source("merge_ens_data.R")
 source("optimize_betas.R")
 source("output2ncdf.R")
 source("utils.R")
 source("create_noaa_data.R")
 source("jobscript.NASA.pbs.R")

 #-- Checking outdir against a few lines in input.geos which MUST match
 geos_inputfile_check(input_geos_file,outdir,prop_model)

 if(que_soft=="sge"){sge.options(sge.use.cluster=TRUE,sge.save.global=TRUE,sge.remove.ﬁles=FALSE)}

 orig_betas_file = "/user1/aschuh/temp/betas.071212.nc"

 #-- for bamboo
 #system(paste("/usr/local/nco4-gcc/bin/ncks -F -d ensemble,0,",ensembles-1,
 #            " ",orig_betas_file," -o ",outdir,"/betas/betas_cycle_prior_000.nc",sep=""))

 #-- For sequoia
 #system(paste("ncks -F -d ensemble,1,",ensembles,
 #            " ",orig_betas_file," -o ",outdir,"/betas/betas_cycle_prior_000.nc",sep=""))


 system(paste("cp ",outdir,"/betas/betas_cycle_prior_000.nc ",outdir,"/betas/betas_cycle_post_000.nc",sep=""))


 #-- *Need to check that output directories are there
 #--  We should check for existence of output files now*

 for(i in 9:length(cycles))
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


     #-- Pull ensemble data from run
         #noaa_dir = "/home/aschuh/carbontracker.obs"
         noaa_dir = "/discover/nobackup/aschuh/carbontracker.obs"
         print("merging data...")

         #-- Check for leap day

         noaa_start_date = as.numeric(julian((i-1)*14*3600*24 + startdate,origin="2000-01-01 00:00:00 UTC"))
         noaa_end_date = as.numeric(julian(i*14*3600*24 + startdate,origin="2000-01-01 00:00:00 UTC"))

         observation.matrix = create_noaa_data(noaa_data_dir=noaa_dir,
                                                min.time=noaa_start_date,
                                                max.time=noaa_end_date)

         fulldat = merge_ens_data(observation.matrix,ensemble.dir=paste(outdir,"/stations/",sep=""))

         #-- Currently optimize_betas takes ens CO2 data w/ 10^-6 multiplier, CHANGE IN FUTURE
         fulldat$fullensdat = fulldat$fullensdat * 10^-6

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
                                         obs_vector=fulldat$obs[brks[k]:brks[k+1]],method=2,
                                         localize=TRUE,diags=TRUE)
                betas_arg = X_post

          	}
          }else
          {
          	X_post = optimize_betas(betas_file=prior_betas_tminus1_file,Rdiag_vector=fulldat$err^2,
          	ens_matrix=fulldat$fullensdat,obs_vector=fulldat$obs,method=2,
          	localize=TRUE,diags=TRUE)
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

 	#-- Remove gosat files, but leave FINALRUN files, this is important
 	#-- because the model assumes that current cycle's
 	#-- data is in this directory, regardless of name.  It simply lists and sorts files.

 	modelfiles = list.files(paste(outdir,"/gosat/",sep=""),full.names=TRUE)

 	#-- We want to keep "FINALRUN" files as well as prior CO2 guesses for each cyle
     first_files  = list.files(paste(outdir,"/gosat/",sep=""),full.names=TRUE,pattern="ens.0001")

     noremove_files  = list.files(paste(outdir,"/gosat/",sep=""),full.names=TRUE,pattern="ens.0001.FINALRUN")

     prior_files = list.files(paste(outdir,"/gosat/",sep=""),full.names=TRUE,pattern="ens.0001.PRIOR.nc")

     move_files = first_files[ (!first_files %in% noremove_files) & (!first_files %in% prior_files)]

     if(length(move_files) > 0)
      {
                 for(k in 1:length(move_files))
             {
                  newfilenm = paste(move_files[k],".PRIOR",sep="")
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
                  newfilenm = gsub(".nc",".PRIOR.nc",move_files[k])
                  file.copy(move_files[k],newfilenm)
                  file.remove(move_files[k])
              }
       }

     remove_files = modelfiles[!modelfiles %in% first_files]

     file.remove(remove_files)


 }