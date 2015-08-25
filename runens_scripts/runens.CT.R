#--- ENSEMBLE SCRIPT FOR GEOS-CHEM

enkf.time.stamp <- "Time-stamp: <fe8.zeus.fairmont.rdhpcs.noaa.gov:/home/Andrew.Schuh/TM5/proj/ctemis/branches/ensemble/bin/enkf: 11 Mar 2015 (Wed) 13:20:21 UTC>"

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
ensembles = 150
cycle_length = 7
lag_window_length = 35
cycles = 1:(14*52)
startcycle = cycles[1]
endcycle = cycles[length(cycles)]
inflation.factor = 1.05
#startdate = as.POSIXlt(strptime('2009-06-01 00:00:00', '%Y-%m-%d %H:%M:%S'),tz="GMT")
startdate = as.POSIXct(strptime('2000-01-01 00:00:00', '%Y-%m-%d %H:%M:%S'),tz="GMT")
#startdate$isdst = 0
estimate_land_ocean_bias = FALSE
pods_numb = 28

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
#run_dir = "/discover/nobackup/aschuh/run.ascends.1RRV_no_obs_biases.fluxbias/"
#outdir = "/discover/nobackup/aschuh/GEOS-CHEM_output/ascends.1RRV_no_obs_biases.fluxbias"
run_dir = "/discover/nobackup/aschuh/run.v9-02_geosfp_CT_4x5/"
outdir = "/discover/nobackup/aschuh/GEOS-CHEM_output/v9.2_geosfp_CT_4x5"
input_geos_file = paste(run_dir,"/input.geos",sep="")
#orig_betas_file = "/discover/nobackup/aschuh/data/betas/betas.061213.nc"

#-- Set working directory and sun grid eng. options
setwd(paste(run_dir,"../run/ENSCODE",sep=""))

#-- Necessary code
source("load.ncdf4.R")
source("CT.covariance.R")
source("merge_ens_gosat_ncdf_data.R")
source("optimize_betas_test.R")
source("output2ncdf.R")
source("utils.R")
source("create_noaa_data.R")
source("create_prior.R")
#source("jobscript.NASA.pbs.R")

#-- Checking outdir against a few lines in input.geos which MUST match
geos_inputfile_check(input_geos_file,outdir,prop_model)

 #-- THIS PART MAKES NEW PBS_NODEFILE W/O HEAD NODE so that GEOS doesn't run on "this" node running R

    Sys.getenv('HOST')

    Sys.getenv()

    #system(paste("sed -i '/",Sys.getenv('HOSTNAME'),"/d' /discover/nobackup/aschuh/reg_folders/node_file",sep=""))
    #system(paste("sed -i '/",Sys.getenv('HOSTNAME'),"/d' ",Sys.getenv('PBS_NODEFILE'),sep=""))

    print(paste("Head node process running on ... ",Sys.getenv('HOSTNAME')))

    print(paste("running pods.sh on :"))

    #system("cat /discover/nobackup/aschuh/reg_folders/node_file")
    system(paste("cat ",Sys.getenv('PBS_NODEFILE')))

#-- *Need to check that output directories are there
#--  We should check for existence of output files now*

   ocn = "c"
   make.covariance <- make.covariance.oif
   if(ocn == "c") {
     make.covariance <- make.covariance.taka
    }

for(i in startcycle:endcycle)
{
	#-- Adjust ensemble start date to cycle start date
	#-- Next time we run, CHANGE BACK THE RUNDATE EQUATION CODE BELOW
    #-- NEED TO CLEAN UP TIME, RIGHT NOW ONLY USING DAY BUT HOURS ARE SORT
    #-- OF SCREWED UP IN STARTDATE AND RUNDATE

     	#rundate = as.POSIXlt(startdate + 3600*24*(cycles[i]-1)*cycle_length)
     	rundate = as.POSIXlt(startdate + 3600*24*((cycles[i]-1)*cycle_length))
        rundate2 = as.POSIXlt(startdate + 3600*24*((cycles[i])*cycle_length))

    	rdate_arg = paste( rundate$year + 1900, pad(rundate$mon+1,width=2,fill="0"),
	                          pad(rundate$mday, width=2, fill="0"),sep="")

        rdate_arg2 = paste( rundate2$year + 1900, pad(rundate2$mon+1,width=2,fill="0"),
                                  pad(rundate2$mday, width=2, fill="0"),sep="")

        rdate_arg_scaling_current = paste( rundate$year + 1900, pad(rundate$mon+1,width=2,fill="0"),
                                  pad(rundate$mday, width=2, fill="0"),sep="-")

        rdate_arg_scaling_next = paste( rundate2$year + 1900, pad(rundate2$mon+1,width=2,fill="0"),
                                  pad(rundate2$mday, width=2, fill="0"),sep="-")

        #prior_betas_tminus1_file = paste(outdir,"/betas/betas_cycle_prior_",pad(cycles[i]-1,width=3,fill="0"),".nc",sep="")
        #post_betas_tminus1_file =  paste(outdir,"/betas/betas_cycle_post_",pad(cycles[i]-1,width=3,fill="0"),".nc",sep="")
        #post_betas_t_file =  paste(outdir,"/betas/betas_cycle_post_",pad(cycles[i],width=3,fill="0"),".nc",sep="")

        prior_betas_tminus1_file = paste("/discover/nobackup/aschuh/GEOS-CHEM_output/v9.2_geosfp_CT_2x25/betas/scaling.factor.240.",rdate_arg_scaling_current,".nc",sep="")
        post_betas_tminus1_file = paste("/discover/nobackup/aschuh/GEOS-CHEM_output/v9.2_geosfp_CT_2x25/betas/scaling.factor.240.",rdate_arg_scaling_current,".nc",sep="")
        post_betas_t_file = paste("/discover/nobackup/aschuh/GEOS-CHEM_output/v9.2_geosfp_CT_2x25/betas/scaling.factor.240.",rdate_arg_scaling_next,".nc",sep="")
        
        print(paste("Working on cycle",i))

        print(paste("running stuff like: ./geos ",1," ",cycles[i],
                               rdate_arg," ",lag_window_length," ",cycle_length," 0",sep=""))


       ###########################################
       #-- NASA (pods.sh), default NASA
       ###########################################

       if(que_soft=="nasa"){

         #working_dir = "/discover/nobackup/aschuh/run"
         reg.folder = "/discover/nobackup/aschuh/reg_folders/my_job_dir40"

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
                               rdate_arg," ",lag_window_length," ",cycle_length," 0 > ",reg.folder,"/outfile.",k,sep=""),con=con)
          }

         close(con)

         #system(paste(" /usr/local/other/PoDS/PoDS/pods.py -x /discover/nobackup/aschuh/reg_folders/my_job_dir41/exec.script -n ",pods_numb))
         system(paste(" /discover/nobackup/aschuh/pods.sh /discover/nobackup/aschuh/reg_folders/my_job_dir40/exec.script ",pods_numb))

         #stop("forced stop")
       }

       ###########################################
       #-- end NASA (pods.sh)
       ###########################################

        #-- Pull ensemble data from run
	   print("merging data...")

#--  This function merges the ensemble realizations
#--  The resolution is needed to adjust model errors based upon
#--  multiple observations per grid cell

#-- Test args
#-- ensemble.dir = "/scratch4/aschuh/GEOS-CHEM_output/scottbetas2/ascends/"
#-- observation.matrix=output
#-- return.landmask = 1 for GEOSCHEM static landmask
#-- return.landmask = 2 return landmask based on ACOS surf/type

rejection.threshold = 3

t0 = rundate
t1 = rundate+60*60*24*7*5

t.in = t1 + seq(-6*24*60*60,0,by=24*60*60)

#-- NEEDS ADJUSTMENT BACK ONE DAY
#t.ins = sapply(t.in,FUN=format,format="%Y%m%d%H")
t.ins = sapply(t.in-24*60*60,FUN=format,format="%Y%m%d%H")
t.outs = sapply(t.in,FUN=format,format="%Y%m%d")

t0s <- format(t0,format="%Y%m%d%H")
t1s <- format(t1,format="%Y%m%d%H")

flasks.in <- sapply(t.ins,FUN=function(x){paste(outdir,"/../../run/prepped_daily_CT2013_20150204195501/flask_input.",x,".nc",sep="")})

flasks.out <- sapply(t.outs,FUN=function(x){paste(outdir,"/CT/CT.log.",x,".ens.0001",sep="")})

merge_ens_obspack_data = function(ensemble.dir=NULL)
{
   
   #require(akima)  #for interp
   require(ncdf4)
   require(abind)

   fls.model  = sort(sapply(t.outs,FUN=function(x){list.files(ensemble.dir,full.names=TRUE,pattern=x)}))
   fls.finalrun.model  = fls.model[grep("FINALRUN",fls.model)]
   fls.priors.model  = fls.model[grep("PRIOR",fls.model)]
   fls.exc = c(fls.finalrun.model,fls.priors.model)
   fls.model = fls.model[!(fls.model %in% fls.exc)]

   fls.short.model     = sort(sapply(t.outs,FUN=function(x){list.files(ensemble.dir,full.names=FALSE,pattern=x)}))
   fls.short.finalrun.model  = fls.short.model[grep("FINALRUN",fls.short.model)]
   fls.short.priors.model  = fls.short.model[grep("PRIOR",fls.short.model)]
   fls.short.exc = c(fls.short.finalrun.model,fls.short.priors.model)
   fls.short.model = fls.short.model[!(fls.short.model %in% fls.short.exc)]


   #-- Drop in the PSEUDO test case
   #pseudo.fls = sort(list.files("/scratch4/aschuh/GEOS-CHEM_output/control/ascends/",full.names=TRUE))
   #pseudo.fls.short = sort(list.files("/scratch4/aschuh/GEOS-CHEM_output/control/ascends/",full.names=FALSE))
   #fls = c(fls.model,pseudo.fls[pseudo.fls.short %in% fls.short.model])
   #fls.short = c(fls.short.model,pseudo.fls.short[pseudo.fls.short %in% fls.short.model])

   fls = fls.model
   fls.short = fls.short.model

   #dts.files = sort(unique(as.vector(unlist(sapply(fls.short,
   #                 FUN=function(x){
   #                       xx = strsplit(x,"\\.")[[1]][3]
   #                       return(xx[length(xx)])
   #                       })))))

   fil_ens_num = sapply(fls.short,FUN=function(x){substring(x,nchar(x)-7,nchar(x))})

      load.ncdfs.in = function(x){
	for(i in 1:length(x))
	{
		print(paste("reading ",x[i]," ....",sep=""))
		fin <- load.ncdf(x[i])
                 print(paste("reading..",length(fin$obspack.id),"obs"))	
	    if(i==1){fin.tot = fin
	    }else{
		fin.tot = list(time=c(fin.tot$time,fin$time),
		               latitude=c(fin.tot$latitude,fin$latitude),
		               longitude=c(fin.tot$longitude,fin$longitude),
		               altitude=c(fin.tot$altitude,fin$altitude),
		               value=c(fin.tot$value,fin$value),
		               time.decimal=c(fin.tot$time.decimal,fin$time.decimal),
		               time.components=cbind(fin.tot$time.components,fin$time.components),
		               obspack.id=c(fin.tot$obspack.id,fin$obspack.id),
		               obs.flag=c(fin.tot$obs.flag,fin$obs.flag),
		               sampling.strategy=c(fin.tot$sampling.strategy,fin$sampling.strategy),
		               mdm=c(fin.tot$mdm,fin$mdm),
		               may.reject=c(fin.tot$may.reject,fin$may.reject),
		               may.localize=c(fin.tot$may.localize,fin$may.localize),
		               obs=c(fin.tot$obs,fin$obs),
		               calendar.components=1:6)
	}
	}
	return(fin.tot)
}


     load.ncdfs.out = function(x){
        unq = unique(fil_ens_num)
	for(i in 1:length(unq))
	{
            tmpx = x[grep(unq[i],x)]
            for(j in 1:length(tmpx))
              {    
		print(paste("reading ",tmpx[j]," ....",sep=""))
		fout <- load.ncdf(tmpx[j])
	        print(paste("reading..",length(fout$obspack.id),"obs"))
	    if(j==1){fout.tot = fout}else{
		fout.tot = list(obspack.id=c(fout.tot$obspack.id,fout$obspack.id),
		               #flask=rbind(fout.tot$flask,fout$flask),  #-- need this for all 6 tracers
                               flask=c(fout.tot$flask,fout$flask),
		               nsamples=c(fout.tot$nsamples,fout$nsamples),
		               tracer.names=fout$tracer.names,
		               averaging.time=c(fout.tot$averaging.time,fout$averaging.time),
		               surface.height=c(fout.tot$surface.height,fout$surface.height),
		               region.name=fout$region.name,
		               region.indices=rbind(fout.tot$region.indices,fout$region.indices),
		               u=c(fout.tot$u,fout$u),
		               v=c(fout.tot$v,fout$v),
		               blh=c(fout.tot$blh,fout$blh),
		               q=c(fout.tot$q,fout$q),
		               pressure=c(fout.tot$pressure,fout$pressure),
		               temperature=c(fout.tot$temperature,fout$temperature),
		               obs=c(fout.tot$obs,fout$obs),
		               tracer=fout$tracer,
		               gridind=fout$gridind)
	                    }
             }
          if(i==1){fout.mat = fout.tot}else{
                   fout.mat = list(obspack.id=cbind(fout.mat$obspack.id,fout.tot$obspack.id),
                               #flask=abind(fout.mat$flask,array(fout.tot$flask,dim=c(dim(fout.tot$flask)[1],dim(fout.tot$flask)[2],1)),along=3),
                               flask=cbind(fout.mat$flask,fout.tot$flask),
                               nsamples=cbind(fout.mat$nsamples,fout.tot$nsamples),
                               tracer.names=fout$tracer.names,
                               averaging.time=cbind(fout.mat$averaging.time,fout.tot$averaging.time),
                               surface.height=cbind(fout.mat$surface.height,fout.tot$surface.height),
                               region.name=fout$region.name,
                               region.indices=abind(fout.mat$region.indices,fout.tot$region.indices,along=3),
                               u=cbind(fout.mat$u,fout.tot$u),
                               v=cbind(fout.mat$v,fout.tot$v),
                               blh=cbind(fout.mat$blh,fout.tot$blh),
                               q=cbind(fout.mat$q,fout.tot$q),
                               pressure=cbind(fout.mat$pressure,fout.tot$pressure),
                               temperature=cbind(fout.mat$temperature,fout.tot$temperature),
                               obs=cbind(fout.mat$obs,fout.tot$obs),
                               tracer=fout$tracer,
                               gridind=fout$gridind)

           }           
	}
	return(fout.mat)
}

  fin <- load.ncdfs.in(flasks.in)

 fout <- load.ncdfs.out(fls.model) 
   
   ensnum = sapply(fls.short,
                    FUN=function(x){
                    	  xx = strsplit(x,"\\.")[[1]]
                          #-- when .nc suffix is added to the CT out files, below should change: "length(xx)" -> "length(xx) -1"
                    	  return(xx[length(xx)])
                    	  })

    ensnum = as.numeric(as.character(ensnum))

 out=list(fin=fin,fout=fout,ensnum=ensnum)
 return(out)
 }


	   out = merge_ens_obspack_data(ensemble.dir=paste(outdir,"/CT/",sep=""))

           fin = out$fin
           fout = out$fout

nobs <- length(fout$obspack.id[,1])
nmemb <- dim(fout$flask)[3]
hqhr <- rep(FALSE,nobs)

nobs.fin <- length(fin$time)
is.dupe <- rep(FALSE,nobs.fin)
n.dupes <- 0

for (iobs in 1:nobs.fin) {

  # Skip this observation if it is already identified as a duplicate
  if(is.dupe[iobs]) {
    next
  }

  # skip this observation if it is not to be assimilated (has missing MDM)
  if(is.na(fin$mdm[iobs])) {
    next
  }
  
  dx2 <- (fin$longitude-fin$longitude[iobs])^2
  dy2 <- (fin$latitude-fin$latitude[iobs])^2
  
  lx <- which((abs(as.numeric(difftime(fin$time[iobs],fin$time,units='secs'))) < 3600) &
              (abs(fin$altitude - fin$altitude[iobs]) < 10) &
              (sqrt(dx2+dy2) < 0.05) &
              !is.na(fin$mdm))

  if(length(lx) < 2 ) {  # had better be at least 1.
    next
  }

  is.dupe[lx] <- TRUE
  fin$mdm[lx] <-   fin$mdm[lx]*sqrt(length(lx))
  
  n.dupes <- n.dupes + 1
}

self = NULL
cat(sprintf("%s %d instances of duplicate obs found.\n",self,n.dupes))

# match input and output observations using obspack.id
# note that there may be more input obs than output
#lx <- match(trim(fout$obspack.id),trim(fin$obspack.id))
lx <- match(fout$obspack.id,fin$obspack.id)

if(any(is.na(lx))) {

#-- Don't think we need this kill switch if R is managing the runs
#  # kill the run
#  tm5.ok <- sprintf("%s/tm5.ok",rundir)
#  if(file.exists(tm5.ok)) {
#    file.remove(tm5.ok)
#  }
#---
  cat(sprintf("\n%s ERROR!  There is a mismatch between flask_input and flask_output obspack IDs.\n",self))
  cat(sprintf("%s  flask input file: %s\n",self,flask.in))
  cat(sprintf("%s  flask output file: %s\n",self,flask.out))
  cat(sprintf("%s Run %s has been stopped.\n\n",self,runid))
  quit(save='no',status=1)
}

is.dupe <- is.dupe[lx]
obs <- 1e6*fin$value[lx]  #  observed values
r <- (fin$mdm[lx])^2 # assumed error on observations
may.localize <- as.logical(fin$may.localize[lx])
may.reject <- as.logical(fin$may.localize[lx])

reject <- rep(FALSE,nobs)

#--- OUTPUT NEEDS TO BE STACKED UP!!!!!!!!(back to fout)
#--  dropped 1e6 scalar for now as well
#dy <- 1e6*fout$flask[,1,]# deviations of sampled observations; dim nmemb, nobs (n.b. change from mol mol-1 to micromol mol-1)
#dy <- fout$flask[,1,] # this is if you're using multiple tracers, 3 dims here
dy <- fout$flask # deviations of sampled observations; dim nmemb, nobs (n.b. change from mol mol-1 to micromol mol-1)

#-- DIFFERENT
dy = t(dy)
#--

# remove ensemble member 1 prediction of co2 observations
# this is the prediction of the "mean" state.
sim <- dy[1,]

for (iobs in 1:nobs) {
  dy[,iobs] <- dy[,iobs] - dy[1,iobs]
}

# transpose dy to nobs, nmemb
dy <- t(dy)

parm.nc.in = paste(outdir,"/betas/scaling.factor.240.",format(t0,format="%Y-%m-%d"),".nc",sep="")
#--
jobstep.step = 7
#---
tnext <- t0 + 86400*jobstep.step
parm.nc.out  = paste(outdir,"/betas/scaling.factor.240.",format(tnext,format="%Y-%m-%d"),".nc",sep="")

# read in parameters
params.in <- load.ncdf(parm.nc.in)

nweeks <- dim(params.in$scaling.factor)[2]  # this makes nweeks = nlag + 1 
nmemb <- dim(params.in$scaling.factor.deviations)[2]
nparms <- dim(params.in$scaling.factor.deviations)[1]

nparms.tot <- nparms*(nweeks-1) # nweeks-1 is nlag

# Here we go from 2:nweeks to drop the first, advance-background week
#  dim(params.in$scaling.factor.deviations) =  nparm, nmemb, nweeks
#  dim(params.in$scaling.factor) = nparm, nweeks
dx <- params.in$scaling.factor.deviations[,,2]
x <-  params.in$scaling.factor[,2]
for (ilag in 3:nweeks) {
  dx <- rbind(dx,params.in$scaling.factor.deviations[,,ilag] )
  x <-  c(x,params.in$scaling.factor[,ilag] )
}

#pb <- progress.bar.start(sprintf("%s considering %d obs",self,nobs),nobs)

n.assim <- 0
n.unassim <- 0
n.rejected <- 0
assim.warnings <- character(0)

for (iobs in 1:nobs) {

  if(is.na(r[iobs])) {
#    cat(sprintf("\n%s WARNING!  Skipping observation %d, with obspack_id %s, due to NA model-data mismatch.\n",self,iobs,fin$obspack.id[iobs]))
    n.unassim <- n.unassim + 1
    next
  }
  
  if(any(is.na(dy[iobs,]))) {  # dy has been transposed to nobs, nmemb

    assim.warnings <- c(assim.warnings,sprintf("Observation %d, with obspack_id %s, has NA samples from model..rejected!\n",iobs,fin$obspack.id[iobs]))

#    cat(sprintf("\n%s WARNING!  Observation %d, with obspack_id %s, has NA samples from model..rejected!\n",self,iobs,fin$obspack.id[iobs]))
#    cat(sprintf("%s  flask output file: %s\n",self,flask.in))
    reject[iobs] <- TRUE
    n.rejected <- n.rejected + 1
    next

  }

  if(sd(dy[iobs,])==0) {
    assim.warnings <- c(assim.warnings,sprintf("Observation %d, with obspack_id %s, has UNIFORM samples from model..rejected!\n",iobs,fin$obspack.id[iobs]))
#    cat(sprintf("\n%s WARNING!  Observation %d, with obspack_id %s, has UNIFORM samples from model..rejected!\n",self,iobs,fin$obspack.id[iobs]))
#    cat(sprintf("%s  flask output file: %s\n",self,flask.in))
    reject[iobs] <- TRUE
    n.rejected <- n.rejected + 1
    next
  }
  
  # make K and rho
  K <- rep(NA,nparms.tot)
  rho <- rep(NA,nparms.tot) # correlation coefficient
  alpha <- NA # Whitaker and Hamill scaling coefficient

  for(iparm in 1:nparms.tot) {
    K[iparm] <- (1/(nmemb-1))*sum(dx[iparm,] * dy[iobs,])  # this is now QH^T
    rho[iparm] <- cor(dx[iparm,],dy[iobs,]) # might not be exactly what Wouter computes
  }
  hqhr[iobs] <- (1/(nmemb-1))*sum(dy[iobs,]*dy[iobs,])+r[iobs] # note r is already squared error

  # check for rejection. 
  chi <- (sim[iobs]-obs[iobs])/sqrt(hqhr[iobs])
  if((abs(chi)>rejection.threshold) & (may.reject[iobs])) {
    assim.warnings <- c(assim.warnings,sprintf("Rejecting observation %d with obspack_id %s: chi is %g.\n",iobs,fin$obspack.id[iobs],chi))
    #cat(sprintf("%s Rejecting observation %d with obspack_id %s: chi is %g.\n",self,iobs,fin$obspack.id[iobs],chi))
    reject[iobs] <- TRUE
    n.rejected <- n.rejected + 1
  }

  if(reject[iobs]) {
    next
  }
  
  alpha <- 1/(1+sqrt(r[iobs]/hqhr[iobs]))
  K <- K/hqhr[iobs]

  # localize
  # Note that 1.97591 = qt(p=0.975,df=150)
  tval <- 1.97591 # for 150 members only
  if(may.localize[iobs]) {
    for(iparm in 1:nparms.tot) {
      prob <- rho[iparm]/sqrt((1-rho[iparm]^2)/(nmemb-2))
      if(abs(prob) < tval) {
#        cat(".")
        K[iparm] <- 0
      } else {
#        cat("x")
      }
    }
  }

  # solve_enkf

  # invert central value
  resid <- obs[iobs] - sim[iobs]
  x <- x+K*resid

  # invert members
  for (imemb in 1:nmemb) {
    res <- dy[iobs,imemb]
    dx[,imemb] <- dx[,imemb] - alpha*K*res
  }

  for (jobs in 1:nobs) {
    if (iobs == jobs) {
      next
    }
    res <- obs[iobs] - sim[iobs]
    fac <- (1/(nmemb-1))*sum(dy[jobs,]*dy[iobs,])/hqhr[iobs]
    sim[jobs] <- sim[jobs]+fac*res
    dy[jobs,] <- dy[jobs,]-alpha*fac*dy[iobs,]
  }
  n.assim <- n.assim + 1
  #pb <- progress.bar.print(pb,iobs)
}
#progress.bar.end(pb)


scaling.factor <- matrix(NA,nrow=nparms,ncol=nweeks)
scaling.factor.deviations <- array(NA,dim=c(nparms,nmemb,nweeks))

cat(sprintf("%s - %d obs, of which %d were assimilated, %d rejected, and %d were not-for-assimilation.\n",
            self,nobs,n.assim,n.rejected,n.unassim))

if(length(assim.warnings) > 0) {
  for (iwarn in 1:(length(assim.warnings))) {
    cat(sprintf("%s [warning %d] %s\n",self,iwarn,assim.warnings[iwarn]))
  }
}
scaling.factor[,1:(nweeks-1)] <- x

for (jweek in 1:(nweeks-1)) {
  scaling.factor.deviations[,,jweek] <- dx[nparms*(jweek-1)+1:nparms,]
}

# add one week at the head of the filter
scaling.factor[,nweeks] <- (1/3)*(rep(1,nparms) + scaling.factor[,nweeks-2] + scaling.factor[,nweeks-1])

#--
bio.cov = "/home/aschuh/covariance_bio_olson19.nc"
ocn.cov = "/home/aschuh/Takahashi2009_qprior_30reg_v2.nc"
ocn.input.dir = ""
land.covariance.inflation.factor = 4
#--

# reproduce a minor bug from CT2010:  code computes covariance for lag+0 week,
# despite the fact that it is creating covariance for lag+4
#cat(sprintf("WARNING: computing covariance for wrong week, to mimic cy2 EnKF\n"))
#covariance.prior <- make.covariance(nparms,bio.cov,ocean.nc.dir,deltapco2.prefix,t0)
if (ocn == 'i') {
  covariance.prior <- make.covariance(nparms,bio.cov,ocean.nc.dir,deltapco2.prefix,t1)
} else {
  covariance.prior <- make.covariance(nparms,bio.cov,ocn.cov)
}
scaling.factor.deviations[,,nweeks] <- compute.prior.dx(covariance=covariance.prior,nmemb=nmemb,iweek=iweek)

# write out parameter results
epoch <- ISOdatetime(2000,1,1,0,0,0,"UTC")
tunits <- "days"
timevals <- as.numeric(difftime(tnext,epoch,units=tunits))

dim.nweeks <- ncdim_def("nweeks","weeks",vals=1:nweeks,unlim=FALSE)
dim.nmemb <- ncdim_def("nmembers","",vals=1:nmemb,unlim=FALSE)
dim.nparms <- ncdim_def("parameters","",vals=1:nparms,unlim=FALSE)

dim.time <- ncdim_def("date",
                      sprintf("%s since %s",tunits,format(epoch,format="%Y-%m-%d %H:%M:%S UTC")),vals=timevals,unlim=TRUE)
    
var.par <- ncvar_def(name="scaling_factor",units="NA",dim=list(dim.nparms,dim.nweeks,dim.time),
                     missval=-1e34,longname="EnKF_scaling_factor")

var.dev <- ncvar_def(name="scaling_factor_deviations",units="NA",dim=list(dim.nparms,dim.nmemb,dim.nweeks,dim.time),
                     missval=-1e34,longname="EnKF_scaling_factor_deviations")

var.cov <- ncvar_def(name="prior_covariance",units="NA",dim=list(dim.nparms,dim.nparms),
                     missval=-1e34,longname="Prior covariance for new week introduced at head of filter")

ncf <- nc_create(parm.nc.out,vars=list(var.par,var.dev,var.cov))
ncvar_put(ncf,var.par,scaling.factor)
ncvar_put(ncf,var.dev,scaling.factor.deviations)
ncvar_put(ncf,var.cov,covariance.prior)

ncatt_put(ncf,0,"Valid",
          attval=sprintf("Parameters valid for the %d-week period starting %s",nweeks,format(tnext,format="%Y-%m-%d %H:%M:%S UTC")),
          prec="text")

ncatt_put(ncf,0,"history",
          attval=sprintf("Created on %s\nby script '%s'",
            format(Sys.time(), "%a %b %d %Y %H:%M:%S %Z"),enkf.time.stamp),
          prec="text")

ncatt_put(ncf,var.cov,"Valid",
          attval=sprintf("Covariance valid for week tarting %s",
            format(tnext + 86400*(nweeks-1)*7),format="%Y-%m-%d %H:%M:%S UTC"),
          prec="text")

ncatt_put(ncf,var.cov,"Comment",
          "Scaling factor deviations for valid week were sampled stochastically from this covariance matrix.",
          prec="text")

nc_close(ncf)

cat(sprintf("%s Wrote %s.\n",self,parm.nc.out))

# write out observations results
obs.out.nc <- sprintf("%s/flasks/flask_enkf.%s_%s.nc",outdir,t0s,t1s)
dim.obs <- ncdim_def("obs","",vals=1:nobs,unlim=TRUE)
dim.char100 <- ncdim_def("char100","",vals=1:100,create_dimvar=FALSE)
var.r <- ncvar_def(name="r",units="ppm^2",dim=list(dim.obs),missval=-1e34,longname='MDM component (R) of obs error variance')
var.obspackid <- ncvar_def(name="obspack_id",units="",dim=list(dim.char100,dim.obs),missval=' ',longname='Observations identifier',prec="char")
var.hqhr <- ncvar_def(name="hqhr",units="ppm^2",dim=list(dim.obs),missval=-1e34,longname='HQH+R total obs error variance')
var.reject <- ncvar_def(name="rejected",units="NA",dim=list(dim.obs),missval=-1,longname='logical-was observation rejected')
var.dupe <- ncvar_def(name="duplicated",units="NA",dim=list(dim.obs),missval=-1,longname='logical-is observation a duplicate')
ncf <- nc_create(obs.out.nc,vars=list(var.reject,var.r,var.hqhr,var.dupe,var.obspackid))
ncvar_put(ncf,var.r,r)
ncvar_put(ncf,var.hqhr,hqhr)
ncvar_put(ncf,var.reject,as.integer(reject))
#-- restricting is.dupe to first ensemble member, not sure AES
ncvar_put(ncf,var.dupe,as.integer(is.dupe)[1:nobs])
ncvar_put(ncf,var.obspackid,sprintf("%100s",fout$obspack.id))
ncatt_put(ncf,0,"history",
          attval=sprintf("Created on %s\nby script '%s'",
            format(Sys.time(), "%a %b %d %Y %H:%M:%S %Z"),enkf.time.stamp),
          prec="text")

nc_close(ncf)

cat(sprintf("%s Wrote %s.\n",self,obs.out.nc))





        #stop("forced stop")

	#-- Single rerun of mean to get restart CO2 for next cycle
	print("launching mean beta rerun to generate new CO2 field ....")

       ###########################################
       #-- NASA (pods.sh), default NASA
       ###########################################

       if(que_soft == "nasa")
       {

 #working_dir = "/discover/nobackup/aschuh/run"
         reg.folder = "/discover/nobackup/aschuh/reg_folders/my_job_dir41"

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
           writeLines(paste("./geos ",k," ",cycles[i],
                               " ",rdate_arg," ",cycle_length," ",cycle_length," 1 > ",reg.folder,"/outfile.",k,sep=""),con=con)
          }

         close(con)

         #system(" /usr/local/other/PoDS/PoDS/pods.py -x /discover/nobackup/aschuh/reg_folders/my_job_dir41/exec.script -n 1")
         system("/discover/nobackup/aschuh/pods.sh /discover/nobackup/aschuh/reg_folders/my_job_dir41/exec.script 1")

        }

       ###########################################
       #-- End  NASA (pods.sh), default NASA
       ###########################################

     #stop("forced stop")

     ##############################################
     #--  End  SGE (sun grid engine)
     ##############################################

     ##############################################
     #-
     #-  Cleanup/removal of files
     #-
     ##############################################

     #-- Remove CT files, but leave FINALRUN files, this is important
     #-- because the model assumes that current cycle's
     #-- data is in this directory, regardless of name.  It simply lists and sorts files.

     modelfiles = list.files(paste(outdir,"/CT/",sep=""),full.names=TRUE)

     #-- We want to keep "FINALRUN" files as well as prior CO2 guesses for each cyle
     first_files  = list.files(paste(outdir,"/CT/",sep=""),full.names=TRUE,pattern="ens.0001")

     noremove_files  = list.files(paste(outdir,"/CT/",sep=""),full.names=TRUE,pattern="ens.0001.FINALRUN")

     prior_files = list.files(paste(outdir,"/CT/",sep=""),full.names=TRUE,pattern="ens.0001.PRIOR")

     move_files = first_files[ (!first_files %in% noremove_files) & (!first_files %in% prior_files)]

     if(length(move_files) > 0)
      {
                 for(k in 1:length(move_files))
             {
                  #newfilenm = paste(move_files[k],".PRIOR.nc",sep="")
                  newfilenm = gsub("ens.0001","ens.0001.PRIOR",move_files[k])
                  file.copy(move_files[k],newfilenm)
                  file.remove(move_files[k])
              }
       }

     remove_files = modelfiles[!modelfiles %in% first_files]

     file.remove(remove_files)

}
