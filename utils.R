############################################################################
# Matrix manipulation methods
#
# For simplicity we have avoided to create generic functions for 'flip' etc.
# and therefore we have to call the corresponding methods coupled to the
# 'matrix' class explicitly, i.e. flip.matrix().
############################################################################
# Flip matrix (upside-down)
flip.matrix <- function(x) {
  mirror.matrix(rotate180.matrix(x))
}

# Mirror matrix (left-right)
mirror.matrix <- function(x) {
  xx <- as.data.frame(x);
  xx <- rev(xx);
  xx <- as.matrix(xx);
  xx;
}

# Rotate matrix 90 clockworks
rotate90.matrix <- function(x) {
  t(mirror.matrix(x))
}

# Rotate matrix 180 clockworks
rotate180.matrix <- function(x) { 
  xx <- rev(x);
  dim(xx) <- dim(x);
  xx;
}

# Rotate matrix 270 clockworks
rotate270.matrix <- function(x) {
  mirror.matrix(t(x))
}
###################
# Utility functions
###################

agg.infl2 <- function(mat,size,norm=T)
  {
    if(size==1)
    { return(as.vector(mat)) }

    else{
        retmat <- matrix(ncol=dim(mat)[2]/size,nrow=dim(mat)[1]/size)
        for(i in 1:(dim(mat)[1]/size))
        {
         for(j in 1:(dim(mat)[2]/size))
         {
             retmat[i,j] <- sum(mat[1:size+(i-1)*size,1:size+(j-1)*size])
         }
       }
     if(norm==T){retarg <- (retmat/size^2)}
     else{
       retarg <- retmat
     }
    return(retarg)
      }
  }
  
pad = function(x,width,fill=" ",left=TRUE)
{
        xneg=FALSE
        nc =  nchar(x)
        if(width>nc)
            {
              if(left){ str = paste(paste(rep(fill,width-nc),collapse=""),x,sep="")}
              else{ str = paste(x,paste(rep(fill,width-nc),collapse=""),sep="") }
            }
          else{str=x}
        return(str)
}

geos_inputfile_check = function(input_geos_file,outdir,prop_model)
{
	llines = readLines(input_geos_file)
	log1 = grepl(outdir,llines[grep("%%% ND50 MENU %%%",llines)+4])
	log2 = grepl(outdir,llines[grep("%%% GOSAT MENU %%%",llines)+3])
    log3 = grepl(outdir,llines[grep("Bias netcdf file",llines)])
    if(any(!c(log1,log2,log3)))
    {stop("'outdir' doesn't match output directories in geos input.geos file, stopping.")}
    log4 = grepl("prior",llines[grep("Bias netcdf file",llines)])    
    #if(prop_model == "ct" & !log4)
    #{stop("need to have 'prior' appended to end of bias file name in input.geos file, stopping.")}
}

betas2x25_to_betas4x5 = function(bmat_2x25)
{
	expand2 = 2
	d2 = aaply(.data=bmat_2x25,.margins=c(3),.fun=function(x){efun(x,expand2)})
	d2 = aperm(d2,c(2,3,1))
	d2 = d2[,-c(1:(expand2/2),((dim(d2)[2]-expand2/2+1):dim(d2)[2])),]

	d4 = aaply(.data=d2[,-c(1:2,179:180),],.margins=c(3),.fun=function(x){agg.infl2(x,4)})
	d4_new = aperm(d4,c(2,3,1))
	d4_part1 = abind(d2[,c(1:2),],d2[,c(1:2),],d2[,c(179:180),],d2[,c(179:180),],along=2)
	dimnames(d4_part1)=NULL

	dpart2 = aaply(.data=d4_part1,.margins=c(3),
	              .fun=function(x){agg.infl2(x,4)})
        dpart2 = aperm(dpart2,c(2,3,1))
    
        final = abind(dpart2[,1,],d4_new,dpart2[,2,],along=2)
        return(final)
}


betas4x5_to_betas2x25 = function(bmat_4x5)
{
        expand2 = 4
        d2 = aaply(.data=bmat_4x5,.margins=c(3),.fun=function(x){efun(x,expand2)})
        d2 = aperm(d2,c(2,3,1))
        
        d2_part1 = abind(d2[,c(1:(expand2)),],d2[,((dim(d2)[2]-expand2+1):dim(d2)[2]),],along=2)
        dimnames(d2_part1)=NULL
        ddd = aaply(.data=d2_part1,.margins=c(3),
                      .fun=function(x){agg.infl2(x,2)})
        d27 = aperm(ddd,c(2,3,1))              
                      
                      
        d2 = d2[,-c(1:(expand2),((dim(d2)[2]-expand2+1):dim(d2)[2])),]
        d2 = abind(d2_part1[,1:2,],d2,d2_part1[,3:4,],along=2)
        d2 = d2[,-c(1,dim(d2)[2]),]
        dimnames(d2)=NULL
        d3 = aaply(.data=d2,.margins=c(3),.fun=function(x){agg.infl2(x,2)})
        
        d33 = abind(ddd[,,1],d3,ddd[,,2],along=3)
        final = aperm(d33,c(2,3,1))

        return(final)
}

post_fluxes=function(dir_4x5,dir_2x25,dir_final,timesteps){

  require(plyr)
  require(abind)

  dir_4x5 = "/discover/nobackup/aschuh/GEOS-CHEM_output/v9.2_test_4x5/betas/"

  dir_2x25 = "/discover/nobackup/aschuh/GEOS-CHEM_output/v9.2_test_2x25/betas/"

  dir_final = "/discover/nobackup/aschuh/GEOS-CHEM_output/v9.2_test_4x5/betas/"

  #window_length = 16

  fls.4x5_last = paste(dir_4x5,"betas_cycle_post_001_",sapply(timesteps,FUN=pad,width=3,fill="0"),".nc",sep="")

  fls.4x5_first = paste(dir_4x5,"betas_cycle_post_regrid_006_",sapply(timesteps,FUN=pad,width=3,fill="0"),".nc",sep="")

  fls.2x25 = paste(dir_2x25,"betas_cycle_post_006_",sapply(timesteps,FUN=pad,width=3,fill="0"),".nc",sep="")

  for(i in 1:length(timesteps))
   {
   	   print(paste("working on step:",i))
   	
        fil = nc_open(fls.4x5_last[i])
        BETAGPP_4x5_last = ncvar_get(fil,"BETAGPP")
        BETARESP_4x5_last = ncvar_get(fil,"BETARESP")
        BETAOCEAN_4x5_last = ncvar_get(fil,"BETAOCEAN")
        nc_close(fil)

        fil = nc_open(fls.4x5_first[i])
        BETAGPP_4x5_first = ncvar_get(fil,"BETAGPP")
        BETARESP_4x5_first = ncvar_get(fil,"BETARESP")
        BETAOCEAN_4x5_first = ncvar_get(fil,"BETAOCEAN")
        nc_close(fil)

        fil = nc_open(fls.2x25[i])
        BETAGPP_2x25 = ncvar_get(fil,"BETAGPP")
        BETARESP_2x25 = ncvar_get(fil,"BETARESP")
        BETAOCEAN_2x25 = ncvar_get(fil,"BETAOCEAN")
        nc_close(fil)

    BETAGPP_2x25_FINAL = betas4x5_to_betas2x25(BETAGPP_4x5_last-BETAGPP_4x5_first) + BETAGPP_2x25
    BETARESP_2x25_FINAL = betas4x5_to_betas2x25(BETARESP_4x5_last-BETARESP_4x5_first) + BETARESP_2x25
    BETAOCEAN_2x25_FINAL = betas4x5_to_betas2x25(BETAOCEAN_4x5_last-BETAOCEAN_4x5_first) + BETAOCEAN_2x25

     #--- CREATE NETCDF FILE

     ens.size = dim(BETAGPP_2x25_FINAL)[3]
     x_res_out = 2.5
     y_res_out = 2

     # Define some straightforward dimensions
      x  = ncdim_def( "lon", "degrees_east", as.double(seq(-180,180-x_res_out,by=x_res_out)))
      y  = ncdim_def( "lat", "degrees_north", as.double(c(-90+0.25*y_res_out,seq((-90+0.25*y_res_out)+0.75*y_res_out,
                                                         90-0.75*y_res_out,by=y_res_out),90-0.25*y_res_out)))
     #z  = ncdim_def("Levels","hPa",as.double(1:47))
     t  = ncdim_def( "time", "hours since 1900-01-01", 1, unlim=TRUE)


     # Create a netCDF file with this variable

    varlist = list()

    nens =  ncdim_def( "ensemble", units="integer", as.double(seq(1:(ens.size))))

    varlist[[1]]  <- ncvar_def(name="BETAGPP",units="",
                                 dim=list(x,y,nens,t), missval=NA,prec="float")

    varlist[[2]]  <- ncvar_def(name="BETARESP",units="",
                                 dim=list(x,y,nens,t), missval=NA,prec="float")

    varlist[[3]]  <- ncvar_def(name="BETAOCEAN",units="",
                                 dim=list(x,y,nens,t), missval=NA,prec="float")

    ncnew <- nc_create(paste(dir_final,"betas_final_",pad(timesteps[i],width=3,fill="0"),".nc",sep=""),varlist)


    ncvar_put(ncnew, varlist[[1]],BETAGPP_2x25_FINAL)
    ncvar_put(ncnew, varlist[[2]],BETARESP_2x25_FINAL)
    ncvar_put(ncnew, varlist[[3]],BETAOCEAN_2x25_FINAL)

    nc_close(ncnew)
   }

}

post_fluxes_old=function(dir_4x5,dir_2x25,dir_final,window_length){

  #dir_4x5 = "/discover/nobackup/aschuh/GEOS-CHEM_output/v9.2_test_4x5/betas/"

  #dir_2x25 = "/discover/nobackup/aschuh/GEOS-CHEM_output/v9.2_test_2x25/betas/"

  #dir_final = "/discover/nobackup/aschuh/GEOS-CHEM_output/v9.2_test_4x5/betas/"

  #window_length = 16

  fls.4x5_last = paste(dir_4x5,"betas_cycle_post_001_",sapply(1:window_length,FUN=pad,width=3,fill="0"),".nc",sep="")

  fls.4x5_first = paste(dir_4x5,"betas_cycle_post_regrid_006_",sapply(1:window_length,FUN=pad,width=3,fill="0"),".nc",sep="")

  fls.2x25 = paste(dir_2x25,"betas_cycle_post_006_",sapply(1:window_length,FUN=pad,width=3,fill="0"),".nc",sep="")

  for(i in 1:length(window_length))
   {
	fil = nc_open(fls.4x5_last[i])
	BETAGPP_4x5_last = ncvar_get(fil,"BETAGPP")
	BETARESP_4x5_last = ncvar_get(fil,"BETARESP")
	BETAOCEAN_4x5_last = ncvar_get(fil,"BETAOCEAN")
	nc_close(fil)
	
	fil = nc_open(fls.4x5_first[i])
	BETAGPP_4x5_first = ncvar_get(fil,"BETAGPP")
	BETARESP_4x5_first = ncvar_get(fil,"BETARESP")
	BETAOCEAN_4x5_first = ncvar_get(fil,"BETAOCEAN")
	nc_close(fil)
	
	fil = nc_open(fls.2x25[i])
	BETAGPP_2x25 = ncvar_get(fil,"BETAGPP")
	BETARESP_2x25 = ncvar_get(fil,"BETARESP")
	BETAOCEAN_2x25 = ncvar_get(fil,"BETAOCEAN")
	nc_close(fil)

    BETAGPP_2x25_FINAL = betas4x5_to_betas2x25(BETAGPP_4x5_last-BETAGPP_4x5_first) + BETAGPP_2x25
    BETARESP_2x25_FINAL = betas4x5_to_betas2x25(BETARESP_4x5_last-BETARESP_4x5_first) + BETARESP_2x25
    BETAOCEAN_2x25_FINAL = betas4x5_to_betas2x25(BETAOCEAN_4x5_last-BETAOCEAN_4x5_first) + BETAOCEAN_2x25

     #--- CREATE NETCDF FILE

     # Define some straightforward dimensions
      x  = ncdim_def( "lon", "degrees_east", as.double(seq(-180,180-x_res_out,by=x_res_out)))
      y  = ncdim_def( "lat", "degrees_north", as.double(c(-90+0.25*y_res_out,seq((-90+0.25*y_res_out)+0.75*y_res_out,
                                                         90-0.75*y_res_out,by=y_res_out),90-0.25*y_res_out)))
     #z  = ncdim_def("Levels","hPa",as.double(1:47))
     t  = ncdim_def( "time", "hours since 1900-01-01", 1, unlim=TRUE)


     # Create a netCDF file with this variable

    varlist = list()

    nens =  ncdim_def( "ensemble", units="integer", as.double(seq(1:(ens.size))))

    varlist[[1]]  <- ncvar_def(name="BETAGPP",units="",
                                 dim=list(x,y,nens,t), missval=NA,prec="float")

    varlist[[2]]  <- ncvar_def(name="BETARESP",units="",
                                 dim=list(x,y,nens,t), missval=NA,prec="float")

    varlist[[3]]  <- ncvar_def(name="BETAOCEAN",units="",
                                 dim=list(x,y,nens,t), missval=NA,prec="float")

    ncnew <- nc_create(paste(dir_final,"betas_final_",pad(i,width=3,fill="0"),".nc",sep=""),varlist)


    ncvar_put(ncnew, varlist[[1]],bgpp_4x5)
    ncvar_put(ncnew, varlist[[2]],bresp_4x5)
    ncvar_put(ncnew, varlist[[3]],bocean_4x5)

    nc_close(ncnew)
   }	
	
}

efun = function(m,expand=2){
	orig.x = dim(m)[1]
	orig.y = dim(m)[2]
	dat = sapply(1:expand,FUN=function(x){as.vector(m)})
	datm = matrix(as.vector(t(dat)),ncol=orig.y,byrow=FALSE)
	datmm = t(datm)
	dat = sapply(1:expand,FUN=function(x){as.vector(datmm)})
	datmmm = matrix(as.vector(t(dat)),ncol=orig.x*expand,byrow=FALSE)	
	datmmm = t(datmmm)
	return(datmmm)
}

regrid_2x25_betas=function(ifile,ofile,x_res_out=5,y_res_out=4)
{
 require(abind)
 require(ncdf4)
 fil = nc_open(ifile)

 bgpp_2x25 = ncvar_get(fil,"BETAGPP")
 bresp_2x25 = ncvar_get(fil,"BETARESP")
 bocean_2x25 = ncvar_get(fil,"BETAOCEAN")  

  ens.size = dim(bgpp_2x25)[3]

 bgpp_4x5 = betas2x25_to_betas4x5(bgpp_2x25)
 bresp_4x5 = betas2x25_to_betas4x5(bresp_2x25)
 bocean_4x5 = betas2x25_to_betas4x5(bocean_2x25)
  
  
 #--- CREATE NETCDF FILE

 # Define some straightforward dimensions
  x  = ncdim_def( "lon", "degrees_east", as.double(seq(-180,180-x_res_out,by=x_res_out)))
  y  = ncdim_def( "lat", "degrees_north", as.double(c(-90+0.25*y_res_out,seq((-90+0.25*y_res_out)+0.75*y_res_out,
                                                         90-0.75*y_res_out,by=y_res_out),90-0.25*y_res_out)))
  #z  = ncdim_def("Levels","hPa",as.double(1:47))
  t  = ncdim_def( "time", "hours since 1900-01-01", 1, unlim=TRUE)


  # Create a netCDF file with this variable

varlist = list()

nens =  ncdim_def( "ensemble", units="integer", as.double(seq(1:(ens.size))))

varlist[[1]]  <- ncvar_def(name="BETAGPP",units="",
                                 dim=list(x,y,nens,t), missval=NA,prec="float")

varlist[[2]]  <- ncvar_def(name="BETARESP",units="",
                                 dim=list(x,y,nens,t), missval=NA,prec="float")

varlist[[3]]  <- ncvar_def(name="BETAOCEAN",units="",
                                 dim=list(x,y,nens,t), missval=NA,prec="float")
                                                             
ncnew <- nc_create(ofile,varlist)


ncvar_put(ncnew, varlist[[1]],bgpp_4x5)
ncvar_put(ncnew, varlist[[2]],bresp_4x5)
ncvar_put(ncnew, varlist[[3]],bocean_4x5)

    nc_close(ncnew)

}

#-- Mem printing 
# improved list of objects
.ls.objects <- function (pos = 1, pattern, order.by,
                        decreasing=FALSE, head=FALSE, n=5) {
    napply <- function(names, fn) sapply(names, function(x)
                                         fn(get(x, pos = pos)))
    names <- ls(pos = pos, pattern = pattern)
    obj.class <- napply(names, function(x) as.character(class(x))[1])
    obj.mode <- napply(names, mode)
    obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
    obj.prettysize <- napply(names, function(x) {
                           capture.output(format(utils::object.size(x), units = "auto")) })
    obj.size <- napply(names, object.size)
    obj.dim <- t(napply(names, function(x)
                        as.numeric(dim(x))[1:2]))
    vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
    obj.dim[vec, 1] <- napply(names, length)[vec]
    out <- data.frame(obj.type, obj.size, obj.prettysize, obj.dim)
    names(out) <- c("Type", "Size", "PrettySize", "Rows", "Columns")
    if (!missing(order.by))
        out <- out[order(out[[order.by]], decreasing=decreasing), ]
    if (head)
        out <- head(out, n)
    out
}

# shorthand
lsos <- function(..., n=10) {
    .ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
}

showMemoryUse <- function(sort="size", decreasing=FALSE, limit) {

  objectList <- ls(parent.frame())

  oneKB <- 1024
  oneMB <- 1048576
  oneGB <- 1073741824

  memoryUse <- sapply(objectList, function(x) as.numeric(object.size(eval(parse(text=x)))))

  memListing <- sapply(memoryUse, function(size) {
        if (size >= oneGB) return(paste(round(size/oneGB,2), "GB"))
        else if (size >= oneMB) return(paste(round(size/oneMB,2), "MB"))
        else if (size >= oneKB) return(paste(round(size/oneKB,2), "kB"))
        else return(paste(size, "bytes"))
      })

  memListing <- data.frame(objectName=names(memListing),memorySize=memListing,row.names=NULL)

  if (sort=="alphabetical") memListing <- memListing[order(memListing$objectName,decreasing=decreasing),] 
  else memListing <- memListing[order(memoryUse,decreasing=decreasing),] #will run if sort not specified or "size"

  if(!missing(limit)) memListing <- memListing[1:limit,]

  print(memListing, row.names=FALSE)
  return(invisible(memListing))
}

cleanMem <- function(n=10) { for (i in 1:n) gc() }
