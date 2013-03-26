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
	log1 = grepl(outdir,llines[grep("%%% ND50 MENU %%%",llines)+3])
	log2 = grepl(outdir,llines[grep("%%% GOSAT MENU %%%",llines)+3])
    log3 = grepl(outdir,llines[grep("Bias netcdf file",llines)])
    if(any(!c(log1,log2,log3)))
    {stop("'outdir' doesn't match output directories in geos input.geos file, stopping.")}
    log4 = grepl("prior",llines[grep("Bias netcdf file",llines)])    
    #if(prop_model == "ct" & !log4)
    #{stop("need to have 'prior' appended to end of bias file name in input.geos file, stopping.")}
}
