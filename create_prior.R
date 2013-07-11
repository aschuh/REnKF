#-- This function creates new a priori betas file from past 2 cycles of betas
#-- Variance is formed from averaging variations, not deviations and mean is currently
#-- just arithmetic mean, NOT variance weighted mean as it should be, let's fix that later
#-- ...

create_prior = function(ifiles,pr_ind,estimate_land_ocean_bias,ensembles)
{
	for(k in 1:length(ifiles))
        {
         if(k==1)
          {
           if(pr_ind[k]){
              #-- currently orig_betas_file has 1000 ensembles in it
              samp = c(1,sample(1:1000,ensembles-1))
              for(j in 1:ensembles){
                     ifiles.fil = nc_open(ifiles[k],readunlim=FALSE)
                     if(j==1){
                       test = ncvar_get(ifiles.fil,"BETAGPP",start=c(1,1,1,1),count=c(-1,-1,1,1))
                       BETAGPP = array(dim=c(dim(test)[1],dim(test)[2],ensembles))
                       BETARESP = array(dim=c(dim(test)[1],dim(test)[2],ensembles))
                       BETAOCEAN = array(dim=c(dim(test)[1],dim(test)[2],ensembles))
                       
                             }
                       BETAGPP[,,j] = ncvar_get(ifiles.fil,"BETAGPP",start=c(1,1,samp[j],1),count=c(-1,-1,1,1))#/length(ifiles)
                       BETARESP[,,j] = ncvar_get(ifiles.fil,"BETARESP",start=c(1,1,samp[j],1),count=c(-1,-1,1,1))#/length(ifiles)
                       BETAOCEAN[,,j] = ncvar_get(ifiles.fil,"BETAOCEAN",start=c(1,1,samp[j],1),count=c(-1,-1,1,1))#/length(ifiles)
                       }
                 #-- Important, we need to convert deviations to variations before we weight, and add
                 BETAGPP[,,2:ensembles]   =  aperm(aaply(BETAGPP[,,2:ensembles],c(3),.fun=function(x){x-BETAGPP[,,1]}),c(2,3,1))
                 BETARESP[,,2:ensembles]  =  aperm(aaply(BETARESP[,,2:ensembles],c(3),.fun=function(x){x-BETARESP[,,1]}),c(2,3,1))
                 BETAOCEAN[,,2:ensembles] =  aperm(aaply(BETAOCEAN[,,2:ensembles],c(3),.fun=function(x){x-BETAOCEAN[,,1]}),c(2,3,1))

                 BETAGPP[,,2:ensembles][BETAGPP[,,2:ensembles] <= 0] = ( -(BETAGPP[,,2:ensembles][BETAGPP[,,2:ensembles] <= 0] ) ^2 ) / length(ifiles)
                 BETAGPP[,,2:ensembles][BETAGPP[,,2:ensembles] > 0] =  ( (BETAGPP[,,2:ensembles][BETAGPP[,,2:ensembles] > 0] ) ^2 ) / length(ifiles)

                 BETARESP[,,2:ensembles][BETARESP[,,2:ensembles] <= 0] = ( -(BETARESP[,,2:ensembles][BETARESP[,,2:ensembles] <= 0] ) ^2 ) / length(ifiles)
                 BETARESP[,,2:ensembles][BETARESP[,,2:ensembles] > 0] =  ( (BETARESP[,,2:ensembles][BETARESP[,,2:ensembles] > 0] ) ^2 ) / length(ifiles)

                 BETAOCEAN[,,2:ensembles][BETAOCEAN[,,2:ensembles] <= 0] = ( -(BETAOCEAN[,,2:ensembles][BETAOCEAN[,,2:ensembles] <= 0] ) ^2 ) / length(ifiles)
                 BETAOCEAN[,,2:ensembles][BETAOCEAN[,,2:ensembles] > 0] =  ( (BETAOCEAN[,,2:ensembles][BETAOCEAN[,,2:ensembles] > 0] ) ^2 ) / length(ifiles)
                 

            }else{
                 ifiles.fil = nc_open(ifiles[k])
                 BETAGPP = ncvar_get(ifiles.fil,"BETAGPP")
                 BETARESP = ncvar_get(ifiles.fil,"BETARESP")
                 BETAOCEAN = ncvar_get(ifiles.fil,"BETAOCEAN")
                 nc_close(ifiles.fil)

                 #-- Important, we need to convert deviations to variations before we weight, and add
                 BETAGPP[,,2:ensembles] = aperm(aaply(BETAGPP[,,2:ensembles],c(3),.fun=function(x){x-BETAGPP[,,1]}),c(2,3,1))
                 BETARESP[,,2:ensembles] = aperm(aaply(BETARESP[,,2:ensembles],c(3),.fun=function(x){x-BETARESP[,,1]}),c(2,3,1))
                 BETAOCEAN[,,2:ensembles] = aperm(aaply(BETAOCEAN[,,2:ensembles],c(3),.fun=function(x){x-BETAOCEAN[,,1]}),c(2,3,1))

                 BETAGPP[,,2:ensembles][BETAGPP[,,2:ensembles] <= 0] = ( -(BETAGPP[,,2:ensembles][BETAGPP[,,2:ensembles] <= 0] ) ^2 ) / length(ifiles)
                 BETAGPP[,,2:ensembles][BETAGPP[,,2:ensembles] > 0] =  ( (BETAGPP[,,2:ensembles][BETAGPP[,,2:ensembles] > 0] ) ^2 ) / length(ifiles)

                 BETARESP[,,2:ensembles][BETARESP[,,2:ensembles] <= 0] = ( -(BETARESP[,,2:ensembles][BETARESP[,,2:ensembles] <= 0] ) ^2 ) / length(ifiles)
                 BETARESP[,,2:ensembles][BETARESP[,,2:ensembles] > 0] =  ( (BETARESP[,,2:ensembles][BETARESP[,,2:ensembles] > 0] ) ^2 ) / length(ifiles)

                 BETAOCEAN[,,2:ensembles][BETAOCEAN[,,2:ensembles] <= 0] = ( -(BETAOCEAN[,,2:ensembles][BETAOCEAN[,,2:ensembles] <= 0] ) ^2 ) / length(ifiles)
                 BETAOCEAN[,,2:ensembles][BETAOCEAN[,,2:ensembles] > 0] =  ( (BETAOCEAN[,,2:ensembles][BETAOCEAN[,,2:ensembles] > 0] ) ^2 ) / length(ifiles)
            }
          }else
          {
           if(pr_ind[k]){
               #-- For now, just use first 'prior' pull for both prior files, usually just first/second cycles
               BETAGPP   = BETAGPP   *2
               BETARESP  = BETARESP  *2
               BETAOCEAN = BETAOCEAN *2
           }else{
             ifiles.fil = nc_open(ifiles[k])
             BETAGPP_NEW = ncvar_get(ifiles.fil,"BETAGPP")#/length(ifiles)
             BETARESP_NEW = ncvar_get(ifiles.fil,"BETARESP")#/length(ifiles)
             BETAOCEAN_NEW = ncvar_get(ifiles.fil,"BETAOCEAN")#/length(ifiles)
             nc_close(ifiles.fil)

             #-- Important, we need to convert deviations to variations before we weight, and add
             BETAGPP_NEW[,,2:ensembles] =  aperm(aaply(BETAGPP_NEW[,,2:ensembles],c(3),.fun=function(x){x-BETAGPP_NEW[,,1]}),c(2,3,1))
             BETARESP_NEW[,,2:ensembles] = aperm(aaply(BETARESP_NEW[,,2:ensembles],c(3),.fun=function(x){x-BETARESP_NEW[,,1]}),c(2,3,1))
             BETAOCEAN_NEW[,,2:ensembles] = aperm(aaply(BETAOCEAN_NEW[,,2:ensembles],c(3),.fun=function(x){x-BETAOCEAN_NEW[,,1]}),c(2,3,1))

             BETAGPP_NEW[,,2:ensembles][BETAGPP_NEW[,,2:ensembles] <= 0] = ( -(BETAGPP_NEW[,,2:ensembles][BETAGPP_NEW[,,2:ensembles] <= 0] ) ^2 ) / length(ifiles)
             BETAGPP_NEW[,,2:ensembles][BETAGPP_NEW[,,2:ensembles] > 0] =  ( (BETAGPP_NEW[,,2:ensembles][BETAGPP_NEW[,,2:ensembles] > 0] ) ^2 ) / length(ifiles)

             BETARESP_NEW[,,2:ensembles][BETARESP_NEW[,,2:ensembles] <= 0] = ( -(BETARESP_NEW[,,2:ensembles][BETARESP_NEW[,,2:ensembles] <= 0] ) ^2 ) / length(ifiles)
             BETARESP_NEW[,,2:ensembles][BETARESP_NEW[,,2:ensembles] > 0] =  ( (BETARESP_NEW[,,2:ensembles][BETARESP_NEW[,,2:ensembles] > 0] ) ^2 ) / length(ifiles)

             BETAOCEAN_NEW[,,2:ensembles][BETAOCEAN_NEW[,,2:ensembles] <= 0] = ( -(BETAOCEAN_NEW[,,2:ensembles][BETAOCEAN_NEW[,,2:ensembles] <= 0] ) ^2 ) / length(ifiles)
             BETAOCEAN_NEW[,,2:ensembles][BETAOCEAN_NEW[,,2:ensembles] > 0] =  ( (BETAOCEAN_NEW[,,2:ensembles][BETAOCEAN_NEW[,,2:ensembles] > 0] ) ^2 ) / length(ifiles)

             BETAGPP = BETAGPP + BETAGPP_NEW
             BETARESP = BETARESP + BETARESP_NEW
             BETAOCEAN = BETAOCEAN + BETAOCEAN_NEW
          }
         }
        }
       #-- Convert back into deviations from variations

       BETAGPP[,,2:ensembles][BETAGPP[,,2:ensembles] <= 0] = ( -(-BETAGPP[,,2:ensembles][BETAGPP[,,2:ensembles] <= 0] ) ^0.5 )
       BETAGPP[,,2:ensembles][BETAGPP[,,2:ensembles] > 0] =  ( (BETAGPP[,,2:ensembles][BETAGPP[,,2:ensembles] > 0] ) ^0.5 )

       BETARESP[,,2:ensembles][BETARESP[,,2:ensembles] <= 0] = ( -(-BETARESP[,,2:ensembles][BETARESP[,,2:ensembles] <= 0] ) ^0.5 )
       BETARESP[,,2:ensembles][BETARESP[,,2:ensembles] > 0] =  ( (BETARESP[,,2:ensembles][BETARESP[,,2:ensembles] > 0] ) ^0.5 )

       BETAOCEAN[,,2:ensembles][BETAOCEAN[,,2:ensembles] <= 0] = ( -(-BETAOCEAN[,,2:ensembles][BETAOCEAN[,,2:ensembles] <= 0] ) ^0.5 )
       BETAOCEAN[,,2:ensembles][BETAOCEAN[,,2:ensembles] > 0] =  ( (BETAOCEAN[,,2:ensembles][BETAOCEAN[,,2:ensembles] > 0] ) ^0.5 )

       #-- Adjust mean (which is currently just added together) by weighting
       BETAGPP[,,1] = BETAGPP[,,1]/length(ifiles)
       BETARESP[,,1] = BETARESP[,,1]/length(ifiles)
       BETAOCEAN[,,1] = BETAOCEAN[,,1]/length(ifiles)

       #-- add back mean into deviations
       BETAGPP[,,2:ensembles] = aperm(aaply(BETAGPP[,,2:ensembles],c(3),.fun=function(x){x + BETAGPP[,,1]}),c(2,3,1))
       BETARESP[,,2:ensembles] = aperm(aaply(BETARESP[,,2:ensembles],c(3),.fun=function(x){x + BETARESP[,,1]}),c(2,3,1))
       BETAOCEAN[,,2:ensembles] = aperm(aaply(BETAOCEAN[,,2:ensembles],c(3),.fun=function(x){x + BETAOCEAN[,,1]}),c(2,3,1))

       return(list(BETAGPP=BETAGPP,BETARESP=BETARESP,BETAOCEAN=BETAOCEAN))
	
}

create_prior_landoceanbias = function(ifiles,pr_ind,ensembles)
{
        divisor = length(ifiles)^2

	for(k in 1:length(ifiles))
        {
         if(k==1)
          {
           if(pr_ind[k]){
              #-- currently orig_betas_file has 1000 ensembles in it
              
              samp = c(1,sample(1:1000,ensembles-1))
              
               for(j in 1:ensembles){
                     ifiles.fil = nc_open(ifiles[k],readunlim=FALSE)
                     if(j==1){
                               OBSLANDBIAS = array(dim=c(ensembles))
                       	       OBSOCEANBIAS =  array(dim=c(ensembles))
                       	      }
                       OBSLANDBIAS[j] = ncvar_get(ifiles.fil,"OBSLANDBIAS",start=c(samp[j],1),count=c(1,1))#/length(ifiles)
                       OBSOCEANBIAS[j] = ncvar_get(ifiles.fil,"OBSOCEANBIAS",start=c(samp[j],1),count=c(1,1))#/length(ifiles)
         
                     #-- I need to update the first vect element to be "zero" for the "grand" priors, "cold ens creation" script
                     if(j==1)
                       {
                          OBSLANDBIAS[1] = 0;OBSOCEANBIAS[1] = 0;
                       }
                 }
                 
                 #-- Important, we need to convert deviations to variations before we weight, and add
                 OBSLANDBIAS[2:ensembles] = OBSLANDBIAS[2:ensembles]-rep(OBSLANDBIAS[1],ensembles-1)
                 OBSOCEANBIAS[2:ensembles] = OBSOCEANBIAS[2:ensembles]-rep(OBSOCEANBIAS[1],ensembles-1)
                                  
                 OBSLANDBIAS[2:ensembles][OBSLANDBIAS[2:ensembles] <= 0] = ( -(OBSLANDBIAS[2:ensembles][OBSLANDBIAS[2:ensembles] <= 0] ) ^2 ) / divisor
                 OBSLANDBIAS[2:ensembles][OBSLANDBIAS[2:ensembles] > 0] =  ( (OBSLANDBIAS[2:ensembles][OBSLANDBIAS[2:ensembles] > 0] ) ^2 ) / divisor

                 OBSOCEANBIAS[2:ensembles][OBSOCEANBIAS[2:ensembles] <= 0] = ( -(OBSOCEANBIAS[2:ensembles][OBSOCEANBIAS[2:ensembles] <= 0] ) ^2 ) / divisor
                 OBSOCEANBIAS[2:ensembles][OBSOCEANBIAS[2:ensembles] > 0] =  ( (OBSOCEANBIAS[2:ensembles][OBSOCEANBIAS[2:ensembles] > 0] ) ^2 ) / divisor

            }else{
                 ifiles.fil = nc_open(ifiles[k])
                 OBSLANDBIAS = ncvar_get(ifiles.fil,"OBSLANDBIAS")
                 OBSOCEANBIAS = ncvar_get(ifiles.fil,"OBSOCEANBIAS")
                 nc_close(ifiles.fil)

                 #-- Important, we need to convert deviations to variations before we weight, and add
                 OBSLANDBIAS[2:ensembles] = OBSLANDBIAS[2:ensembles]-rep(OBSLANDBIAS[1],ensembles-1)
                 OBSOCEANBIAS[2:ensembles] = OBSOCEANBIAS[2:ensembles]-rep(OBSOCEANBIAS[1],ensembles-1)
                                  
                 OBSLANDBIAS[2:ensembles][OBSLANDBIAS[2:ensembles] <= 0] = ( -(OBSLANDBIAS[2:ensembles][OBSLANDBIAS[2:ensembles] <= 0] ) ^2 ) / divisor
                 OBSLANDBIAS[2:ensembles][OBSLANDBIAS[2:ensembles] > 0] =  ( (OBSLANDBIAS[2:ensembles][OBSLANDBIAS[2:ensembles] > 0] ) ^2 ) / divisor

                 OBSOCEANBIAS[2:ensembles][OBSOCEANBIAS[2:ensembles] <= 0] = ( -(OBSOCEANBIAS[2:ensembles][OBSOCEANBIAS[2:ensembles] <= 0] ) ^2 ) / divisor
                 OBSOCEANBIAS[2:ensembles][OBSOCEANBIAS[2:ensembles] > 0] =  ( (OBSOCEANBIAS[2:ensembles][OBSOCEANBIAS[2:ensembles] > 0] ) ^2 ) / divisor
            }
          }else
          {
           if(pr_ind[k]){
               #-- For now, just use first 'prior' pull for both prior files, usually just first/second cycles
               OBSLANDBIAS   = OBSLANDBIAS   *2
               OBSOCEANBIAS  = OBSOCEANBIAS  *2
           }else{
             ifiles.fil = nc_open(ifiles[k])
             OBSLANDBIAS_NEW = ncvar_get(ifiles.fil,"OBSLANDBIAS")#/length(ifiles)
             OBSOCEANBIAS_NEW = ncvar_get(ifiles.fil,"OBSOCEANBIAS")#/length(ifiles)
             nc_close(ifiles.fil)

             #-- Important, we need to convert deviations to variations before we weight, and add
                 OBSLANDBIAS_NEW[2:ensembles] = OBSLANDBIAS_NEW[2:ensembles]-rep(OBSLANDBIAS_NEW[1],ensembles-1)
                 OBSOCEANBIAS_NEW[2:ensembles] = OBSOCEANBIAS_NEW[2:ensembles]-rep(OBSOCEANBIAS_NEW[1],ensembles-1)
                                  
                 OBSLANDBIAS_NEW[2:ensembles][OBSLANDBIAS_NEW[2:ensembles] <= 0] = ( -(OBSLANDBIAS_NEW[2:ensembles][OBSLANDBIAS_NEW[2:ensembles] <= 0] ) ^2 ) / divisor
                 OBSLANDBIAS_NEW[2:ensembles][OBSLANDBIAS_NEW[2:ensembles] > 0] =  ( (OBSLANDBIAS_NEW[2:ensembles][OBSLANDBIAS_NEW[2:ensembles] > 0] ) ^2 ) / divisor

                 OBSOCEANBIAS_NEW[2:ensembles][OBSOCEANBIAS_NEW[2:ensembles] <= 0] = ( -(OBSOCEANBIAS_NEW[2:ensembles][OBSOCEANBIAS_NEW[2:ensembles] <= 0] ) ^2 ) / divisor
                 OBSOCEANBIAS_NEW[2:ensembles][OBSOCEANBIAS_NEW[2:ensembles] > 0] =  ( (OBSOCEANBIAS_NEW[2:ensembles][OBSOCEANBIAS_NEW[2:ensembles] > 0] ) ^2 ) / divisor

             OBSLANDBIAS = OBSLANDBIAS + OBSLANDBIAS_NEW
             OBSOCEANBIAS = OBSOCEANBIAS + OBSOCEANBIAS_NEW
          }
         }
        }
       #-- Convert back into deviations from variations

       OBSLANDBIAS[2:ensembles][OBSLANDBIAS[2:ensembles] <= 0] = ( -(-OBSLANDBIAS[2:ensembles][OBSLANDBIAS[2:ensembles] <= 0] ) ^0.5 )
       OBSLANDBIAS[2:ensembles][OBSLANDBIAS[2:ensembles] > 0] =  ( (OBSLANDBIAS[2:ensembles][OBSLANDBIAS[2:ensembles] > 0] ) ^0.5 )

       OBSOCEANBIAS[2:ensembles][OBSOCEANBIAS[2:ensembles] <= 0] = ( -(-OBSOCEANBIAS[2:ensembles][OBSOCEANBIAS[2:ensembles] <= 0] ) ^0.5 )
       OBSOCEANBIAS[2:ensembles][OBSOCEANBIAS[2:ensembles] > 0] =  ( (OBSOCEANBIAS[2:ensembles][OBSOCEANBIAS[2:ensembles] > 0] ) ^0.5 )

       #-- Adjust mean (which is currently just added together) by weighting
       OBSLANDBIAS[1] = OBSLANDBIAS[1]/length(ifiles)
       OBSOCEANBIAS[1] = OBSOCEANBIAS[1]/length(ifiles)

       #-- add back mean into deviations
       OBSLANDBIAS[2:ensembles] = OBSLANDBIAS[2:ensembles] + rep(OBSLANDBIAS[1],ensembles-1)
       OBSOCEANBIAS[2:ensembles] = OBSOCEANBIAS[2:ensembles] + rep(OBSOCEANBIAS[1],ensembles-1)
       
       return(list(OBSLANDBIAS=OBSLANDBIAS,OBSOCEANBIAS=OBSOCEANBIAS))
	
}

