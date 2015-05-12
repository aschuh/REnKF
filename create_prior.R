#-- This function creates new a priori betas file from past 2 cycles of betas
#-- Variance is formed from averaging variations, not deviations and mean is currently
#-- just arithmetic mean, NOT variance weighted mean as it should be, let's fix that later
#-- ...

create_prior = function(ifiles,pr_ind,estimate_land_ocean_bias,ensembles,land_prior_scaling=1,ocean_prior_scaling=1)
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

                 #-- Scale down initial deviations
                 BETAGPP = BETAGPP*land_prior_scaling
                 BETARESP = BETARESP*land_prior_scaling
                 BETAOCEAN = BETAOCEAN*ocean_prior_scaling

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
          
               OBSLANDHBIAS = array(dim=c(ensembles))
               OBSLANDMBIAS = array(dim=c(ensembles))
               OBSOCEANBIAS =  array(dim=c(ensembles))
               OBSS32COEF = array(dim=c(ensembles))
               OBSB1_OFFSETCOEF =  array(dim=c(ensembles))
               OBSALBEDO_2_H_COEF =  array(dim=c(ensembles))
               OBSALBEDO_2_M_COEF =  array(dim=c(ensembles))
               OBSDP_CLD_COEF =  array(dim=c(ensembles))

           if(pr_ind[k]){
              #-- currently orig_betas_file has 1000 ensembles in it
              
              samp = c(1,sample(1:1000,ensembles-1))

               for(j in 1:ensembles){
                     ifiles.fil = nc_open(ifiles[k],readunlim=FALSE)
                     if(j==1){
                               #OBSLANDHBIAS = array(dim=c(ensembles))
                               #OBSLANDMBIAS = array(dim=c(ensembles))
                       	       #OBSOCEANBIAS =  array(dim=c(ensembles))
                               #OBSS32COEF = array(dim=c(ensembles))
                               #OBSB1_OFFSETCOEF =  array(dim=c(ensembles))
                               #OBSALBEDO_2_H_COEF =  array(dim=c(ensembles))
                               #OBSALBEDO_2_M_COEF =  array(dim=c(ensembles))
                               #OBSDP_CLD_COEF =  array(dim=c(ensembles))
                       	      }
                       OBSLANDHBIAS[j] = ncvar_get(ifiles.fil,"OBSLANDHBIAS",start=c(samp[j],1),count=c(1,1))#/length(ifiles)
                       OBSLANDMBIAS[j] = ncvar_get(ifiles.fil,"OBSLANDMBIAS",start=c(samp[j],1),count=c(1,1))#/length(ifiles)
                       OBSOCEANBIAS[j] = ncvar_get(ifiles.fil,"OBSOCEANBIAS",start=c(samp[j],1),count=c(1,1))#/length(ifiles)
                       OBSS32COEF[j] = ncvar_get(ifiles.fil,"OBSS32_COEF",start=c(samp[j],1),count=c(1,1))#/length(ifiles)
                       OBSB1_OFFSETCOEF[j] = ncvar_get(ifiles.fil,"OBSB1_OFFSETCOEF",start=c(samp[j],1),count=c(1,1))#/length(ifiles)
                       OBSALBEDO_2_H_COEF[j] = ncvar_get(ifiles.fil,"OBSALBEDO_2_H_COEF",start=c(samp[j],1),count=c(1,1))#/length(ifiles)
                       OBSALBEDO_2_M_COEF[j] = ncvar_get(ifiles.fil,"OBSALBEDO_2_M_COEF",start=c(samp[j],1),count=c(1,1))#/length(ifiles)
                       OBSDP_CLD_COEF[j] = ncvar_get(ifiles.fil,"OBSDP_CLD_COEF",start=c(samp[j],1),count=c(1,1))#/length(ifiles)

                     #-- I need to update the first vect element to be "zero" for the "grand" priors, "cold ens creation" script
                     if(j==1)
                       {
                          OBSLANDMBIAS[1] = 0;OBSLANDHBIAS[1] = 0;OBSOCEANBIAS[1] = 0;
                          OBSB1_OFFSETCOEF[1]=0;OBSALBEDO_2_M_COEF[1]=0;OBSS32COEF[1]=0;
                          OBSALBEDO_2_H_COEF[1]=0;OBSDP_CLD_COEF[1]=0;
                       }
                 }
                 
                 #-- Important, we need to convert deviations to variations before we weight, and add
                 OBSLANDHBIAS[2:ensembles] = OBSLANDHBIAS[2:ensembles]-rep(OBSLANDHBIAS[1],ensembles-1)
                 OBSLANDMBIAS[2:ensembles] = OBSLANDMBIAS[2:ensembles]-rep(OBSLANDMBIAS[1],ensembles-1)
                 OBSOCEANBIAS[2:ensembles] = OBSOCEANBIAS[2:ensembles]-rep(OBSOCEANBIAS[1],ensembles-1)
                 OBSS32COEF[2:ensembles] = OBSS32COEF[2:ensembles]-rep(OBSS32COEF[1],ensembles-1)
                 OBSB1_OFFSETCOEF[2:ensembles] = OBSB1_OFFSETCOEF[2:ensembles]-rep(OBSB1_OFFSETCOEF[1],ensembles-1)
                 OBSALBEDO_2_M_COEF[2:ensembles] = OBSALBEDO_2_M_COEF[2:ensembles]-rep(OBSALBEDO_2_M_COEF[1],ensembles-1)
                 OBSALBEDO_2_H_COEF[2:ensembles] = OBSALBEDO_2_H_COEF[2:ensembles]-rep(OBSALBEDO_2_H_COEF[1],ensembles-1)
                 OBSDP_CLD_COEF[2:ensembles] = OBSDP_CLD_COEF[2:ensembles]-rep(OBSDP_CLD_COEF[1],ensembles-1)

                 OBSLANDHBIAS = dev2var(OBSLANDHBIAS,divisor=divisor)
                 OBSLANDMBIAS = dev2var(OBSLANDMBIAS,divisor=divisor)
                 #OBSLANDBIAS[2:ensembles][OBSLANDBIAS[2:ensembles] <= 0] = ( -(OBSLANDBIAS[2:ensembles][OBSLANDBIAS[2:ensembles] <= 0] ) ^2 ) / divisor
                 #OBSLANDBIAS[2:ensembles][OBSLANDBIAS[2:ensembles] > 0] =  ( (OBSLANDBIAS[2:ensembles][OBSLANDBIAS[2:ensembles] > 0] ) ^2 ) / divisor

                 OBSOCEANBIAS = dev2var(OBSOCEANBIAS,divisor=divisor)
                 #OBSOCEANBIAS[2:ensembles][OBSOCEANBIAS[2:ensembles] <= 0] = ( -(OBSOCEANBIAS[2:ensembles][OBSOCEANBIAS[2:ensembles] <= 0] ) ^2 ) / divisor
                 #OBSOCEANBIAS[2:ensembles][OBSOCEANBIAS[2:ensembles] > 0] =  ( (OBSOCEANBIAS[2:ensembles][OBSOCEANBIAS[2:ensembles] > 0] ) ^2 ) / divisor

                 OBSS32COEF = dev2var(OBSS32COEF,divisor=divisor)
                 #OBSS32COEF[2:ensembles][OBSS32COEF[2:ensembles] <= 0] = ( -(OBSS32COEF[2:ensembles][OBSS32COEF[2:ensembles] <= 0] ) ^2 ) / divisor
                 #OBSS32COEF[2:ensembles][OBSS32COEF[2:ensembles] > 0] =  ( (OBSS32COEF[2:ensembles][OBSS32COEF[2:ensembles] > 0] ) ^2 ) / divisor

                 OBSB1_OFFSETCOEF = dev2var(OBSB1_OFFSETCOEF,divisor=divisor)
                 #OBSB1_OFFSETCOEF[2:ensembles][OBSB1_OFFSETCOEF[2:ensembles] <= 0] = ( -(OBSB1_OFFSETCOEF[2:ensembles][OBSB1_OFFSETCOEF[2:ensembles] <= 0] ) ^2 ) / divisor
                 #OBSB1_OFFSETCOEF[2:ensembles][OBSB1_OFFSETCOEF[2:ensembles] > 0] =  ( (OBSB1_OFFSETCOEF[2:ensembles][OBSB1_OFFSETCOEF[2:ensembles] > 0] ) ^2 ) / divisor

                 OBSALBEDO_2_H_COEF = dev2var(OBSALBEDO_2_H_COEF,divisor=divisor)
                 #OBSALBEDO_2COEF[2:ensembles][OBSALBEDO_2COEF[2:ensembles] <= 0] = ( -(OBSALBEDO_2COEF[2:ensembles][OBSALBEDO_2COEF[2:ensembles] <= 0] ) ^2 ) / divisor
                 #OBSALBEDO_2COEF[2:ensembles][OBSALBEDO_2COEF[2:ensembles] > 0] =  ( (OBSALBEDO_2COEF[2:ensembles][OBSALBEDO_2COEF[2:ensembles] > 0] ) ^2 ) / divisor

                 OBSALBEDO_2_M_COEF = dev2var(OBSALBEDO_2_M_COEF,divisor=divisor)
                 OBSDP_CLD_COEF = dev2var(OBSDP_CLD_COEF,divisor=divisor)
            }else{
                 ifiles.fil = nc_open(ifiles[k])
                 OBSLANDHBIAS = ncvar_get(ifiles.fil,"OBSLANDHBIAS")
                 OBSLANDMBIAS = ncvar_get(ifiles.fil,"OBSLANDMBIAS")
                 OBSOCEANBIAS = ncvar_get(ifiles.fil,"OBSOCEANBIAS")
                 OBSS32COEF = ncvar_get(ifiles.fil,"OBSS32_COEF")
                 OBSB1_OFFSETCOEF = ncvar_get(ifiles.fil,"OBSB1_OFFSETCOEF")
                 OBSALBEDO_2_H_COEF = ncvar_get(ifiles.fil,"OBSALBEDO_2_H_COEF")
                 OBSALBEDO_2_M_COEF = ncvar_get(ifiles.fil,"OBSALBEDO_2_M_COEF")
                 OBSDP_CLD_COEF = ncvar_get(ifiles.fil,"OBSDP_CLD_COEF")
                 nc_close(ifiles.fil)

                 #-- Important, we need to convert deviations to variations before we weight, and add
                 OBSLANDHBIAS[2:ensembles] = OBSLANDHBIAS[2:ensembles]-rep(OBSLANDHBIAS[1],ensembles-1)
                 OBSLANDMBIAS[2:ensembles] = OBSLANDMBIAS[2:ensembles]-rep(OBSLANDMBIAS[1],ensembles-1)
                 OBSOCEANBIAS[2:ensembles] = OBSOCEANBIAS[2:ensembles]-rep(OBSOCEANBIAS[1],ensembles-1)
                 OBSS32COEF[2:ensembles] = OBSS32COEF[2:ensembles]-rep(OBSS32COEF[1],ensembles-1)
                 OBSB1_OFFSETCOEF[2:ensembles] = OBSB1_OFFSETCOEF[2:ensembles]-rep(OBSB1_OFFSETCOEF[1],ensembles-1)
                 OBSALBEDO_2_H_COEF[2:ensembles] = OBSALBEDO_2_H_COEF[2:ensembles]-rep(OBSALBEDO_2_H_COEF[1],ensembles-1)
                 OBSALBEDO_2_M_COEF[2:ensembles] = OBSALBEDO_2_M_COEF[2:ensembles]-rep(OBSALBEDO_2_M_COEF[1],ensembles-1)
                 OBSDP_CLD_COEF[2:ensembles] = OBSDP_CLD_COEF[2:ensembles]-rep(OBSDP_CLD_COEF[1],ensembles-1)

                 OBSLANDHBIAS = dev2var(OBSLANDHBIAS,divisor=divisor)
                 OBSLANDMBIAS = dev2var(OBSLANDMBIAS,divisor=divisor)                                  
                 #OBSLANDBIAS[2:ensembles][OBSLANDBIAS[2:ensembles] <= 0] = ( -(OBSLANDBIAS[2:ensembles][OBSLANDBIAS[2:ensembles] <= 0] ) ^2 ) / divisor
                 #OBSLANDBIAS[2:ensembles][OBSLANDBIAS[2:ensembles] > 0] =  ( (OBSLANDBIAS[2:ensembles][OBSLANDBIAS[2:ensembles] > 0] ) ^2 ) / divisor

                 OBSOCEANBIAS = dev2var(OBSOCEANBIAS,divisor=divisor)
                 #OBSOCEANBIAS[2:ensembles][OBSOCEANBIAS[2:ensembles] <= 0] = ( -(OBSOCEANBIAS[2:ensembles][OBSOCEANBIAS[2:ensembles] <= 0] ) ^2 ) / divisor
                 #OBSOCEANBIAS[2:ensembles][OBSOCEANBIAS[2:ensembles] > 0] =  ( (OBSOCEANBIAS[2:ensembles][OBSOCEANBIAS[2:ensembles] > 0] ) ^2 ) / divisor

                 OBSS32COEF = dev2var(OBSS32COEF,divisor=divisor)
                 #OBSS32COEF[2:ensembles][OBSS32COEF[2:ensembles] <= 0] = ( -(OBSS32COEF[2:ensembles][OBSS32COEF[2:ensembles] <= 0] ) ^2 ) / divisor
                 #OBSS32COEF[2:ensembles][OBSS32COEF[2:ensembles] > 0] =  ( (OBSS32COEF[2:ensembles][OBSS32COEF[2:ensembles] > 0] ) ^2 ) / divisor

                 OBSB1_OFFSETCOEF = dev2var(OBSB1_OFFSETCOEF,divisor=divisor)
                 #OBSB1_OFFSETCOEF[2:ensembles][OBSB1_OFFSETCOEF[2:ensembles] <= 0] = ( -(B1_OFFSETCOEF[2:ensembles][B1_OFFSETCOEF[2:ensembles] <= 0] ) ^2 ) / divisor
                 #OBSB1_OFFSETCOEF[2:ensembles][OBSB1_OFFSETCOEF[2:ensembles] > 0] =  ( (B1_OFFSETCOEF[2:ensembles][B1_OFFSETCOEF[2:ensembles] > 0] ) ^2 ) / divisor

                 OBSALBEDO_2_H_COEF = dev2var(OBSALBEDO_2_H_COEF,divisor=divisor)
                 OBSALBEDO_2_M_COEF = dev2var(OBSALBEDO_2_M_COEF,divisor=divisor)
                 #OBSALBEDO_2COEF[2:ensembles][OBSALBEDO_2COEF[2:ensembles] <= 0] = ( -(OBSALBEDO_2COEF[2:ensembles][OBSALBEDO_2COEF[2:ensembles] <= 0] ) ^2 ) / divisor
                 #OBSALBEDO_2COEF[2:ensembles][OBSALBEDO_2COEF[2:ensembles] > 0] =  ( (OBSALBEDO_2COEF[2:ensembles][OBSALBEDO_2COEF[2:ensembles] > 0] ) ^2 ) / divisor
 
                 OBSDP_CLD_COEF = dev2var(OBSDP_CLD_COEF,divisor=divisor)
            }
          }else
          {
           if(pr_ind[k]){
               #-- For now, just use first 'prior' pull for both prior files, usually just first/second cycles
               OBSLANDHBIAS   = OBSLANDHBIAS   *2
               OBSLANDMBIAS   = OBSLANDMBIAS *2
               OBSOCEANBIAS  = OBSOCEANBIAS  *2
               OBSS32COEF   = OBSS32COEF   *2
               OBSB1_OFFSETCOEF  = OBSB1_OFFSETCOEF  *2
               OBSALBEDO_2_H_COEF   = OBSALBEDO_2_H_COEF   *2
               OBSALBEDO_2_M_COEF   = OBSALBEDO_2_M_COEF   *2
               OBSDP_CLD_COEF  = OBSDP_CLD_COEF * 2
           }else{
             ifiles.fil = nc_open(ifiles[k])
             OBSLANDHBIAS_NEW = ncvar_get(ifiles.fil,"OBSLANDHBIAS")#/length(ifiles)
             OBSLANDMBIAS_NEW = ncvar_get(ifiles.fil,"OBSLANDMBIAS")#/length(ifiles)
             OBSOCEANBIAS_NEW = ncvar_get(ifiles.fil,"OBSOCEANBIAS")#/length(ifiles)
             OBSS32COEF_NEW = ncvar_get(ifiles.fil,"OBSS32_COEF")#/length(ifiles)
             OBSB1_OFFSETCOEF_NEW = ncvar_get(ifiles.fil,"OBSB1_OFFSETCOEF")#/length(ifiles)
             OBSALBEDO_2_H_COEF_NEW = ncvar_get(ifiles.fil,"OBSALBEDO_2_H_COEF")#/length(ifiles)
             OBSALBEDO_2_M_COEF_NEW = ncvar_get(ifiles.fil,"OBSALBEDO_2_M_COEF")#/length(ifiles)
             OBSDP_CLD_COEF_NEW = ncvar_get(ifiles.fil,"OBSDP_CLD_COEF")#/length(ifiles)
             nc_close(ifiles.fil)

             #-- Important, we need to convert deviations to variations before we weight, and add

                 OBSLANDHBIAS_NEW[2:ensembles] = OBSLANDHBIAS_NEW[2:ensembles]-rep(OBSLANDHBIAS_NEW[1],ensembles-1)
                 OBSLANDMBIAS_NEW[2:ensembles] = OBSLANDMBIAS_NEW[2:ensembles]-rep(OBSLANDMBIAS_NEW[1],ensembles-1)
                 OBSOCEANBIAS_NEW[2:ensembles] = OBSOCEANBIAS_NEW[2:ensembles]-rep(OBSOCEANBIAS_NEW[1],ensembles-1)
                 OBSS32COEF_NEW[2:ensembles] = OBSS32COEF_NEW[2:ensembles]-rep(OBSS32COEF_NEW[1],ensembles-1)
                 OBSB1_OFFSETCOEF_NEW[2:ensembles] = OBSB1_OFFSETCOEF_NEW[2:ensembles]-rep(OBSB1_OFFSETCOEF_NEW[1],ensembles-1)
                 OBSALBEDO_2_H_COEF_NEW[2:ensembles] = OBSALBEDO_2_H_COEF_NEW[2:ensembles]-rep(OBSALBEDO_2_H_COEF_NEW[1],ensembles-1)
                 OBSALBEDO_2_M_COEF_NEW[2:ensembles] = OBSALBEDO_2_M_COEF_NEW[2:ensembles]-rep(OBSALBEDO_2_M_COEF_NEW[1],ensembles-1)
                 OBSDP_CLD_COEF_NEW[2:ensembles] = OBSDP_CLD_COEF_NEW[2:ensembles]-rep(OBSDP_CLD_COEF_NEW[1],ensembles-1)
                                  
                 OBSLANDHBIAS_NEW = dev2var(OBSLANDHBIAS_NEW,divisor=divisor)
                 OBSLANDMBIAS_NEW = dev2var(OBSLANDMBIAS_NEW,divisor=divisor)
                 #OBSLANDBIAS_NEW[2:ensembles][OBSLANDBIAS_NEW[2:ensembles] <= 0] = ( -(OBSLANDBIAS_NEW[2:ensembles][OBSLANDBIAS_NEW[2:ensembles] <= 0] ) ^2 ) / divisor
                 #OBSLANDBIAS_NEW[2:ensembles][OBSLANDBIAS_NEW[2:ensembles] > 0] =  ( (OBSLANDBIAS_NEW[2:ensembles][OBSLANDBIAS_NEW[2:ensembles] > 0] ) ^2 ) / divisor

                 OBSOCEANBIAS_NEW = dev2var(OBSOCEANBIAS_NEW,divisor=divisor)
                 #OBSOCEANBIAS_NEW[2:ensembles][OBSOCEANBIAS_NEW[2:ensembles] <= 0] = ( -(OBSOCEANBIAS_NEW[2:ensembles][OBSOCEANBIAS_NEW[2:ensembles] <= 0] ) ^2 ) / divisor
                 #OBSOCEANBIAS_NEW[2:ensembles][OBSOCEANBIAS_NEW[2:ensembles] > 0] =  ( (OBSOCEANBIAS_NEW[2:ensembles][OBSOCEANBIAS_NEW[2:ensembles] > 0] ) ^2 ) / divisor

                 OBSS32COEF_NEW = dev2var(OBSS32COEF_NEW,divisor=divisor)
                 #OBSS32COEF_NEW[2:ensembles][OBSS32COEF_NEW[2:ensembles] <= 0] = ( -(OBSS32COEF_NEW[2:ensembles][OBSS32COEF[2:ensembles] <= 0] ) ^2 ) / divisor
                 #OBSS32COEF_NEW[2:ensembles][OBSS32COEF_NEW[2:ensembles] > 0] =  ( (OBSS32COEF_NEW[2:ensembles][OBSS32COEF[2:ensembles] > 0] ) ^2 ) / divisor

                 OBSB1_OFFSETCOEF_NEW = dev2var(OBSB1_OFFSETCOEF_NEW,divisor=divisor)
                 #OBSB1_OFFSETCOEF_NEW[2:ensembles][OBSB1_OFFSETCOEF_NEW[2:ensembles] <= 0] = ( -(B1_OFFSETCOEF_NEW[2:ensembles][B1_OFFSETCOEF_NEW[2:ensembles] <= 0] ) ^2 ) / divisor
                 #OBSB1_OFFSETCOEF_NEW[2:ensembles][OBSB1_OFFSETCOEF_NEW[2:ensembles] > 0] =  ( (B1_OFFSETCOEF_NEW[2:ensembles][B1_OFFSETCOEF_NEW[2:ensembles] > 0] ) ^2 ) / divisor

                 OBSALBEDO_2_H_COEF_NEW = dev2var(OBSALBEDO_2_H_COEF_NEW,divisor=divisor)
                 OBSALBEDO_2_M_COEF_NEW = dev2var(OBSALBEDO_2_M_COEF_NEW,divisor=divisor)
                 #OBSALBEDO_2COEF_NEW[2:ensembles][OBSALBEDO_2COEF_NEW[2:ensembles] <= 0] = ( -(OBSALBEDO_2COEF_NEW[2:ensembles][OBSALBEDO_2COEF_NEW[2:ensembles] <= 0] ) ^2 ) / divisor
                 #OBSALBEDO_2COEF_NEW[2:ensembles][OBSALBEDO_2COEF_NEW[2:ensembles] > 0] =  ( (OBSALBEDO_2COEF_NEW[2:ensembles][OBSALBEDO_2COEF_NEW[2:ensembles] > 0] ) ^2 ) / divisor

                 OBSDP_CLD_COEF_NEW = dev2var(OBSDP_CLD_COEF_NEW,divisor=divisor)

             OBSLANDHBIAS = OBSLANDHBIAS + OBSLANDHBIAS_NEW
             OBSLANDMBIAS = OBSLANDMBIAS + OBSLANDMBIAS_NEW
             OBSOCEANBIAS = OBSOCEANBIAS + OBSOCEANBIAS_NEW
             OBSS32COEF = OBSS32COEF + OBSS32COEF_NEW
             OBSB1_OFFSETCOEF = OBSB1_OFFSETCOEF + OBSB1_OFFSETCOEF_NEW
             OBSALBEDO_2_H_COEF = OBSALBEDO_2_H_COEF + OBSALBEDO_2_H_COEF_NEW
             OBSALBEDO_2_M_COEF = OBSALBEDO_2_M_COEF + OBSALBEDO_2_M_COEF_NEW
             OBSDP_CLD_COEF = OBSDP_CLD_COEF + OBSDP_CLD_COEF_NEW
          }
         }
        }
       #-- Convert back into deviations from variations

       OBSLANDHBIAS = var2dev(OBSLANDHBIAS)
       OBSLANDMBIAS = var2dev(OBSLANDMBIAS)
       #OBSLANDBIAS[2:ensembles][OBSLANDBIAS[2:ensembles] <= 0] = ( -(-OBSLANDBIAS[2:ensembles][OBSLANDBIAS[2:ensembles] <= 0] ) ^0.5 )
       #OBSLANDBIAS[2:ensembles][OBSLANDBIAS[2:ensembles] > 0] =  ( (OBSLANDBIAS[2:ensembles][OBSLANDBIAS[2:ensembles] > 0] ) ^0.5 )

       OBSOCEANBIAS = var2dev(OBSOCEANBIAS)
       #OBSOCEANBIAS[2:ensembles][OBSOCEANBIAS[2:ensembles] <= 0] = ( -(-OBSOCEANBIAS[2:ensembles][OBSOCEANBIAS[2:ensembles] <= 0] ) ^0.5 )
       #OBSOCEANBIAS[2:ensembles][OBSOCEANBIAS[2:ensembles] > 0] =  ( (OBSOCEANBIAS[2:ensembles][OBSOCEANBIAS[2:ensembles] > 0] ) ^0.5 )

       OBSS32COEF = var2dev(OBSS32COEF)
       #OBSS32COEF[2:ensembles][OBSS32COEF[2:ensembles] <= 0] = ( -(-OBSS32COEF[2:ensembles][OBSS32COEF[2:ensembles] <= 0] ) ^0.5 )
       #OBSS32COEF[2:ensembles][OBSS32COEF[2:ensembles] > 0] =  ( (OBSS32COEF[2:ensembles][OBSS32COEF[2:ensembles] > 0] ) ^0.5 )

       OBSB1_OFFSETCOEF = var2dev(OBSB1_OFFSETCOEF)
       #OBSB1_OFFSETCOEF[2:ensembles][OBSB1_OFFSETCOEF[2:ensembles] <= 0] = ( -(-OBSB1_OFFSETCOEF[2:ensembles][OBSB1_OFFSETCOEF[2:ensembles] <= 0] ) ^0.5 )
       #OBSB1_OFFSETCOEF[2:ensembles][OBSB1_OFFSETCOEF[2:ensembles] > 0] =  ( (OBSB1_OFFSETCOEF[2:ensembles][OBSB1_OFFSETCOEF[2:ensembles] > 0] ) ^0.5 )

       OBSALBEDO_2_H_COEF = var2dev(OBSALBEDO_2_H_COEF)
       OBSALBEDO_2_M_COEF = var2dev(OBSALBEDO_2_M_COEF)
       #OBSALBEDO_2COEF[2:ensembles][OBSALBEDO_2COEF[2:ensembles] <= 0] = ( -(-OBSALBEDO_2COEF[2:ensembles][OBSALBEDO_2COEF[2:ensembles] <= 0] ) ^0.5 )
       #OBSALBEDO_2COEF[2:ensembles][OBSALBEDO_2COEF[2:ensembles] > 0] =  ( (OBSALBEDO_2COEF[2:ensembles][OBSALBEDO_2COEF[2:ensembles] > 0] ) ^0.5 )

       OBSDP_CLD_COEF = var2dev(OBSDP_CLD_COEF)

       #-- Adjust mean (which is currently just added together) by weighting
       OBSLANDHBIAS[1] = OBSLANDHBIAS[1]/length(ifiles)
       OBSLANDMBIAS[1] = OBSLANDMBIAS[1]/length(ifiles)
       OBSOCEANBIAS[1] = OBSOCEANBIAS[1]/length(ifiles)
       OBSS32COEF[1] = OBSS32COEF[1]/length(ifiles)
       OBSB1_OFFSETCOEF[1] = OBSB1_OFFSETCOEF[1]/length(ifiles)
       OBSALBEDO_2_H_COEF[1] = OBSALBEDO_2_H_COEF[1]/length(ifiles)
       OBSALBEDO_2_M_COEF[1] = OBSALBEDO_2_M_COEF[1]/length(ifiles)
       OBSDP_CLD_COEF[1] = OBSDP_CLD_COEF[1]/length(ifiles)

       #-- add back mean into deviations
       OBSLANDHBIAS[2:ensembles] = OBSLANDHBIAS[2:ensembles] + rep(OBSLANDHBIAS[1],ensembles-1)
       OBSLANDMBIAS[2:ensembles] = OBSLANDMBIAS[2:ensembles] + rep(OBSLANDMBIAS[1],ensembles-1)
       OBSOCEANBIAS[2:ensembles] = OBSOCEANBIAS[2:ensembles] + rep(OBSOCEANBIAS[1],ensembles-1)
       OBSS32COEF[2:ensembles] = OBSS32COEF[2:ensembles] + rep(OBSS32COEF[1],ensembles-1)
       OBSB1_OFFSETCOEF[2:ensembles] = OBSB1_OFFSETCOEF[2:ensembles] + rep(OBSB1_OFFSETCOEF[1],ensembles-1)
       OBSALBEDO_2_H_COEF[2:ensembles] = OBSALBEDO_2_H_COEF[2:ensembles] + rep(OBSALBEDO_2_H_COEF[1],ensembles-1)
       OBSALBEDO_2_M_COEF[2:ensembles] = OBSALBEDO_2_M_COEF[2:ensembles] + rep(OBSALBEDO_2_M_COEF[1],ensembles-1)
       OBSDP_CLD_COEF[2:ensembles] = OBSDP_CLD_COEF[2:ensembles] + rep(OBSDP_CLD_COEF[1],ensembles-1)       

       return(list(OBSLANDHBIAS=OBSLANDHBIAS, OBSLANDMBIAS=OBSLANDMBIAS,OBSOCEANBIAS=OBSOCEANBIAS,
                   OBSS32COEF=OBSS32COEF,OBSB1_OFFSETCOEF=OBSB1_OFFSETCOEF,OBSALBEDO_2_H_COEF=OBSALBEDO_2_H_COEF,
                   OBSALBEDO_2_M_COEF=OBSALBEDO_2_M_COEF,OBSDP_CLD_COEF=OBSDP_CLD_COEF))
	
}


dev2var = function(vect,divisor)
                 {
                         vect[2:ensembles][vect[2:ensembles] <= 0] = ( -(vect[2:ensembles][vect[2:ensembles] <= 0] ) ^2 ) / divisor
                         vect[2:ensembles][vect[2:ensembles] > 0] =  ( (vect[2:ensembles][vect[2:ensembles] > 0] ) ^2 ) / divisor
                         return(vect)
                 }

var2dev = function(vect)
                 {
                   vect[2:ensembles][vect[2:ensembles] <= 0] = ( -(-vect[2:ensembles][vect[2:ensembles] <= 0] ) ^0.5 )
                   vect[2:ensembles][vect[2:ensembles] > 0] =  ( (vect[2:ensembles][vect[2:ensembles] > 0] ) ^0.5 )
                   return(vect)
                 }
