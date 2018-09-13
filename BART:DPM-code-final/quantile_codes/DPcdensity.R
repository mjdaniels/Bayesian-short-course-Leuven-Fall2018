##################################################################
# This program has been modified to estimate 
#    marginal distribution of potential outcomes
##################################################################

### DPcdencity.R                   
### Fit a Dirichlet Process Mixture of Normals model for conditional density
### estimation.
###
### Copyright: Alejandro Jara, 2008-2012.
###
### Last modification: 27-08-2010.
###
### This program is free software; you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation; either version 2 of the License, or (at
### your option) any later version.
###
### This program is distributed in the hope that it will be useful, but
### WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
### General Public License for more details.
###
### You should have received a copy of the GNU General Public License
### along with this program; if not, write to the Free Software
### Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
###
### The author's contact information:
###
###      Alejandro Jara
###      Department of Statistics
###      Facultad de Matematicas
###      Pontificia Universidad Catolica de Chile
###      Casilla 306, Correo 22 
###      Santiago
###      Chile
###      Voice: +56-2-3544506  URL  : http://www.mat.puc.cl/~ajara
###      Fax  : +56-2-3547729  Email: atjara@uc.cl
###


"DPcdensity"<-
function(y,x,xpred,ngrid=100,grid=NULL,compute.band=FALSE,type.band="PD",prior,mcmc,state,status,data=sys.frame(sys.parent()),work.dir=NULL)
UseMethod("DPcdensity")

DPcdensity.default<-
function(y,
         x,
         xpred,
         ngrid=100,
		 grid=NULL,
		 compute.band=FALSE,
		 type.band="PD",
         prior,
         mcmc,
         state,
         status,
         data=sys.frame(sys.parent()),
         work.dir=NULL)
{
       #########################################################################################
       # call parameters
       #########################################################################################
         cl <- match.call()

       #########################################################################################
       # response
       #########################################################################################
         nrec <- length(y)
          
       #########################################################################################
       # add covariates to the data matrix
       #########################################################################################
         x <- as.matrix(x)
         nx <- dim(x)[2]
         z <- cbind(y,x)
         nvar <- nx + 1  

       #########################################################################################
       # identify missing values to be imputed
       #########################################################################################
       nmiss <- sum(is.na(z))
       nmissi <- 1
       missp <- NULL
       for(i in 1:nrec)
       {
         nmi <- sum(is.na(z[i,]))
         nmi2 <- seq(1,nvar)[is.na(z[i,])] 
         if(nmi>0)
         {
           for(j in 1:nmi)
           { 
             missp <- rbind(missp,c(i,nmi2[j]))  
           }
         }  
       }
       
       if(nmiss==0)
       {
         nmissi <- 0
         nmiss <- 1
         missp <- matrix(0,nrow=nmiss,2)
       }
	   if(nmissi==1){
	   z[,2]<-rep(0,length(z[,2]))
	   }

       #########################################################################################
       # change working directory (if requested..)
       #########################################################################################
         
		 if(!is.null(work.dir)){
            cat("\n Changing working directory to ",work.dir,"\n")
            old.dir <- getwd()  # by default work in current working directory
            setwd(work.dir)
         }

       #########################################################################################
       # prediction
       #########################################################################################
         xpred <- as.matrix(xpred)
         npred <- nrow(xpred)
         pz <- ncol(xpred)  

         if(pz != nx)
         {
            cat("\n *** Error: ncol(xpred) != ncol(x). Need to match.\n")
            return(-1)
         }
 
		 if(is.null(grid))
		 { 
			yy <- na.omit(y)  
			miny <- min(yy)
			maxy <- max(yy)
			vary <- var(yy)
			grid <- seq(from=miny-0.25*sqrt(vary),to=maxy+0.25*sqrt(vary),length.out=ngrid)
		 }
		 else
		 {
			grid <- as.vector(grid)
			ngrid <- length(grid)
		 }

		 cband <- 0
		 if(compute.band)
		 {
			cband <- 1
		 }
			
		 tband <- 1
		 if(type.band!="HPD")
		 {
			tband <- 2
		 }


       #########################################################################################
       # prior information
       #########################################################################################

         if(is.null(prior$a0))      
		 {
            a0 <--1
			b0 <--1 
            alpha<-prior$alpha
            alpharand<-0
         }
         else
         {
            a0 <- prior$a0
            b0 <- prior$b0
			alpha <- 1
			alpharand <- 1
		 }
         a0b0 <- c(a0,b0)
  	 
		 if(is.null(prior$nu2))
		 {
            psiinv1 <- matrix(prior$psiinv1,nvar,nvar)
            psiinv2 <- psiinv1
            psi1 <- matrix(solve(psiinv1),nvar,nvar)
            nuvec <- c(prior$nu1,-1)
            psi1rand <- 0
		 }
		 else
		 {
			psiinv1 <- prior$s2*2
            psi1 <- matrix(solve(psiinv1),nvar,nvar)
            psiinv2 <- matrix(prior$psiinv2,nvar,nvar)
			nuvec <- c(prior$nu1,prior$nu2)
			psi1rand <- 1
		 }
  	 
		if(is.null(prior$m2) && is.null(prior$s2))
		{
			s2inv <- matrix(0,nrow=nvar,ncol=nvar) 
			s2invm2 <- matrix(0,nrow=nvar,ncol=1)
			m1 <- prior$m1
			m1rand <- 0
		}
		else
		{
            s2inv <- solve(prior$s2)
            s2invm2 <- s2inv%*%prior$m2
			m1 <- prior$m2

            m1rand <- 1
         }     

         if(is.null(prior$tau1) && is.null(prior$tau2)) 
         {
            tau <- c(-2,-2)
            k0 <- prior$k0
            k0rand <- 0
         }
         else
         {
            tau <- c(prior$tau1,prior$tau2)
            k0 <- rgamma(1,shape=prior$tau1,scale=prior$tau2)
            k0rand <- 1
         }

       #########################################################################################
       # mcmc specification
       #########################################################################################
         mcmcvec <- c(mcmc$nburn,mcmc$nskip,mcmc$ndisplay,cband,tband)
         nsave <- mcmc$nsave

       #########################################################################################
       # output
       #########################################################################################
         cpo <- matrix(0,nrow=nrec,ncol=2) 
         thetasave <- matrix(0,nrow=nsave,ncol=nvar+nvar*(nvar+1)/2+3)
		 randsave <- matrix(0,nrow=nsave,ncol=(nrec+2)*nvar+(nrec+1)*nvar*(nvar+1)/2)
         denspm <- matrix(0,nrow=nsave,ncol=ngrid)
         denspl <- matrix(0,nrow=nsave,ncol=ngrid)
         densph <- matrix(0,nrow=npred,ncol=ngrid)
		densp0 <- matrix(0,nrow=nsave,ncol=npred)
		     meanfpm <- matrix(0,nrow=npred,ncol=ngrid)
		     meanfpl <- matrix(0,nrow=npred,ncol=ngrid)
		     meanfph <- matrix(0,nrow=npred,ncol=ngrid)
		     meanfpm1 <- matrix(0,nrow=npred,ncol=ngrid)
		     meanfpl1 <- matrix(0,nrow=npred,ncol=ngrid)
		     meanfph1 <- matrix(0,nrow=npred,ncol=ngrid)
 

       #########################################################################################
       # parameters depending on status
       #########################################################################################
         nuniq <- nvar*(nvar+1)/2

         if(status==TRUE)
		 {
            muclus <- matrix(0,nrow=nrec+100,ncol=nvar)
            sigmaclus <- matrix(0,nrow=nrec+100,ncol=nuniq)
            for(i in 1:1)
            {
                counter <- 0
                for(j in 1:nvar)
                {
                    muclus[i,j] <- m1[j]
                    for(k in j:nvar)
                    {
          		counter<-counter+1
          		sigmaclus[i,counter]<-psiinv1[j,k]
                    }
                }
             }
             ncluster <- 1
             ss <- rep(1,nrec)
		 }
	 
      	 if(status==FALSE)
		 {
			alpha <- state$alpha
            m1 <- state$m1
            muclus <- state$muclus 
			ncluster <- state$ncluster
			psi1 <- state$psi1
			psiinv1 <- solve(psi1)
			k0 <- state$k0
			sigmaclus <- state$sigmaclus 
			ss <- state$ss
            if(nmissi==1) z <- state$z
		 }    

       #########################################################################################
       # working space
       #########################################################################################
         ccluster <- rep(0,nrec)
         cstrt <- matrix(0,nrow=nrec,ncol=nrec) 
         num <- matrix(0,nrow=npred,ncol=ngrid)
         denom <- rep(0,npred)  
         iflag <- rep(0,nvar) 
         muwork <- rep(0,nvar) 
         muwork2 <- rep(0,nvar) 
         prob <- rep(0,(nrec+100))
         s1 <- matrix(0,nvar,nvar)
         seed1 <- sample(1:29000,1)
         seed2 <- sample(1:29000,1)
         seed <- c(seed1,seed2)
         sigmawork <- matrix(0,nrow=nvar,ncol=nvar)
         sigmawork2 <- matrix(0,nrow=nvar,ncol=nvar)
         sigworkinv <- matrix(0,nrow=nvar,ncol=nvar)
         theta <- rep(0,nvar)
         workm1 <- matrix(0,nrow=nvar,ncol=nvar)
         workm2 <- matrix(0,nrow=nvar,ncol=nvar)
         workm3 <- matrix(0,nrow=nvar,ncol=nvar)
         workmh1 <- rep(0,nvar*(nvar+1)/2) 
         workmh2 <- rep(0,nvar*(nvar+1)/2) 
         workv1 <- rep(0,nvar) 
         workv2 <- rep(0,nvar) 
         workv3 <- rep(0,nvar) 
		 ywork <- rep(0,nvar)

         iflagx <- rep(0,nx)
         workvx <- rep(0,nx) 
         workmx <- matrix(0,nrow=nx,ncol=nx) 

         fs <- rep(0,ngrid) 
         fm <- rep(0,npred) 

         worksam <- rep(0,nsave) 

         
         numcpo <- rep(0,nrec)
         denomcpo <- rep(0,nrec)
         
         nuvec <- c(nuvec,m1rand)

       #########################################################################################
       # calling the fortran code
       #########################################################################################


         foo <- .Fortran("dpdenregr",
				nrec       =as.integer(nrec),
                nx         =as.integer(nx),
				nvar       =as.integer(nvar),
				z          =as.double(z),
				npred      =as.integer(npred),
				xpred      =as.double(xpred),
				ngrid      =as.integer(ngrid),
				grid       =as.double(grid),
				a0b0       =as.double(a0b0),
				k0         =as.double(k0),
				nuvec      =as.integer(nuvec),
				s2inv      =as.double(s2inv),
				s2invm2    =as.double(s2invm2),
				psiinv2    =as.double(psiinv2),
				tau        =as.double(tau),
				mcmc       =as.integer(mcmcvec),
				nsave      =as.integer(nsave),
				cpo        =as.double(cpo),
				thetasave  =as.double(thetasave),
				randsave   =as.double(randsave),
        densp0     =as.double(densp0),
				denspm     =as.double(denspm),
				denspl     =as.double(denspl),
				meanfpm    =as.double(meanfpm),
				meanfpl    =as.double(meanfpl),
				meanfph    =as.double(meanfph),
				meanfpm1    =as.double(meanfpm1),
				meanfpl1    =as.double(meanfpl1),
				meanfph1    =as.double(meanfph1),
				alpha      =as.double(alpha),		
				m1         =as.double(m1),		
                muclus     =as.double(muclus),		 		
				ncluster   =as.integer(ncluster),
				psi1       =as.double(psi1),
				psiinv1    =as.double(psiinv1),
				s1         =as.double(s1),
				sigmaclus  =as.double(sigmaclus),
				ss         =as.integer(ss),
				ccluster   =as.integer(ccluster),
                cstrt      =as.integer(cstrt), 
				iflag      =as.integer(iflag),
				num        =as.double(num),
				denom      =as.double(denom),
				fs         =as.double(fs),
				fm         =as.double(fm),
				muwork     =as.double(muwork),
				prob       =as.double(prob),
				seed       =as.integer(seed),
				sigmawork  =as.double(sigmawork),
				sigworkinv =as.double(sigworkinv),
				theta      =as.double(theta),
				workm1     =as.double(workm1),
				workm2     =as.double(workm2),
				workm3     =as.double(workm3),
				workmh1    =as.double(workmh1),
				workmh2    =as.double(workmh2),
				ywork      =as.double(ywork),
				iflagx     =as.integer(iflagx),
				workvx     =as.double(workvx),
				workmx     =as.double(workmx),
				worksam    =as.double(worksam),
				numcpo     =as.double(numcpo),
				denomcpo   =as.double(denomcpo),
				PACKAGE    ="DPpackage")

       #########################################################################################
       # save state
       #########################################################################################

         model.name <- "Bayesian Semiparametric Density Regression"		

		cpom<-matrix(foo$cpo,nrow=nrec,ncol=2)         
		cpo<-cpom[,1]         
		fso<-cpom[,2]
       
         varnames <- colnames(as.matrix(x))
         if(is.null(varnames))
         {
            varnames<-all.vars(cl)[2]
         }
         varnames <- c("y",varnames)
	
         state <- list(
                  alpha=foo$alpha,
                  m1=matrix(foo$m1,nrow=nvar,ncol=1),
                  muclus=matrix(foo$muclus,nrow=nrec+100,ncol=nvar),
                  ncluster=foo$ncluster,
                  psi1=matrix(foo$psi1,nrow=nvar,ncol=nvar),
                  k0=foo$k0,
                  sigmaclus=matrix(foo$sigmaclus,nrow=nrec+100,ncol=nuniq),
                  ss=foo$ss,
                  z=matrix(foo$z,nrow=nrec,ncol=nvar) 
                  )

         pnames1 <- NULL
         for(i in 1:nvar)
		 {
			pnames1 <- c(pnames1,paste("m1",varnames[i],sep=":"))
		 }
         pnames2 <- "k0"
       
         pnames3 <- NULL
		 for(i in 1:nvar)
		 {
			for(j in i:nvar)
			{
				if(i==j)
				{
					tmp<-varnames[i]
				}
				else
				{
					tmp <- paste(varnames[i],varnames[j],sep="-")
				}   
				pnames3 <- c(pnames3,paste("psi1",tmp,sep=":"))
			}	
		 }

         pnames4 <- c("ncluster","alpha")
         pnames <- c(pnames1,pnames2,pnames3,pnames4)

         thetasave <- matrix(foo$thetasave,nrow=nsave,ncol=nvar+nvar*(nvar+1)/2+3)
		 randsave <- matrix(foo$randsave,nrow=nsave,ncol=(nrec+2)*nvar+(nrec+1)*nvar*(nvar+1)/2)
		densp0 <- matrix(foo$densp0,nrow=nsave,ncol=npred)     
    densp.m <- matrix(foo$denspm,nrow=nsave,ncol=ngrid)
         densp.l <- matrix(foo$denspl,nrow=nsave,ncol=ngrid)

         meanfp.m <- matrix(foo$meanfpm,nrow=npred,ncol=ngrid)
		     meanfp.l <- matrix(foo$meanfpl,nrow=npred,ncol=ngrid)
		     meanfp.h <- matrix(foo$meanfph,nrow=npred,ncol=ngrid)
    
		meanfp.m1 <- matrix(foo$meanfpm1,nrow=npred,ncol=ngrid)
		meanfp.l1 <- matrix(foo$meanfpl1,nrow=npred,ncol=ngrid)
		meanfp.h1 <- matrix(foo$meanfph1,nrow=npred,ncol=ngrid)

         colnames(thetasave) <- pnames 

         coeff <- apply(thetasave,2,mean)
    

 
         save.state <- list(thetasave=thetasave,randsave=randsave,densp0=densp0,denspm=densp.m,denspl=densp.l,
                            meanfpm=meanfp.m,meanfpl=meanfp.l,meanfph=meanfp.h,meanfpm1=meanfp.m1,
                            meanfpl1=meanfp.l1,meanfph1=meanfp.h1,grid=grid,cpo=cpo,
                            fso=fso)
    
 
         z <- list(save.state=save.state)
                 
         cat("\n\n")
         class(z)<-"DPcdensity"
  	 return(z)
}


###                    
### Tools
###
### Copyright: Alejandro Jara, 2008
### Last modification: 25-05-2008.
###
