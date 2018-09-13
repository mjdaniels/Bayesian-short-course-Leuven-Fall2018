library(MASS)
library(Rlab)
library(msm)
library(splines)
library(BayesTree)
library(Hmisc)  
library(gtools)
devtools::install_github("mjdaniels/SequentialBART",ref="travis")
library(sbart)
library(DPpackage)

## change the directory
dat = read.csv("D:\\Downloads\\EHR_data.csv", header = TRUE,na.strings="NA")
y<-dat$scr_after-dat$scr
t<-dat$exposure


race1<-ifelse(dat$race==1,1,0)
race2<-ifelse(dat$race==2,1,0)

x1<-cbind(dat$scr,dat$age_n,dat$CrCl_CG,dat$vanco_daily_dose,dat$SBP,dat$MAP,
          dat$bun_scr_ratio,dat$albumin_cts,dat$TEMP_C)
x2<-cbind(dat$current_sepsis,dat$current_HF,dat$current_pneumonia,dat$current_dm,
          dat$current_ami,dat$current_anemia,dat$current_HTN,dat$current_liver,
          dat$current_CKD,dat$CURRENT_SI,dat$current_bacteremia,dat$current_uti,
          dat$current_arth_inf,dat$contrast,dat$cardiac,
          dat$diuretics,dat$ACEIsARBs,dat$cephalsporins,dat$others,
          dat$albumin, dat$nsaids_baseline,dat$LOCATION,race1,race2,dat$male)



##########  IMPUTING missing data using sequential bart
xx<-cbind(x1,x2)
yy<-t
x.type=c(rep(0,9),rep(1,25))
y.type=4

seqBART_r<-sbart::seqBART(x=xx,y=yy,x.type=x.type,y.type=y.type,prediction=TRUE)
pp_all<-seqBART_r$prediction[c(200,400,600,800,1000),]

pdft2<-NULL
cdft2<-NULL
pdfc2<-NULL
cdfc2<-NULL
nsave<-200
nskip<-1
nburn<-500
yt<-y[t==1]
yc<-y[t==0]
ygrid0<-seq(-0.3,0.8,length=100)
ngrid<-length(ygrid0)
for (ii in 1:5) 
  {
pp<-pp_all[ii,]
xt<-pp[t==1]
xc<-pp[t==0]
xxpred<-pp


## BART-DPM
#treated
mcmc <- list(nburn=nburn,nsave=nsave,nskip=nskip,ndisplay=100)
w<-cbind(yt,xt)
m2<-(apply(w,2,max)+apply(w,2,min))/2
s2<-diag(((apply(w,2,max)-apply(w,2,min))/4)^2)

prior <- list(a0=10,b0=1,nu1=4,nu2=4,m2=m2,s2=s2,
              psiinv2=4*solve(s2),tau1=2.01,tau2=2.01)
fit1 <- DPcdensity(yt,xt,xpred=xxpred,ngrid=ngrid,
                   grid=ygrid0,prior=prior,
                   mcmc=mcmc,state=NULL,status=TRUE)

pdft2<-rbind(pdft2,fit1$save.state$denspm)
cdft2<-rbind(cdft2,fit1$save.state$denspl)

#control
w<-cbind(yc,xc)
m2<-(apply(w,2,max)+apply(w,2,min))/2
s2<-diag(((apply(w,2,max)-apply(w,2,min))/4)^2)

prior <- list(a0=10,b0=1,nu1=4,nu2=4,m2=m2,s2=s2,
              psiinv2=4*solve(s2),tau1=2.01,tau2=2.01)
fit2 <- DPcdensity(yc,xc,xpred=xxpred,ngrid=ngrid,grid=ygrid0,
                   prior=prior,
                   mcmc=mcmc,state=NULL,status=TRUE)

pdfc2<-rbind(pdfc2,fit2$save.state$denspm)
cdfc2<-rbind(cdfc2,fit2$save.state$denspl)
}


quantile.fun<-function(cdf,aa) {
  nq<-length(aa)
  nl<-nrow(cdf)
  bb<-matrix(NA,nq,nl)
  for (l in 1:nq) {
    for (j in 1:nl) {
      if(aa[l]>cdf[j,ngrid]) {bb[l,j]=ygrid0[ngrid]}
      else {
        if(aa[l]<=cdf[j,1]) {bb[l,j]=ygrid0[1]}
        else {
          tar0<-cdf[j,]-aa[l]
          tar1<-tar0
          tar1[tar0<0]<-1
          i<-which.min(tar1)-1 
          bb[l,j]<-ygrid0[i]+(ygrid0[i+1]-ygrid0[i])*(aa[l]-cdf[j,i])/(cdf[j,i+1]-cdf[j,i])
        }}
    }}
  return(bb)
}

quan2<-quantile.fun(1-cdft2,c(0.1,0.25,0.5,0.75,0.9))-quantile.fun(1-cdfc2,c(0.1,0.25,0.5,0.75,0.9))
quantile_bnp<-cbind(rowMeans(quan2),t(apply(quan2,1,function (x) quantile(x,c(0.025,0.975,0.05,0.95,0.1,0.9)))))

