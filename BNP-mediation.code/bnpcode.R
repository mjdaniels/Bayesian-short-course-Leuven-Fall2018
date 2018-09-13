install.packages("BNPMediation")
library(DPpackage)
library(BNPMediation)
stride<-read.csv("STRIDE_data.csv", sep=",", header=TRUE, na.strings=c("","#NULL!","NA"))
#outcomes, etc. for stride data
out<-log(stride$WKTOT12+1) # outcome
m<-stride$BEHPRM6 # mediator
x1<-stride$AGE
x2<-stride$BBMI

trt<-ifelse(stride$GROUP==3, 0, 1)

#data structure for the function 
w<-cbind(out,m,x1,x2)
w1<-w[which(trt==1),] # dataTreatment
w0<-w[which(trt==0),] # dataControl

n0<-length(which(complete.cases(w0))) # observed data
n1<-length(which(complete.cases(w1)))

RANGE<-function(x){
  (summary(x)[6]-summary(x)[1])^2/16
}

#------ Prior information
wbar <- apply(w,2,function(x) median(x,na.rm=TRUE))
wcov <- diag(c(RANGE(out), RANGE(m), RANGE(x1), RANGE(x2)),4,4)
prior <- list(a0=1,b0=1,nu1=16,nu2=16,s2=wcov,m2=wbar,psiinv2=2*solve(wcov),tau1=10,tau2=5)

#------ Initial State
state <- NULL

#------ MCMC parameters
mcmc <- list(nburn=1000,nsave=2000,nskip=4,ndisplay=500)

#------ Conditional Effects
#bnpconmediation(w1, w0, prior=prior, mcmc=mcmc, state=state, status = TRUE,na.action, q = dim(w)[2], NN = 10, n1 = n1, n0 = n0, extra.thin = 0,cond.values = c(2.5,2), col.values = c(1,2),seed = 12345)


#------ Marginal Effects
# type help(bnpconmediation) for details on the arguments
results = bnpmediation(w1, w0, prior=prior, mcmc=mcmc, state=state, status = TRUE,na.action, q = dim(w)[2], NN = 10, n1 = n1, n0 = n0, extra.thin = 0,seed = 12345)

# output > names(results)
# [1] "Y11"    "Y00"    "Y10"    "ENIE"   "ENDE"   "ETE"    "TE.c.i" "IE.c.i"
# [9] "DE.c.i" "call"  

