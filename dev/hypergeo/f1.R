setwd('~/marystemp/hgdr.new')
dyn.load('hypg.dll')

Probinfbpa <- function(a,b,d)
  .C("dr1f1",
     as.double(a), as.double(b), as.double(d),
     pinf = double(1))$pinf

reps<-10
#transformed parameters for Noro dose-resposne
ot=rbinorm(n=reps,mean1 = -.608, mean2 =.194, var1 = 1.79, var2 = 2.54, cov12 = -1.03)#GI SE+
illIo=rbinorm(n=reps,mean1 = 1.74, mean2 =1.82, var1 = 5.55, var2 = 4.64, cov12 = -7.08E-1)
#best estimate parameters for Noro dose-response
Cj<-data.frame(p1=.393, p2=.767,p4=3.19,p3=0.801, min=0.001,mode=0.002,max=0.003)#N
Daly.case<-rtriang(reps, a = Cj$min, b = Cj$max, c = Cj$mode)

#Campy inputs
#transparinf=ot=rbinorm(n=100000,mean1 = -.177, mean2 =.054, var1 = 1.303, var2 = 1.070, cov12 = -.041)#Challenge
#transparill=illIo=rbinorm(n=100000,mean1 = -2.744, mean2 =-4.89E-3, var1 = 1.337, var2 = .993, cov12 = 0.01)#Challenge

test<-c(8.9,9.0,9.1)
init1 <- vector(mode = 'list', length = 3)

for(t in 1:3){
  #Irrigation
  tase<-0
  arisk8<-arisk7<-arisk6<-arisk5<-arisk4<-arisk2<-arisk3<-arisk<-rep(0,reps)
  V2<-.001
  LRT<-test[t]
  
  for(ii in 1:reps){
    N<-50
    pathden<-rnorm(N,4.7,1.5)-1
    noro<-function(X){
      risk<-ill<-rep(0,reps)
      for (i in 1:reps){
        risk[i]=Probinfbpa(a=((exp(ot[i,1]+ot[i,2])/(1+exp(ot[i,1])))),b=exp(ot[i,2])* (1-(exp(ot[i,1])/(1+exp(ot[i,1])))),X)
        r<-((exp(illIo[i,1]+illIo[i,2])/(1+exp(illIo[i,1]))))
        n<-exp(illIo[i,2])* (1-(exp(illIo[i,1])/(1+exp(illIo[i,1]))))
        ill[i]=risk[i]*(1-(1+X/n)^-r)#se+
      }
      return(c(median(ill),median(risk)))
    }
    
    #noro<-function(X){Probinfbpa(a=Cj$p1,b=Cj$p2,V2*10^(-LRT+X))*(1-(1+(V2*10^(-LRT+X))/Cj$p3)^-Cj$p4)}
    tase<-sapply(X=V2*10^(-LRT+pathden),FUN=noro)
    
    arisk4[ii]<-(1-prod(1-tase[2,]))#*Daly.case[i]
    arisk5[ii]<-(1-prod(1-tase[1,]))*Daly.case[ii]
  }
  
  #print(quantile(arisk4,c(.05,.5,.95,.99)))
  init1[[t]] <- arisk4
}
  