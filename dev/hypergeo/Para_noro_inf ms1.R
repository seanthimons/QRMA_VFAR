#automatic install of packages if they are not installed already
{
list.of.packages <- c(
  "foreach",
  "doParallel",
  "ranger",
  "extraDistr",
  "tidyverse",
  "VGAM",
  'future',
  'tictoc',
  'pbapply',
  'furrr',
  'tidyverse'
  
)

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages) > 0){
  install.packages(new.packages, dep=TRUE)
}

#loading packages
for(package.i in list.of.packages){
  suppressPackageStartupMessages(
    library(
      package.i, 
      character.only = TRUE
    )
  )
}
}
#Cores----
{
tic()
plan(multisession)
cores <- as.numeric(availableCores()) #grabs number of cores available to user
RNGkind("L'Ecuyer-CMRG")

#Total number of reps to be split across # of available cores
total_reps <- 5000

ftr1 <- vector('list',length = cores)

for (y in 1:cores) {
  names(ftr1)[y] <- paste0('t',y)
  }

ftr2 <- vector('list',length = cores)

for (y in 1:cores) {
  names(ftr2)[y] <- paste0('t',y)
  }
}
#loop 1-----
{
  tic()
for (z in 1:cores) {
  
  ftr1[[z]] <- future({
    setwd('~/marystemp/hgdr.new')
    dyn.load('hypg.dll')
    
    Probinfbpa <- function(a,b,d)
      .C("dr1f1",
         as.double(a), as.double(b), as.double(d),
         pinf = double(1))$pinf
    
    reps<- ceiling(total_reps/cores) #evenly splits loads across cores
    #transformed parameters for Noro dose-resposne
    ot=rbinorm(n=reps,mean1 = -.608, mean2 =.194, var1 = 1.79, var2 = 2.54, cov12 = -1.03)#GI SE+
    illIo=rbinorm(n=reps,mean1 = 1.74, mean2 =1.82, var1 = 5.55, var2 = 4.64, cov12 = -7.08E-1)
    #best estimate parameters for Noro dose-response
    Cj<-data.frame(p1=.393, p2=.767,p4=3.19,p3=0.801, min=0.001,mode=0.002,max=0.003)#N
    Daly.case<-rtriang(reps, a = Cj$min, b = Cj$max, c = Cj$mode)
    
    #Campy inputs
    #transparinf=ot=rbinorm(n=100000,mean1 = -.177, mean2 =.054, var1 = 1.303, var2 = 1.070, cov12 = -.041)#Challenge
    #transparill=illIo=rbinorm(n=100000,mean1 = -2.744, mean2 =-4.89E-3, var1 = 1.337, var2 = .993, cov12 = 0.01)#Challenge
    
test1<-c(8.8,8.9,9.0)       ################################################11/9change
    f1_results <- vector(mode = 'list', length = 3)
    names(f1_results) <- c('percent_0.1_und', 'target','percent_0.1_over')
    for(t in 1:3){
      #Irrigation
      tase<-0
      arisk8<-arisk7<-arisk6<-arisk5<-arisk4<-arisk2<-arisk3<-arisk<-rep(0,reps)
      V2<-.001
      LRT<-test1[t]
      
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
      f1_results[[t]] <- arisk4 
      }
   return(f1_results)
  })
}
  
loop1 <- lapply(ftr1, FUN = value) #Resolves the future

saveRDS(loop1, file = paste0('~/marystemp/loop1_noro__daly_inf',as.numeric(Sys.time()),'.RDS'))
toc()
rm(ftr1)
}


#loop 2----

{
  tic()

for (z in 1:cores) {
  
  ftr2[[z]] <- future({
    setwd('~/marystemp/hgdr.new')
    dyn.load('hypg.dll')
    
    Probinfbpa <- function(a,b,d)
      .C("dr1f1",
         as.double(a), as.double(b), as.double(d),
         pinf = double(1))$pinf
    
    reps<- ceiling(total_reps/cores)
    #transformed parameters for Noro dose-resposne
    ot=rbinorm(n=reps,mean1 = -.608, mean2 =.194, var1 = 1.79, var2 = 2.54, cov12 = -1.03)#GI SE+
    illIo=rbinorm(n=reps,mean1 = 1.74, mean2 =1.82, var1 = 5.55, var2 = 4.64, cov12 = -7.08E-1)
    #best estimate parameters for Noro dose-response
    Cj<-data.frame(p1=.393, p2=.767,p4=3.19,p3=0.801, min=0.001,mode=0.002,max=0.003)#N
    Daly.case<-rtriang(reps, a = Cj$min, b = Cj$max, c = Cj$mode)

test2<-c(9.0,9.1,9.2) #####################11/9 change
f2_results <- vector(mode = 'list', length = 3)
names(f2_results) <- c('percent_0.1_und', 'target','percent_0.1_over')
for(t in 1:3){
  #Domestic
  tase<-0
  arisk8<-arisk7<-arisk6<-arisk5<-arisk4<-arisk2<-arisk9<-arisk<-rep(0,reps)
  V2<-.00004
  V<-2
  LRT<-test2[t]
  
  for(ii in 1:reps){
    N<-365
    pathden<-rnorm(N,4.7,1.5)-1
    
    N<-1
    pathden1<-rnorm(N,4.7,1.5)-1
    
    
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
    tase2<-sapply(X=V*10^(-LRT+pathden1),FUN=noro)
    
    
    arisk4[ii]<-(1-prod(1-tase[2,]))#*Daly.case[i]
    arisk5[ii]<-(1-prod(1-tase2[2,]))#*Daly.case[i]          ############################11/9 change
    
    arisk6[ii]<-(.9*(arisk4[ii])+.1*(1-(1-arisk4[ii])*(1-arisk5[ii])))#*Daly.case[i]
    
    arisk7[ii]<-(1-prod(1-tase[1,]))#*Daly.case[i]
    arisk8[ii]<-(1-prod(1-tase[1,]))#*Daly.case[i]
    
    arisk9[ii]<-(.9*(arisk7[ii])+.1*(1-(1-arisk7[ii])*(1-arisk8[ii])))*Daly.case[ii]
  }
  f2_results[[t]] <- arisk6
  }
return(f2_results)
  })
}

loop2 <-lapply(ftr2, FUN = value) #Resolves the future

saveRDS(loop2, file = paste0('~/marystemp/loop2_noro__daly_inf',as.numeric(Sys.time()),'.RDS'))
toc()
rm(ftr2)
}

toc()

#results----

{
  
  #loop1
  
  loop1_under_1 <- map_dfr(loop1_1, 1) %>% unlist()
  loop1_target_1 <- map_dfr(loop1_1, 2) %>% unlist()
  loop1_over_1 <- map_dfr(loop1_1, 3) %>% unlist()
  
  loop1_under_2 <- map_dfr(loop1_2, 1) %>% unlist()
  loop1_target_2 <- map_dfr(loop1_2, 2) %>% unlist()
  loop1_over_2 <- map_dfr(loop1_2, 3) %>% unlist()
  
  loop1_under <- c(loop1_under_1, loop1_under_2)
  loop1_target <- c(loop1_target_1, loop1_target_2)
  loop1_over <- c(loop1_over_1, loop1_over_2)
  
#loop2

loop2_under_1 <- map_dfr(loop2_1, 1) %>% unlist()
loop2_target_1 <- map_dfr(loop2_1, 2) %>% unlist()
loop2_over_1 <- map_dfr(loop2_1, 3) %>% unlist()

loop2_under_2 <- map_dfr(loop2_2, 1) %>% unlist()
loop2_target_2 <- map_dfr(loop2_2, 2) %>% unlist()
loop2_over_2 <- map_dfr(loop2_2, 3) %>% unlist()

loop2_under <- c(loop2_under_1, loop2_under_2)
loop2_target <- c(loop2_target_1, loop2_target_2)
loop2_over <- c(loop2_over_1, loop2_over_2)

setwd('~/marystemp/')
fileConn <- file('noro_daly_inf.txt')
{cat('5%              50%             95%         99%','\n', file = "noro_daly_inf.txt", append = TRUE)
  cat('Loop 1\n',file = "noro_daly_inf.txt", append = TRUE)
  cat(test1,'\n',file = "noro_daly_inf.txt", append = TRUE)
  cat(quantile(loop1_under, c(0.05, .5, .95, .99)), '\n', file = "noro_daly_inf.txt", append = TRUE)
  cat(quantile(loop1_target, c(0.05, .5, .95, .99)),'\n', file = "noro_daly_inf.txt", append = TRUE)
  cat(quantile(loop1_over, c(0.05, .5, .95, .99)),'\n', file = "noro_daly_inf.txt", append = TRUE)
  cat('Loop 2\n',file = "noro_daly_inf.txt", append = TRUE)
  cat(test2, '\n',file = "noro_daly_inf.txt", append = TRUE)
  cat(quantile(loop2_under, c(0.05, .5, .95, .99)), '\n', file = "noro_daly_inf.txt", append = TRUE)
  cat(quantile(loop2_target, c(0.05, .5, .95, .99)),'\n', file = "noro_daly_inf.txt", append = TRUE)
  cat(quantile(loop2_over, c(0.05, .5, .95, .99)),'\n', file = "noro_daly_inf.txt", append = TRUE)
}
close(fileConn)
}




