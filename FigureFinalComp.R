library(ggplot2)

#Model
model <- function(t, state, parms){
  state <- ifelse(state < 0, 0, state)
  with(as.list(c(state,parms)),{
    #parameter re-definitions
    ACT <- ifelse(t<tau[1],1,ifelse(((t>tau[2])&(t<tau[3])),1,0))
    lMA <- ACT*lA; lEA <- lMA; r <- ACT*r; dC <- ACT*dC
    alpha <- 1/(1+(seq(n)/cutof)^slop)
    treset <- ifelse(t<86,t,t-86); delMQ <- ACT*delMQ*exp(-(aq*treset)^5)
    dA <- (1-ACT)*dA; dMA <- dA; dEA <- dA; ae <- ifelse(t<86,ae,2*ae)
    delEA <- ACT*delEA*(1-exp(-(ae*treset)^5)); dMC <- lMC; dEC <- lEC
    muA <- ifelse(t<86,dA/19,dA/4); muMA <- muA; muEA <- muA; muMQ <- (1-ACT)*10
    
    #variable name definitions
    MAG <- state[1:n]; MAGb <- c(0,MAG[1:(n-1)])
    MAR <- state[(n+1):(2*n)]; MARb <- c(0,MAR[1:(n-1)])
    EAG <- state[(2*n+1):(3*n)]; EAGb <- c(0,EAG[1:(n-1)])
    EAR <- state[(3*n+1):(4*n)]; EARb <- c(0,EAR[1:(n-1)])
    MCG <- state[(4*n+1):(5*n)]; MCGb <- c(0,MCG[1:(n-1)])
    MCR <- state[(5*n+1):(6*n)]; MCRb <- c(0,MCR[1:(n-1)])
    ECG <- state[(6*n+1):(7*n)]; ECGb <- c(0,ECG[1:(n-1)])
    ECR <- state[(7*n+1):(8*n)]; ECRb <- c(0,ECR[1:(n-1)])
    MQnG <- state[(8*n+1):(9*n)]
    MQnR <- state[(9*n+1):(10*n)]
    MQG <- state[(10*n+1):(11*n)]
    MQR <- state[(11*n+1):(12*n)]
    
    #system of equations
    dMAG <- rep(0,n); dMAR <- rep(0,n); dEAG <- rep(0,n); dEAR <- rep(0,n)
    dMCG <- rep(0,n); dMCR <- rep(0,n); dECG <- rep(0,n); dECR <- rep(0,n)
    dMQnG <- rep(0,n); dMQnR <- rep(0,n); dMQG <- rep(0,n); dMQR <- rep(0,n)
    if(t>ton){
      dMAG <- r*alpha*(MCG+MQG) + (2-sp)*lMA*MAGb -
        (lMA+delEA+delMQ+muMA+dMA)*MAG
      dMAR <- r*alpha*(MCR+MQR) + sp*lMA*MAGb + 2*lMA*MARb -
        (lMA+delEA+delMQ+muMA+dMA)*MAR
      dEAG <- delEA*MAG + (2-sp)*lEA*EAGb -
        (lEA+muEA+dEA)*EAG
      dEAR <- delEA*MAR + sp*lEA*EAGb + 2*lEA*EARb -
        (lEA+muEA+dEA)*EAR
      dMCG <- muMA*MAG + (2-sp)*lMC*MCGb -
        (lMC+dMC+r*alpha+dC*(1-alpha))*MCG
      dMCR <- muMA*MAR + sp*lMC*MCGb + 2*lMC*MCRb -
        (lMC+dMC+r*alpha+dC*(1-alpha))*MCR
      dECG <- muEA*EAG + (2-sp)*lEC*ECGb - (lEC+dEC+dC)*ECG
      dECR <- muEA*EAR + sp*lEC*ECGb + 2*lEC*ECRb - (lEC+dEC+dC)*ECR
      dMQnG <- delMQ*MAG - muMQ*MQnG
      dMQnR <- delMQ*MAR - muMQ*MQnR
      dMQG <- muMQ*MQnG - r*alpha*MQG
      dMQR <- muMQ*MQnR - r*alpha*MQR
    }
    return(list(c(dMAG,dMAR,dEAG,dEAR,dMCG,dMCR,
                  dECG,dECR,dMQnG,dMQnR,dMQG,dMQR)))
  })
}

#Compute percentages in different populations 
tweak <- function(nsol,parms) {
  with(as.list(parms),{
    effs <- c(seq((2*n+1),(4*n)),seq((6*n+1),(8*n)))+1
    mems <- setdiff(seq(2,ncol(nsol)),effs)
    quis <- seq((8*n+1),(12*n))+1
    rquis <- c(seq((9*n+1),(10*n)),seq((11*n+1),(12*n)))+1
    reffs <- c(seq((3*n+1),(4*n)),seq((7*n+1),(8*n)))+1
    rmems <- c(seq((n+1),(2*n)),seq((5*n+1),(6*n)),
               seq((9*n+1),(10*n)),seq((11*n+1),(12*n)))+1
    
    Tmult <- apply(nsol[,mems],1,sum); RTmult <- apply(nsol[,rmems],1,sum)
    QTmult <- apply(nsol[,quis],1,sum); RQTmult <- apply(nsol[,rquis],1,sum)
    Tterm <- apply(nsol[,effs],1,sum); RTterm <- apply(nsol[,reffs],1,sum)
    
    nsol$Tmult <- log10(Tmult); nsol$FTmult <- 100*(RTmult)/Tmult
    nsol$Tterm <- log10(Tterm); nsol$FTterm <- 100*(RTterm)/Tterm
    nsol$QTmult <- log10(QTmult); nsol$FQTmult <- 100*(RQTmult)/QTmult
    
    nsol$RTmult <- log10(RTmult); nsol$RTterm <- log10(RTterm)
    nsol$RQTmult <- log10(RQTmult)
    
    nsol$QTmult[nsol$time==0] <- 0; nsol$FQTmult[nsol$time==0] <- 0
    nsol$RQTmult[nsol$time==0] <- 0
    nsol$Tterm[nsol$time==0] <- 0; nsol$FTterm[nsol$time==0] <- 0
    nsol$RTterm[nsol$time==0] <- 0
    nsol$RTmult[nsol$time==0] <- 0
    
    nsol$time <- nsol$time - 86
    return(nsol)
  })
}

#Plot the distribution of the populations at time tt
dist <- function(tt){
  x <- nsolS
  mAG <- x[which(x$time==tt),1+seq(n)]
  mAR <- x[which(x$time==tt),1+n+seq(n)]
  mCG <- x[which(x$time==tt),1+4*n+seq(n)]
  mCR <- x[which(x$time==tt),1+5*n+seq(n)]
  mQnG <- x[which(x$time==tt),1+8*n+seq(n)]
  mQnR <- x[which(x$time==tt),1+9*n+seq(n)]
  mQG <- x[which(x$time==tt),1+10*n+seq(n)]
  mQR <- x[which(x$time==tt),1+11*n+seq(n)]
  mQGR <- mQnG+mQnR+mQG+mQR; muQ <- sum(mQGR*seq(n))/sum(mQGR)
  mCGR <- mAG+mAR+mCG+mCR; muC <- sum(mCGR*seq(n))/sum(mCGR)
  mGR <- mCGR+mQGR; muM <- sum(mGR*seq(n))/sum(mGR)
  eAG <- x[which(x$time==tt),1+2*n+seq(n)]
  eAR <- x[which(x$time==tt),1+3*n+seq(n)]
  eCG <- x[which(x$time==tt),1+6*n+seq(n)]
  eCR <- x[which(x$time==tt),1+7*n+seq(n)]
  eGR <- eAG+eAR+eCG+eCR; muE <- sum(eGR*seq(n))/sum(eGR)
  alph <- 1/(1+(seq(n)/p['cutof'])^p['slop'])
  alph <- alph*max(c(unlist(unname(eGR)),unlist(unname(mGR))))
  print(paste('Qmult/Cmult: ',sum(mQGR)/sum(mCGR),sep=""))
  
  if(cmnd != 'prnt'){
    Alph <<- c(Alph,alph)
  }
  else{
    Alph <<- c(Alph,alph)
    dataDist <- data.frame(div=rep(seq(n),5),data=c(unlist(unname(eGR)),
                                                    unlist(unname(mGR)),Alph),
                           col=rep(c('T_TERM','T_MULT',
                                     'func1','func2','func3'),each=n))
    print(ggplot(dataDist,aes(x=div,y=data,col=col))+
            geom_line(size=1)+theme_test()+labs(x='Generation number')+
            geom_vline(xintercept=muQ,linetype='dashed',color='green')+
            geom_vline(xintercept=muM,linetype='dashed',color='black')+
            geom_vline(xintercept=muC,linetype='dashed',color='blue')+
            geom_vline(xintercept=muE,linetype='dashed',color='red')+
            scale_y_continuous('Cell numbers',
            sec.axis=sec_axis(~./max(dataDist$data),
                              name='Fraction re-activated')))
  }
}

#Read data
dataB <- read.table('Blood.txt',header=T)
dataB$Day <- dataB$Day - 86
TNaiveB <- mean(dataB$Naive)
dataB$Tmult <- dataB$GFPpCD27pKLRG1n*TNaiveB/dataB$Naive
dataB$Tterm <- dataB$GFPpCD27nKLRG1p*TNaiveB/dataB$Naive
dataB$FTmult <- dataB$tdTomatopCD27pKLRG1n/dataB$GFPpCD27pKLRG1n
dataB$FTterm <- dataB$tdTomatopCD27nKLRG1p/dataB$GFPpCD27nKLRG1p
dataB <- data.frame(time=dataB$Day,Tmult=dataB$Tmult,FTmult=dataB$FTmult,
                    Tterm=dataB$Tterm,FTterm=dataB$FTterm)

TNaiveS <- 2*10^7
dataS <- read.table('Spleen.txt',header=T)
dataS$Day <- dataS$Day - 86
dataS$Tmult <- dataS$GFPpCD27pKLRG1n*TNaiveS/dataS$Naive
dataS$Tterm <- dataS$GFPpCD27nKLRG1p*TNaiveS/dataS$Naive
dataS$FTmult <- dataS$tdTomatopCD27pKLRG1n/dataS$GFPpCD27pKLRG1n
dataS$FTterm <- dataS$tdTomatopCD27nKLRG1p/dataS$GFPpCD27nKLRG1p
dataS <- data.frame(time=dataS$Day,Tmult=dataS$Tmult,FTmult=dataS$FTmult,
                    Tterm=dataS$Tterm,FTterm=dataS$FTterm)

dataB$Tmult <- log10(dataB$Tmult); dataB$Tterm <- log10(dataB$Tterm)
dataS$Tmult <- log10(dataS$Tmult); dataS$Tterm <- log10(dataS$Tterm)
dataB$FTmult <- 100*dataB$FTmult; dataB$FTterm <- 100*dataB$FTterm
dataS$FTmult <- 100*dataS$FTmult; dataS$FTterm <- 100*dataS$FTterm

#initial value definitions
n <- 100; tend <- 111; tau <- c(6,86,90) #n <- 300; tend <- 200
MAG <- rep(0,n); names(MAG) <- paste("MAG",seq(1,n),sep="")
MAR <- rep(0,n); names(MAR) <- paste("MAR",seq(1,n),sep="")
EAG <- rep(0,n); names(EAG) <- paste("EAG",seq(1,n),sep="")
EAR <- rep(0,n); names(EAR) <- paste("EAR",seq(1,n),sep="")
MCG <- rep(0,n); names(MCG) <- paste("MCG",seq(1,n),sep="")
MCR <- rep(0,n); names(MCR) <- paste("MCR",seq(1,n),sep="")
ECG <- rep(0,n); names(ECG) <- paste("ECG",seq(1,n),sep="")
ECR <- rep(0,n); names(ECR) <- paste("ECR",seq(1,n),sep="")
MQnG <- rep(0,n); names(MQnG) <- paste("MQnG",seq(1,n),sep="")
MQnR <- rep(0,n); names(MQnR) <- paste("MQnR",seq(1,n),sep="")
MQG <- rep(0,n); names(MQG) <- paste("MQG",seq(1,n),sep="")
MQR <- rep(0,n); names(MQR) <- paste("MQR",seq(1,n),sep="")
sB <- c(MAG,MAR,EAG,EAR,MCG,MCR,ECG,ECR,MQnG,MQnR,MQG,MQR); sS <- sB

#Simulate to make desired data for final plots
MakeData <- function(initcell){
  sS['MAG1'] <- initcell[2]
  
  nsolS <<- run(tend,tweak="nsol<-tweak(nsol,parms)",parms=p,
                state=sS,arrest=c(3,"ton",tau),table=T,timeplot=F)
  
  tempM <- data.frame(time=nsolS$time,data=nsolS$FTmult,pop=popn)
  PlotM <<- rbind(PlotM,tempM)
  tempT <- data.frame(time=nsolS$time,data=nsolS$FTterm,pop=popn)
  PlotT <<- rbind(PlotT,tempT)
  
  dist(tm)
}

PlotM <- data.frame(); PlotT <- data.frame()

pdf('Figure4.pdf',useDingbats=F)
#parameter definitions
tm <- 86-86; cmnd <- 'dntprnt'; Alph <- c()
#O(3) quiescent cells
p <- c(lA=0.88,ae=0.15,delEA=2,aq=0.25,delMQ=0.01,dA=0.3,lMC=0.15,lEC=0.04,
       r=1,dC=0.3,cutof=10,slop=30,sp=0.0053,ton=1)
popn <- 'O(3):func1'; MakeData(c(0.01335458,1.49088696)*1e+4)

p <- c(lA=0.88,ae=0.15,delEA=2,aq=0.25,delMQ=0.01,dA=0.3,lMC=0.15,lEC=0.04,
       r=1,dC=0.3,cutof=10,slop=1,sp=0.0053,ton=1)
popn <- 'O(3):func2'; MakeData(c(0.01335458,1.49088696)*1e+4)

cmnd <- 'prnt'
p <- c(lA=0.88,ae=0.15,delEA=2,aq=0.25,delMQ=0.01,dA=0.3,lMC=0.15,lEC=0.04,
       r=1,dC=0.3,cutof=25,slop=30,sp=0.0053,ton=1)
popn <- 'O(3):func3'; MakeData(c(0.01335458,1.49088696)*1e+4)

tm <- 86-86; cmnd <- 'dntprnt'; Alph <- c()
#O(4) quiescent cells
p <- c(lA=0.89,ae=0.15,delEA=2,aq=0.25,delMQ=0.1,dA=0.3,lMC=0.30,lEC=0.04,
       r=1,dC=0.3,cutof=10,slop=30,sp=0.0053,ton=1)
popn <- 'O(4):func1'; MakeData(c(0.0149298,1.6748715)*1e+4)

p <- c(lA=0.89,ae=0.15,delEA=2,aq=0.25,delMQ=0.1,dA=0.3,lMC=0.30,lEC=0.04,
       r=1,dC=0.3,cutof=10,slop=1,sp=0.0053,ton=1)
popn <- 'O(4):func2'; MakeData(c(0.0149298,1.6748715)*1e+4)

cmnd <- 'prnt'
p <- c(lA=0.89,ae=0.15,delEA=2,aq=0.25,delMQ=0.1,dA=0.3,lMC=0.30,lEC=0.04,
       r=1,dC=0.3,cutof=25,slop=30,sp=0.0053,ton=1)
popn <- 'O(4):func3'; MakeData(c(0.0149298,1.6748715)*1e+4)
dev.off()

pdf('Figure5.pdf',useDingbats=F)
par(mar=c(2.6,2.6,1.6,0.2),mgp=c(1.5,0.5,0),mfrow=c(1,1))
dataFM <- data.frame(time=dataS$time,data=dataS$FTmult)
print(ggplot(dataFM,aes(x=time,y=data))+geom_point()+
        theme_classic()+geom_line(data=PlotM,aes(x=time,y=data,col=pop),size=1)+
        theme(plot.title=element_text(hjust=0.5),legend.text.align=0.3,
              strip.background=element_blank())+ggtitle('A')+
        labs(x='Time(in days)',y='Fraction labelled')+ylim(0,10)+xlim(0,25))

dataFT <- data.frame(time=dataS$time,data=dataS$FTterm)
print(ggplot(dataFT,aes(x=time,y=data))+geom_point()+
        theme_classic()+geom_line(data=PlotT,aes(x=time,y=data,col=pop),size=1)+
        theme(plot.title=element_text(hjust=0.5),legend.text.align=0.3,
              strip.background=element_blank())+ggtitle('B')+
        labs(x='Time(in days)',y='Fraction labelled')+ylim(0,10)+xlim(0,25))
dev.off()
