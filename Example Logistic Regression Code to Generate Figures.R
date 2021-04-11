# Analysis of Cp Values ---------------------------------------------------

MetricSummary<-ddply(SDYdataMedsPooled[,c(KeepCols1,"Cp")],KeepCols1,colwise(mean,na.rm=TRUE))
ScalingFactor<-0.05

PosNegCols<-c("blue","red")
PosNegPch<-c(21,22)
SpecMatrix<-data.frame(LineWidths=c(1,0.5,0.5),TCKvals=c(-0.07,-0.06,-0.05),Xlabels=c("1X LoD\n","10X LoD\n","100X LoD\n"),
                       XlabelLines=c(1,0.65,0.5),XcexVals=c(0.9,0.8,0.8))

#Set up empty lists and perform analysis
AssayResults<-vector(mode="list",length=length(unique(SDYdataRepeated$AssayName)))
AssayResultTrim<-vector(mode="list",length=length(unique(SDYdataRepeated$AssayName)))
SDYdataRepeated.glm<-vector(mode="list",length=length(unique(SDYdataRepeated$AssayName)))
ConcPredict<-vector(mode="list",length=length(unique(SDYdataRepeated$AssayName)))
PredVals<-vector(mode="list",length=length(unique(SDYdataRepeated$AssayName)))
CpLM<-vector(mode="list",length=length(unique(SDYdataRepeated$AssayName)))
LODnew<-c()
AssayVals<-c()
for (assays in 1:length(unique(SDYdataRepeated$AssayName)))
{
  AssayVals<-which(SDYdataRepeated$AssayName==unique(SDYdataRepeated$AssayName)[assays])
  AssayResults[[assays]]<-data.frame(SDYdataRepeated[AssayVals,],log10(SDYdataRepeated[AssayVals,"Dilution"]))
  colnames(AssayResults[[assays]])<-c("Assay","Call","Conc","value","logConc")
  AssayResultTrim[[assays]]<-unique(AssayResults[[assays]])
  SDYdataRepeated.glm[[assays]]<-glm(Call~logConc,data=AssayResults[[assays]],family=binomial("logit"))
  ConcPredict[[assays]]<-data.frame(logConc=seq(min(AssayResults[[assays]][,"logConc"]),max(AssayResults[[assays]][,"logConc"]),length.out=1001))
  PredVals[[assays]]<-predict(SDYdataRepeated.glm[[assays]],newdata=data.frame(logConc=ConcPredict[[assays]]),type="response")
  LODnew[assays]<-10^unlist(ConcPredict[[assays]])[which.min(abs(PredVals[[assays]]-LODlim/100))]
  rm(AssayVals)
  
  TempCps<-SDYdataMedsPooled[which(SDYdataMedsPooled$AssayName==AssayResults[[assays]][1,"Assay"]),]
  ScaledMin<-(1-ScalingFactor)*min(TempCps$Cp,na.rm=TRUE)
  ScaledMax<-(1+ScalingFactor)*max(TempCps$Cp,na.rm=TRUE)
  TempCps$TempCpsNorm<-(TempCps$Cp-ScaledMin)/(ScaledMax-ScaledMin)
  #TempCutoff<-(Cutoffs$Cutoff[which(Cutoffs$Assay==unique(SDYdataRepeated$AssayName)[assays])]-ScaledMin)/(ScaledMax-ScaledMin)
  
  CpLM[[assays]]<-lm(TempCpsNorm~log10(Dilution),data=TempCps)
  LoDFactors<-c(log10(LODnew[assays]),log10(10*LODnew[assays]),log10(100*LODnew[assays]))
  LoDValuesCpNorm<-coef(CpLM[[assays]])[2]*LoDFactors+coef(CpLM[[assays]])[1]
  LoDValuesCp<-round(LoDValuesCpNorm*(ScaledMax-ScaledMin)+ScaledMin,1)
  png(paste(SaveDir,"DF-SDY-030217-Global Fever Raptor Sub-Mix LoD ",AssayResults[[assays]]$Assay[1]," Assay Logit with Cp Values.png",sep=""),
      pointsize=6.5,res=600,height=8.5,width=8.5*1.5,units="cm")
  layout(mat=matrix(c(1,2),1,2),widths=c(9,3))
  par(mar=c(5,4,1,8))
  plot(AssayResults[[assays]][,"logConc"],unlist(lapply(AssayResults[[assays]][,"Call"],as.character)),type="n",ylab="",xlab="",
       axes=FALSE,ylim=c(-0.05,1.05),panel.last={
         set.seed(60617)
         abline(h=c(0,1),lty=2,col="gray75",lwd=0.5)
         
         lines(data.frame(ConcPredict[[assays]],PredVals[[assays]]),col=Colors[2],lwd=2)
         abline(a=coef(CpLM[[assays]])[1],b=coef(CpLM[[assays]])[2],lty=2)
         abline(h=0.95)
         points(AssayResults[[assays]][,"logConc"]+rnorm(length(AssayResults[[assays]][,"logConc"]),0,0.03),
                unlist(lapply(AssayResults[[assays]][,"Call"],as.numeric))+rnorm(length(AssayResults[[assays]][,"logConc"]),0,0.007),
                lwd=0.4,cex=1)
         if(any(AssayResultTrim[[assays]]$Call==0)) {
           TempVals<-AssayResultTrim[[assays]][which(AssayResultTrim[[assays]]$Call==0),]
           text(TempVals$logConc,TempVals$Call-0.05,TempVals$value)
         }
         if(any(AssayResultTrim[[assays]]$Call==1)) {
           TempVals<-AssayResultTrim[[assays]][which(AssayResultTrim[[assays]]$Call==1),]
           text(TempVals$logConc,TempVals$Call+0.05,TempVals$value)
         }
         
         for (lodVals in 1:length(LoDFactors))
         {
           segments(LoDFactors[lodVals],-1,LoDFactors[lodVals],LoDValuesCpNorm[lodVals],lwd=SpecMatrix$LineWidths[lodVals],col=Colors[2],lty=3)
           segments(LoDFactors[lodVals],LoDValuesCpNorm[lodVals],100,LoDValuesCpNorm[lodVals],lwd=SpecMatrix$LineWidths[lodVals],col=Colors[2],lty=3)
           axis(1,at=LoDFactors[lodVals],labels=FALSE,las=2,lty=3,tck=SpecMatrix$TCKvals[lodVals],col.ticks=Colors[2],lwd=SpecMatrix$LineWidths[lodVals])
           axis(1,at=LoDFactors[lodVals],labels=paste(SpecMatrix$Xlabels[lodVals],signif(10^LoDFactors[lodVals],3)),las=2,cex.axis=SpecMatrix$XcexVals[lodVals],
                line=SpecMatrix$XlabelLines[lodVals],tick=FALSE,col.axis=Colors[2])
           axis(4,at=LoDValuesCpNorm[lodVals],las=2,col=Colors[2],col.axis=Colors[2],lwd=SpecMatrix$LineWidths[lodVals],tck=-0.1,lty=3,labels=FALSE)
           axis(4,at=LoDValuesCpNorm[lodVals],las=2,col=Colors[2],col.axis=Colors[2],line=1.8,tick=FALSE,cex.axis=0.8,
                lwd=SpecMatrix$LineWidths[lodVals],labels=paste(SpecMatrix$Xlabels[lodVals],"Cp ~ ",LoDValuesCp[lodVals]))
         }
         segments(LoDFactors[1],-1,LoDFactors[1],0.95,lwd=SpecMatrix$LineWidths[1],col=Colors[2],lty=3)
         
         axis(1,at=seq(-10,10,1),labels=signif(10^(seq(-10,10,1)),3),las=2,tck=-0.025,line=-0.75)
         axis(1,at=log10(TickVals),labels=FALSE,las=2,tck=-0.01,line=-0.75)
         
         axis(2,las=2,at=seq(0,1,0.25))
         axis(2,las=2,at=0.95,cex.axis=0.9,font=2)
         axis(4,at=seq(0,0.9,length.out=5),las=2,col=PosNegCols[1],col.axis=PosNegCols[1],tick=FALSE,line=-0.2,
              labels=round(seq(min(TempCps$Cp,na.rm=TRUE),max(TempCps$Cp,na.rm=TRUE),length.out=5),1))
         axis(4,at=seq(0,0.9,length.out=5),las=2,col=PosNegCols[1],col.axis=PosNegCols[1],labels=FALSE)
         title(xlab="Dilution",line=4)
         title(ylab="Probability")
         #
         mtext("Cp Values",4,line=7,col=PosNegCols[1])
         for (result in c(0,1))
           points(log10(TempCps$Dilution[which(TempCps$Call==result)])+
                    rnorm(length(TempCps$Dilution[which(TempCps$Call==result)]),0,0.015),
                  TempCps$TempCpsNorm[which(TempCps$Call==result)],
                  bg=PosNegCols[c(0,1)==result],pch=PosNegPch[c(0,1)==result],lwd=0.5,cex=1.2)
       })
  par(mar=c(1,0,3,0))
  plot.new()
  legend("left",c("Jittered Data","Logit Best Fit Line","LoD95 Values","Cp Best Fit Line"),title="Logit Curve",
         lty=c(NA,1,3,2),col=c("black",Colors[c(2,2)],"black"),
         pch=c(1,rep(NA,3)),bty="n",lwd=c(NA,2,0.8,1))
  legend(0,0.44,legend="Numbers represent\nthe number of\nobservations of\nPositive (1) or\nNegative (0)",bty="n")
  legend("topleft",c("Negative Assay Call","Positive Assay Call","Points were jittered"),title="Cp Points",
         pch=c(PosNegPch,NA),pt.bg=c(PosNegCols,NA),pt.lwd=0.5,pt.cex=1.2,bty="n")
  dev.off()
}
