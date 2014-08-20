## Australian SPF simulation code. - Parallel computing version
## This code cycles through all species, then through Fs and produces 
## plots of depletion vs. F for all the stocks.
## Code written by F. Hurtado-Ferro. Please reference the author if you use this code.

##############################################################################################
# Libraries needed

require(snow)
require(snowfall)
require(RColorBrewer)

##############################################################################################
# Control

# Switches
run.again      <- F # Run the operating model?
read.again.R   <- F # Read the output using R?
read.again.for <- T # Read the output using fortran?
plot.how       <- 0 # Plot? 0 = No, 1 = Together, 2 = Sp by Sp

# Parallel stuff
Ncpus          <- 7 # How many CPUs to use for parallel computing

# Plotting options
# What to plot
NoUncPlot      <- F # Plot the figures without uncertainty?
B1pPlot        <- F # Plot the 1+ Biomass profiles?
# Colors
ColPal         <- brewer.pal(5,"RdYlGn")
EColLine       <- ColPal[1]
EColPol        <- ColPal[2]
WColLine       <- ColPal[5]  
WColPol        <- ColPal[4]
JackCol        <- ColPal[5]
# Window dimensions
Wheight        <- 4.5
Wwidth         <- 5.4

# File names and important paths
Outpath <- "C:/Users/Felipe/Desktop/Chile_MSE/FvD"
OMpath <- "C:/Users/Felipe/Desktop/Chile_MSE/V1.1"

DATfile <- "Chile_OM.DAT"
#Stocks <- c("Ered","Wred","Jack","Eblu","Wblu","Esar","Wsar")
#SppN <- c(1,1,2,3,3,4,4)
#StockN <- c(0,1,0,0,1,0,1)
dir.create(Outpath, showWarnings = FALSE)

# Limits for the Fs
MaxF <- 1.0
Fstep <- 0.025
FFs <- seq(0,MaxF,by=Fstep)
#FFs <- c(seq(0,0.25,by=0.01),seq(0.3,MaxF,by=Fstep))

##############################################################################################
# The big loop for running models

if(run.again==T){
  
  # A function to parallelize over Fs
  RunOM <- function(FFs,Outpath,OMpath,DATfiles){
    # Create the folders needed
    Fpath <- paste0(Outpath,"/F_",FFs)
    dir.create(Fpath)
    setwd(Fpath)
    # Copy the appropriate files
    file.copy(from=paste(OMpath,DATfile,sep="/"),to=paste(Fpath,DATfile,sep="/"),overwrite=T)
    file.copy(from=paste(OMpath,"SPFOM.SPEC",sep="/"),to=paste(Fpath,"SPFOM.SPEC",sep="/"),overwrite=T)
    file.copy(from=paste(OMpath,"SPFspp.SPEC",sep="/"),to=paste(Fpath,"SPFspp.SPEC",sep="/"),overwrite=T)
    # Modify the DAT file
    t.datfile <- readLines(DATfile)
    # Modify the SPEC file
    t.specfile <- readLines("SPFOM.SPEC")
    t.specfile[6] <- paste(FFs,0.7,sep="\t")
    writeLines(t.specfile,"SPFOM.SPEC")
    # Modify the spp file
#     t.sppfile <- readLines("SPFspp.SPEC")
#     t.sppfile[2] <- SppN[stock]
#     t.sppfile[4] <- StockN[stock]
#     writeLines(t.sppfile,"SPFspp.SPEC")
    # Run the model 
    shell(paste0(OMpath,"/SPFOM.exe"),intern=FALSE, wait=TRUE)
  }
  
  
  # Loop over stocks and run the function above
#  for(stock in 1:length(Stocks)){
    #stockpath <- paste(Outpath,Stocks[stock],sep="/")
    #dir.create(stockpath)
    
    #Parallel part
    sfInit(parallel=TRUE, cpus=Ncpus, type="SOCK")
    sfSapply(FFs,RunOM,Outpath,OMpath,DATfiles,SppN,StockN,stock)
    sfStop()   
    
#    cat("Finished ",Stocks[stock],"\n")
#  } #Stock loop close
}



##############################################################################################
# Read the output using R
if(read.again.R==T){
  # Create a table to store the results
  final.table <- data.frame(Performance_Measure=c("Mean_B1","Sigma_B1","Q05_B1","Q25_B1",
                                                  "Q50_B1","Q75_B1","Q95_B1",
                                                  "Mean_B15","Sigma_B15","Q05_B15","Q25_B15",
                                                  "Q50_B15","Q75_B15","Q95_B15",
                                                  "Mean_SSB","Sigma_SSB","Q05_SSB","Q25_SSB",
                                                  "Q50_SSB","Q75_SSB","Q95_SSB",
                                                  "Mean_SSB5","Sigma_SSB5","Q05_SSB5","Q25_SSB5",
                                                  "Q50_SSB5","Q75_SSB5","Q95_SSB5",
                                                  "Mean_C","Sigma_C","Q05_C","Q25_C",
                                                  "Q50_C","Q75_C","Q95_C",
                                                  "Mean_C5","Sigma_C5","Q05_C5","Q25_C5",
                                                  "Q50_C5","Q75_C5","Q95_C5"))
  
  # Populate the table
  for(stock in 1:length(Stocks)){
    stockpath <- paste(Outpath,Stocks[stock],sep="/")
    # Loop over Fs
    for(ff in 1:length(FFs)){
      # Get the right folder
      Fpath <- paste0(stockpath,"/F_",FFs[ff])
      setwd(Fpath)
      # Read the output
      outs <- read.table("SUMMARY.OUT")
      tt <- length(outs[,1])
      sims <- (length(outs[1,])-1)/9
      totts <- (tt-1)*sims
      
      #create a new matrix with the right result format (all results)
      outsnew <- as.matrix(outs[-tt,2:10])
      for(it in 2:sims)
        outsnew <- rbind(outsnew,as.matrix(outs[-tt,((it-1)*9+2):(it*9+1)]))
      #create a new matrix with the right result format (last 5 years)
      outs5 <- as.matrix(outs[(tt-5):(tt-1),2:10])
      for(it in 2:sims)
        outs5 <- rbind(outs5,as.matrix(outs[(tt-5):(tt-1),((it-1)*9+2):(it*9+1)]))
      
      #Calculate means and all that (all years)
      means <- apply(outsnew,2,mean,na.rm=T)
      #medians <- apply(outsnew,2,median,na.rm=T)
      SDs <- apply(outsnew,2,sd,na.rm=T)
      Quants <- apply(outsnew,2,quantile,na.rm=T,probs=c(0.05,0.25,0.5,0.75,0.95))
      #Calculate means and all that (last 5 years)
      means5 <- apply(outs5,2,mean,na.rm=T)
      #medians5 <- apply(outs5,2,median,na.rm=T)
      SDs5 <- apply(outs5,2,sd,na.rm=T)
      Quants5 <- apply(outs5,2,quantile,na.rm=T,probs=c(0.05,0.25,0.5,0.75,0.95))
      
      #Save to the final table
      final.table <- cbind(final.table,c(means[1],SDs[1],Quants[,1],
                                         means5[1],SDs5[1],Quants5[,1],
                                         means[2],SDs[2],Quants[,2],
                                         means5[2],SDs5[2],Quants5[,2],
                                         means[4],SDs[4],Quants[,4],
                                         means5[4],SDs5[4],Quants5[,4]))
      cat(Stocks[stock],"F =",FFs[ff],"\n")
    }
  }
  
  # Save the results
  dir.create(paste0(Outpath,"/AAA-Summary"))
  setwd(paste0(Outpath,"/AAA-Summary"))
  write.csv(final.table,"final.csv", row.names=F)
}

##############################################################################################
# Read the output using fortran
if(read.again.for==T){
  # Create a table to store the results
  plot.table <- data.frame(Performance_Measure=c("Mean_B1","Sigma_B1","Q05_B1","Q25_B1",
                                                 "Q50_B1","Q75_B1","Q95_B1",
                                                 "Mean_B15","Sigma_B15","Q05_B15","Q25_B15",
                                                 "Q50_B15","Q75_B15","Q95_B15",
                                                 "Mean_SSB","Sigma_SSB","Q05_SSB","Q25_SSB",
                                                 "Q50_SSB","Q75_SSB","Q95_SSB",
                                                 "Mean_SSB5","Sigma_SSB5","Q05_SSB5","Q25_SSB5",
                                                 "Q50_SSB5","Q75_SSB5","Q95_SSB5",
                                                 "Mean_C","Sigma_C","Q05_C","Q25_C",
                                                 "Q50_C","Q75_C","Q95_C",
                                                 "Mean_C5","Sigma_C5","Q05_C5","Q25_C5",
                                                 "Q50_C5","Q75_C5","Q95_C5"))
  
  report.table <- data.frame(Performance_Measure=c("Mean_B1","Sigma_B1","Mean_B15","Sigma_B15",
                                                   "Mean_SSB","Sigma_SSB","Mean_SSB5","Sigma_SSB5",
                                                   "Mean_C","Sigma_C","Mean_C5","Sigma_C5",
                                                   "Mean_Depl",
                                                   "Prob_Depl_75","Prob_Depl_50","Prob_Depl_40",
                                                   "Prob_Depl_30","Prob_Depl_25","Prob_Depl_20"))
  
  # Populate the table
  for(stock in 1:length(Stocks)){
    stockpath <- paste(Outpath,Stocks[stock],sep="/")
    # Loop over Fs
    for(scenario in 1:length(FFs)){
      # Get the right folder
      Fpath <- paste0(stockpath,"/F_",FFs[scenario])
      setwd(Fpath)
      # Run the summarizing fortran program
      shell(paste0(OMpath,"/Sort/Sort.exe"),intern=FALSE, wait=TRUE)
      # Get those results
      outs <- read.table("PerfInd.OUT")
      outs.vec <- as.vector(as.matrix(outs[1:7,c(2,5,3,6,4,7)]))
      # Calculate probabilities of falling below thresholds
      if(scenario==1) B0 <- outs.vec[8]
      mean.depl  <- outs.vec[8]/B0
      depl.thres <- B0*c(0.75,0.5,0.4,0.3,0.25,0.2)
      #depl.probs <- approx(x=outs[8:107,6],y=outs[8:107,1],xout=depl.thres)$y
      depl.probs <- approx(x=outs[8:106,6],y=1:99,xout=depl.thres)$y
      depl.probs <- ifelse(depl.thres<min(outs[8:107,6]),0,depl.probs)
      depl.probs <- ifelse(depl.thres>max(outs[8:107,6]),100,depl.probs)
      # Save to final tables
      plot.table <- cbind(plot.table,outs.vec)
      report.table <- cbind(report.table,c(outs.vec[c(1,2,8,9,15,16,22,23,29,30,36,37)],
                                           mean.depl,depl.probs))
      cat(Stocks[stock],"F =",FFs[scenario],"\n")
    }
  }
  
  # Add names to the final table
  final.names <- character(length(Stocks)*length(FFs))
  for(stock in 1:length(Stocks))
    final.names[((stock-1)*length(FFs)+1):(stock*length(FFs))] <- paste0(Stocks[stock],"_",FFs)
  names(plot.table) <- c("Performance_Measure",final.names)
  names(report.table) <- c("Performance_Measure",final.names)
  # Save the results
  dir.create(paste0(Outpath,"/AAA-Summary"))
  setwd(paste0(Outpath,"/AAA-Summary"))
  write.csv(plot.table,"plot.csv", row.names=F)
  write.csv(report.table,"report.csv", row.names=F)
}

##############################################################################################
# Load the final table
setwd(paste0(Outpath,"/AAA-Summary"))
plot.table <- read.csv("plot.csv")

##############################################################################################
# Summarize and plot

if(plot.how==1|plot.how==2){
  # Calculate depletion B1
  RightCol <- c(0,2,3,4,5,6)
  
  DepletionsB1 <- array(0,dim=c(length(FFs),8,6))
  DepletionsB1[,1,] <- FFs
  for(qua in 1:6){
    for(stock in 1:length(Stocks)){
      for(ff in 1:length(FFs)){
        DepletionsB1[ff,stock+1,qua] <- plot.table[8+RightCol[qua],(stock-1)*length(FFs)+ff+1]/
          plot.table[8,(stock-1)*length(FFs)+2]
      }
    }
  }
  
  # Calculate depletion SSB
  DepletionsSSB <- array(0,dim=c(length(FFs),8,6))
  DepletionsSSB[,1,] <- FFs
  for(qua in 1:6){
    for(stock in 1:length(Stocks)){
      for(ff in 1:length(FFs)){
        DepletionsSSB[ff,stock+1,qua] <- plot.table[22+RightCol[qua],(stock-1)*length(FFs)+ff+1]/plot.table[22,(stock-1)*length(FFs)+2]
      }
    }
  }
  #Catch and msy
  Totcatch <- array(0,dim=c(length(FFs),8,6))
  Totcatch[,1,] <- FFs
  Relcatch <- array(0,dim=c(length(FFs),8,6))
  Relcatch[,1,] <- FFs
  for(qua in 1:6){
    for(stock in 1:length(Stocks)){
      for(ff in 1:length(FFs)){
        Relcatch[ff,stock+1,qua] <- plot.table[36+RightCol[qua],(stock-1)*length(FFs)+ff+1]/
          max(plot.table[36,((stock-1)*length(FFs)+2):(stock*length(FFs)+1)])
        Totcatch[ff,stock+1,qua] <- plot.table[36+RightCol[qua],(stock-1)*length(FFs)+ff+1]
      }
    }
  }
  
  ##################################################################################################
  # Plot all species together
  
  if(plot.how==1){
    #####################################################################################
    # Plot everything /// Just means, no uncertainty
    if(NoUncPlot==T){
      
      #B1+
      if(B1pPlot==T){
        windows(height=Wheight, width=Wwidth)
        par(mfrow=c(2,2),mar=c(0,0,0,0),oma=c(5,4,0,0))
        
        spps <- matrix(c(1:3,3:7),ncol=4)
        spnames <- c("Redbait", "Jack Mackerel", "Blue Mackerel", "Sardine")
        
        for(spp in 1:4){
          plot(1,1, xlim=c(0,MaxF),ylim=c(0,1.1),axes=F,type='n')
          box()
          lines(DepletionsB1[,1,1],DepletionsB1[,spps[1,spp]+1,1],col="gray")
          lines(DepletionsB1[,1,1],DepletionsB1[,spps[2,spp]+1,1])
          if(spp==1|spp==3) axis(2, at=c(0,0.5,1))
          if(spp==3|spp==4) axis(1)
          text(x=-0.05,y=1.1,labels=spnames[spp],pos=4)
          if(spp==1)
            legend("topright",legend=c("East","West"),lwd=1,col=c("gray","black"),bty='n')
        }
        mtext(side=1, text="Exploitation rate", outer=T, line=2.5)
        mtext(side=2, text=expression("Depletion (B1+/B1+"[0]*")"),outer=T, line=2.5)
      }
      
      #SSB
      windows(height=Wheight, width=Wwidth)
      par(mfrow=c(2,2),mar=c(0,0,0,0),oma=c(5,4,0,0))
      
      spps <- matrix(c(1:3,3:7),ncol=4)
      spnames <- c("Redbait", "Jack Mackerel", "Blue Mackerel", "Sardine")
      
      for(spp in 1:4){
        plot(1,1, xlim=c(0,MaxF),ylim=c(0,1.1),axes=F,type='n')
        box()
        lines(DepletionsSSB[,1,1],DepletionsSSB[,spps[1,spp]+1,1],col="gray")
        lines(DepletionsSSB[,1,1],DepletionsSSB[,spps[2,spp]+1,1])
        if(spp==1|spp==3) axis(2, at=c(0,0.5,1))
        if(spp==3|spp==4) axis(1)
        text(x=-0.05,y=1.1,labels=spnames[spp],pos=4)
        if(spp==1)
          legend("topright",legend=c("East","West"),lwd=1,col=c("gray","black"),bty='n')
      }
      mtext(side=1, text="Exploitation rate", outer=T, line=2.5)
      mtext(side=2, text=expression("Depletion (SSB/SSB"[0]*")"),outer=T, line=2.5)
      
      
      #Catch
      windows(height=Wheight, width=Wwidth)
      par(mfrow=c(2,2),mar=c(0,0,0,0),oma=c(5,4,0,0))
      
      spps <- matrix(c(1:3,3:7),ncol=4)
      spnames <- c("Redbait", "Jack Mackerel", "Blue Mackerel", "Sardine")
      
      for(spp in 1:4){
        plot(1,1, xlim=c(0,MaxF),ylim=c(0,1.1),axes=F,type='n')
        box()
        lines(Relcatch[,1,1],Relcatch[,spps[1,spp]+1,1],col="gray")
        lines(Relcatch[,1,1],Relcatch[,spps[2,spp]+1,1])
        if(spp==1|spp==3) axis(2, at=c(0,0.5,1))
        if(spp==3|spp==4) axis(1)
        text(x=-0.05,y=1.1,labels=spnames[spp],pos=4)
        if(spp==1)
          legend("topright",legend=c("East","West"),lwd=1,col=c("gray","black"),bty='n')
      }
      mtext(side=1, text="Exploitation rate", outer=T, line=2.5)
      mtext(side=2, text="Catch",outer=T, line=2.5)
    }
    
    
    ##################################################################################
    # Plots with uncertainty
    
    # B1+
    if(B1pPlot==T){
      windows(height=Wheight, width=Wwidth)
      par(mfrow=c(2,2),mar=c(0,0,0,0),oma=c(5,4,0,0))
      
      spps <- matrix(c(1:3,3:7),ncol=4)
      spnames <- c("Redbait", "Jack Mackerel", "Blue Mackerel", "Sardine")
      maxY <- max(DepletionsB1[,,6])
      
      for(spp in 1:4){
        Wcol <- WColLine
        if(spp==2) Wcol <- JackCol
        plot(1,1, xlim=c(0,MaxF),ylim=c(0,maxY),axes=F,type='n')
        box()
        polygon(x=c(DepletionsB1[,1,1],rev(DepletionsB1[,1,1])),
                y=c(DepletionsB1[,spps[1,spp]+1,2],rev(DepletionsB1[,spps[1,spp]+1,6])),
                col=adjustcolor(EColPol,alpha.f=0.5),border=NA)
        polygon(x=c(DepletionsB1[,1,1],rev(DepletionsB1[,1,1])),
                y=c(DepletionsB1[,spps[2,spp]+1,2],rev(DepletionsB1[,spps[2,spp]+1,6])),
                col=adjustcolor(WColPol,alpha.f=0.5),border=NA)
        lines(DepletionsB1[,1,1],DepletionsB1[,spps[1,spp]+1,1],col=EColLine,lwd=2)
        lines(DepletionsB1[,1,1],DepletionsB1[,spps[2,spp]+1,1],col=Wcol,lwd=2)
        if(spp==1|spp==3) axis(2, at=c(0,0.5,1,1.5))
        if(spp==3|spp==4) axis(1)
        text(x=MaxF+0.05,y=1.5,labels=spnames[spp],pos=2)
        if(spp==1)
          legend("topleft",legend=c("East","West"),lwd=2,col=c(EColLine,WColLine),bty='n')
      }
      mtext(side=1, text="Exploitation rate", outer=T, line=2.5)
      mtext(side=2, text=expression("Depletion (B1+/B1+"[0]*")"),outer=T, line=2.5)    
    }
    
    # SSB
    windows(height=Wheight, width=Wwidth)
    par(mfrow=c(2,2),mar=c(0,0,0,0),oma=c(5,4,0,0))
    
    spps <- matrix(c(1:3,3:7),ncol=4)
    spnames <- c("Redbait", "Jack Mackerel", "Blue Mackerel", "Sardine")
    maxY <- max(DepletionsSSB[,,6])
    
    for(spp in 1:4){
      Wcol <- WColLine
      if(spp==2) Wcol <- JackCol
      plot(1,1, xlim=c(0,MaxF),ylim=c(0,maxY),axes=F,type='n')
      box()
      polygon(x=c(DepletionsSSB[,1,1],rev(DepletionsSSB[,1,1])),
              y=c(DepletionsSSB[,spps[1,spp]+1,2],rev(DepletionsSSB[,spps[1,spp]+1,6])),
              col=adjustcolor(EColPol,alpha.f=0.5),border=NA)
      polygon(x=c(DepletionsSSB[,1,1],rev(DepletionsSSB[,1,1])),
              y=c(DepletionsSSB[,spps[2,spp]+1,2],rev(DepletionsSSB[,spps[2,spp]+1,6])),
              col=adjustcolor(WColPol,alpha.f=0.5),border=NA)
      lines(DepletionsSSB[,1,1],DepletionsSSB[,spps[1,spp]+1,1],col=EColLine,lwd=2)
      lines(DepletionsSSB[,1,1],DepletionsSSB[,spps[2,spp]+1,1],col=Wcol,lwd=2)
      if(spp==1|spp==3) axis(2, at=c(0,0.5,1,1.5))
      if(spp==3|spp==4) axis(1)
      text(x=MaxF+0.05,y=1.5,labels=spnames[spp],pos=2)
      if(spp==1)
        legend("topleft",legend=c("East","West"),lwd=2,col=c(EColLine,WColLine),bty='n')
    }
    mtext(side=1, text="Exploitation rate", outer=T, line=2.5)
    mtext(side=2, text=expression("Depletion (SSB/SSB"[0]*")"),outer=T, line=2.5)
    
    # Catch
    windows(height=Wheight, width=Wwidth)
    par(mfrow=c(2,2),mar=c(0,0,0,0),oma=c(5,4,0,0))
    
    spps <- matrix(c(1:3,3:7),ncol=4)
    spnames <- c("Redbait", "Jack Mackerel", "Blue Mackerel", "Sardine")
    maxY <- max(Relcatch[,,6])
    
    for(spp in 1:4){
      Wcol <- WColLine
      if(spp==2) Wcol <- JackCol
      plot(1,1, xlim=c(0,MaxF),ylim=c(0,maxY),axes=F,type='n')
      box()
      polygon(x=c(Relcatch[,1,1],rev(Relcatch[,1,1])),
              y=c(Relcatch[,spps[1,spp]+1,2],rev(Relcatch[,spps[1,spp]+1,6])),
              col=adjustcolor(EColPol,alpha.f=0.5),border=NA)
      polygon(x=c(Relcatch[,1,1],rev(Relcatch[,1,1])),
              y=c(Relcatch[,spps[2,spp]+1,2],rev(Relcatch[,spps[2,spp]+1,6])),
              col=adjustcolor(WColPol,alpha.f=0.5),border=NA)
      lines(Relcatch[,1,1],Relcatch[,spps[1,spp]+1,1],col=EColLine,lwd=2)
      lines(Relcatch[,1,1],Relcatch[,spps[2,spp]+1,1],col=Wcol,lwd=2)
      if(spp==1|spp==3) axis(2, at=c(0,0.5,1,1.5,2))
      if(spp==3|spp==4) axis(1)
      text(x=-0.05,y=2,labels=spnames[spp],pos=4)
      if(spp==1)
        legend("topright",legend=c("East","West"),lwd=2,col=c(EColLine,WColLine),bty='n')
    }
    mtext(side=1, text="Exploitation rate", outer=T, line=2.5)
    mtext(side=2, text="Catch",outer=T, line=2.5)
  }
  
  ##################################################################################################
  # Plot species separately
  
  if(plot.how==2){
    
    spps <- matrix(c(1:3,3:7),ncol=4)
    spnames <- c("Redbait", "Jack Mackerel", "Blue Mackerel", "Sardine")
    
    for(spp in 1:4){
      windows(height=Wheight*0.7, width=Wwidth*1.2)
      par(mfrow=c(1,2),mar=c(0,4,0,0),oma=c(4,0,0,0.2))    
      
      # SSB
      maxY <- max(DepletionsSSB[,,6])
      Wcol <- WColLine
      if(spp==2) Wcol <- JackCol
      plot(1,1, xlim=c(0,MaxF),ylim=c(0,maxY),axes=F,type='n', ylab=NA)
      box()
      polygon(x=c(DepletionsSSB[,1,1],rev(DepletionsSSB[,1,1])),
              y=c(DepletionsSSB[,spps[1,spp]+1,2],rev(DepletionsSSB[,spps[1,spp]+1,6])),
              col=adjustcolor(EColPol,alpha.f=0.5),border=NA)
      polygon(x=c(DepletionsSSB[,1,1],rev(DepletionsSSB[,1,1])),
              y=c(DepletionsSSB[,spps[2,spp]+1,2],rev(DepletionsSSB[,spps[2,spp]+1,6])),
              col=adjustcolor(WColPol,alpha.f=0.5),border=NA)
      lines(DepletionsSSB[,1,1],DepletionsSSB[,spps[1,spp]+1,1],col=EColLine,lwd=2)
      lines(DepletionsSSB[,1,1],DepletionsSSB[,spps[2,spp]+1,1],col=Wcol,lwd=2)
      axis(2, at=c(0,0.5,1,1.5))
      axis(1)
      if(spp!=2)
        legend("topright",legend=c("East","West"),lwd=2,col=c(EColLine,WColLine),bty='n')
      mtext(side=2, text=expression("Depletion (SSB/SSB"[0]*")"), line=2.0)
      
      # Catch
      maxY <- max(Relcatch[,,6])
      Wcol <- WColLine
      if(spp==2) Wcol <- JackCol
      plot(1,1, xlim=c(0,MaxF),ylim=c(0,maxY),axes=F,type='n', ylab=NA)
      box()
      polygon(x=c(Relcatch[,1,1],rev(Relcatch[,1,1])),
              y=c(Relcatch[,spps[1,spp]+1,2],rev(Relcatch[,spps[1,spp]+1,6])),
              col=adjustcolor(EColPol,alpha.f=0.5),border=NA)
      polygon(x=c(Relcatch[,1,1],rev(Relcatch[,1,1])),
              y=c(Relcatch[,spps[2,spp]+1,2],rev(Relcatch[,spps[2,spp]+1,6])),
              col=adjustcolor(WColPol,alpha.f=0.5),border=NA)
      lines(Relcatch[,1,1],Relcatch[,spps[1,spp]+1,1],col=EColLine,lwd=2)
      lines(Relcatch[,1,1],Relcatch[,spps[2,spp]+1,1],col=Wcol,lwd=2)
      axis(2, at=c(0,0.5,1,1.5,2))
      axis(1)
      mtext(side=1, text="Exploitation rate", outer=T, line=2.5)
      mtext(side=2, text="Catch", line=2.0)
    }
  }
  
  ##########################################################
  #Fmsy
  Emsy <- data.frame(Stock=Stocks, 
                     Emsy=Relcatch[,1,1][which(Relcatch[,,1]==1,arr.ind=T)[2:8,1]],
                     Emsy.catch=Totcatch[,,1][which(Relcatch[,,1]==1,arr.ind=T)[2:8,]],
                     Emsy.depl=round(DepletionsSSB[,,1][which(Relcatch[,,1]==1,arr.ind=T)[2:8,]],digits=2)
  )
  print("Emsy values for all stocks")
  print(Emsy)
  
  ############################################################
  # Depletions
  
  DepObj <- c(0.9,0.8,0.75,seq(0.7,0.4,by=-0.1))
  EspDepls <- data.frame(Targ_Depl = DepObj)
  for(spp in 1:7){
    EspDepls <- cbind(EspDepls,
                      round(approx(x=DepletionsSSB[,spp+1,1],y=DepletionsSSB[,1,1],xout=DepObj)$y,3))
  }
  names(EspDepls) <- c("Targ_Depl",Stocks)
  print("Harvest rates to achieve depletion objectives")
  print(EspDepls)
}
