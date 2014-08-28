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
read.again.for <- F # Read the output using fortran?
plot.how       <- 2 # Plot? 0 = No, 1 = Together, 2 = Sp by Sp

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
Outpath <- "C:/Users/Felipe/Desktop/Chile_MSE/FvD_B1_Surv"
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
  RunOM <- function(FFs,Outpath,OMpath,DATfile){
    # Create the folders needed
    Fpath <- paste0(Outpath,"/F_",FFs)
    dir.create(Fpath)
    setwd(Fpath)
    # Copy the appropriate files
    file.copy(from=paste(OMpath,DATfile,sep="/"),to=paste(Fpath,DATfile,sep="/"),overwrite=T)
    file.copy(from=paste(OMpath,"Chile_OM.SPEC",sep="/"),to=paste(Fpath,"Chile_OM.SPEC",sep="/"),overwrite=T)
    file.copy(from=paste(OMpath,"SPFspp.SPEC",sep="/"),to=paste(Fpath,"SPFspp.SPEC",sep="/"),overwrite=T)
    # Modify the DAT file
    t.datfile <- readLines(DATfile)
    # Modify the SPEC file
    t.specfile <- readLines("Chile_OM.SPEC")
    t.specfile[6] <- paste(FFs,0.7,sep="\t")
    writeLines(t.specfile,"Chile_OM.SPEC")
    # Modify the spp file
#     t.sppfile <- readLines("SPFspp.SPEC")
#     t.sppfile[2] <- SppN[stock]
#     t.sppfile[4] <- StockN[stock]
#     writeLines(t.sppfile,"SPFspp.SPEC")
    # Run the model 
    shell(paste0(OMpath,"/Chile_OM.exe"),intern=FALSE, wait=TRUE)
  }
  
  
  # Loop over stocks and run the function above
#  for(stock in 1:length(Stocks)){
    #stockpath <- paste(Outpath,Stocks[stock],sep="/")
    #dir.create(stockpath)
    
    #Parallel part
    sfInit(parallel=TRUE, cpus=Ncpus, type="SOCK")
    sfSapply(FFs,RunOM,Outpath,OMpath,DATfile)
    sfStop()   
    
#    cat("Finished ",Stocks[stock],"\n")
#  } #Stock loop close
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
#  for(stock in 1:length(Stocks)){
#    stockpath <- paste(Outpath,Stocks[stock],sep="/")
    # Loop over Fs
    for(scenario in 1:length(FFs)){
      # Get the right folder
      Fpath <- paste0(Outpath,"/F_",FFs[scenario])
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
      cat("F =",FFs[scenario],"\n")
    }
#  }
  
  # Add names to the final table
  final.names <- paste0("F_",FFs)
  names(plot.table) <- c("Performance_Measure",final.names)
  names(report.table) <- c("Performance_Measure",final.names)
  # Save the results
  dir.create(paste0(Outpath,"/AAA-Summary"), showWarnings = FALSE)
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
  
  DepletionsB1 <- array(0,dim=c(7,length(FFs)))
  DepletionsB1[1,] <- FFs
  for(qua in 1:6){
    DepletionsB1[qua+1,] <- as.numeric(plot.table[8+RightCol[qua],2:(length(FFs)+1)]/
      plot.table[8,2])
  }
  
  # Calculate depletion SSB
  DepletionsSSB <- array(0,dim=c(7,length(FFs)))
  DepletionsSSB[1,] <- FFs
  for(qua in 1:6){
    for(ff in 1:length(FFs)){
      DepletionsSSB[qua+1,] <- as.numeric(plot.table[22+RightCol[qua],2:(length(FFs)+1)]/
        plot.table[22,2])
    }
  }

  #Catch and msy
  Totcatch <- array(0,dim=c(7,length(FFs)))
  Totcatch[1,] <- FFs
  Relcatch <- array(0,dim=c(7,length(FFs)))
  Relcatch[1,] <- FFs
  for(qua in 1:6){
    Relcatch[qua+1,] <- as.numeric(plot.table[36+RightCol[qua],2:(length(FFs)+1)]/
      max(plot.table[36,2:(length(FFs)+1)]))
    Totcatch[qua+1,] <- as.numeric(plot.table[36+RightCol[qua],2:(length(FFs)+1)])
  }
  
  ##################################################################################################
  # Plot species separately
  
  if(plot.how==2){
    
#     spps <- matrix(c(1:3,3:7),ncol=4)
#     spnames <- c("Redbait", "Jack Mackerel", "Blue Mackerel", "Sardine")
#     
#     for(spp in 1:4){
      windows(height=Wheight*0.7, width=Wwidth*1.2)
      par(mfrow=c(1,2),mar=c(0,4,0,0),oma=c(4,0,0,0.2))    
      
      # B1+
      maxY <- max(DepletionsB1[7,])
      Wcol <- WColLine
      plot(1,1, xlim=c(0,MaxF),ylim=c(0,maxY),axes=F,type='n', ylab=NA, yaxs='i')
      box()
      polygon(x=c(DepletionsB1[1,],rev(DepletionsB1[1,])),
              y=c(DepletionsB1[3,],rev(DepletionsB1[7,])),
              col=adjustcolor(WColPol,alpha.f=0.5),border=NA)
      polygon(x=c(DepletionsB1[1,],rev(DepletionsB1[1,])),
              y=c(DepletionsB1[4,],rev(DepletionsB1[6,])),
              col=adjustcolor(WColPol,alpha.f=0.8),border=NA)
      lines(DepletionsB1[1,],DepletionsB1[2,],col=Wcol,lwd=2)
      axis(2, at=c(0,0.5,1,1.5))
      axis(1)
      mtext(side=2, text=expression("Depletion (B1+/B1+"[0]*")"), line=2.0)
      
      # Catch
      maxY <- max(Relcatch[7,])
      Wcol <- WColLine
      plot(1,1, xlim=c(0,MaxF),ylim=c(0,maxY),axes=F,type='n', ylab=NA, yaxs='i')
      box()
      polygon(x=c(Relcatch[1,],rev(Relcatch[1,])),
              y=c(Relcatch[3,],rev(Relcatch[7,])),
              col=adjustcolor(WColPol,alpha.f=0.5),border=NA)
      polygon(x=c(Relcatch[1,],rev(Relcatch[1,])),
              y=c(Relcatch[4,],rev(Relcatch[6,])),
              col=adjustcolor(WColPol,alpha.f=0.8),border=NA)
      lines(Relcatch[1,],Relcatch[2,],col=Wcol,lwd=2)
      axis(2, at=c(0,0.5,1,1.5,2))
      axis(1)
      mtext(side=1, text="Exploitation rate", outer=T, line=2.5)
      mtext(side=2, text="Catch", line=2.0)
#    }
  }
  
  ##########################################################
  #Fmsy
  Emsy <- data.frame(Emsy=Relcatch[1,][which(Relcatch[2,]==1,arr.ind=T)],
                     Emsy.catch=Totcatch[2,][which(Relcatch[2,]==1,arr.ind=T)],
                     Emsy.depl=round(DepletionsB1[2,][which(Relcatch[2,]==1,arr.ind=T)],digits=2)
  )
  print("Emsy values for all stocks")
  print(Emsy)
  
  ############################################################
  # Depletions
  
  DepObj <- c(seq(0.9,0.1,by=-0.1))
  EspDepls <- data.frame(Targ_Depl = DepObj)
  EspDepls <- cbind(EspDepls,
                      round(approx(x=DepletionsSSB[2,],y=DepletionsSSB[1,],xout=DepObj)$y,3))

  names(EspDepls) <- c("Targ_Depl","Harvest_Rate")
  print("Harvest rates to achieve depletion objectives")
  print(EspDepls)
}
