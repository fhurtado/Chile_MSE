      INTEGER MaxYr,MaxSim,MaxDim,MaxVec
      PARAMETER(MaxYr=1000,MaxSim=1000,MaxDim=10000,MaxVec=60000)
	  
	  
      COMMON / Results / Resmat,Biomass,SSB,Catches,Bio5,SSB5,Catch5,
     +                   Tres,ADim,Nsim,Nyear
      REAL*8 Resmat(0:MaxYr,0:MaxDim),Biomass(0:MaxVec),SSB(0:MaxVec)
      REAL*8 Catches(0:MaxVec),Bio5(0:MaxVec),SSB5(0:MaxVec)
      REAL*8 Catch5(0:MaxVec),Tres(0:MaxDim)
      INTEGER ADim,Nsim,Nyear
	  