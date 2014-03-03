      INTEGER MaxYr,MaxAge,MaxSim
      PARAMETER(MaxYr=1000,MaxAge=25,MaxSim=1000)
C
      COMMON / General / Nsim,Spp,EorW
      INTEGER Nsim,Spp,EorW
C      Nsim                  - Number of simulations  
C      Spp                   - Which species parameters file is used 
C 	   EorW					 - Use East or West stock parameters?
      
      COMMON / BiolPar / Nage,Nyear,Weight,Fecu,Matu,Selex,Niprop,NatM,
     +     Rstar,R0,Ninit,DeplI,SSB0,Steepness,Alpha,Beta1,Beta2,Rerr,
     +     Rtype,AlogRPS,BlogRPS,PhilogRPS 
      INTEGER Nyear,Nage
      REAL*8 Weight(0:MaxAge),Fecu(0:MaxAge),Matu(0:MaxAge)
      REAL*8 Selex(0:MaxAge)
      REAL*8 SelxPNW(0:MaxAge),Niprop(0:MaxAge)
      REAL*8 NatM,Rstar,R0,Ninit,DeplI,SSB0,Steepness,Alpha,Beta1,Beta2
      REAL*8 Rerr,Rtype,AlogRPS,BlogRPS,PhilogRPS	  
C      Nyear                 - Number of projection years
C      Nage                  - Number of age-classes year      
C      Weight(a)             - Weight-at-age      
C      Matu(a)               - Maturity
C      Fecu(a)               - Fecundity
C      Selex(a)              - Selectivity
C      Niprop(a)             - Initial numbers-at-age (proportions)
C      NatM                  - Natural mortality
C      R0                    - Virgin recruitment
C      Rstar                 - Maximum recruitment
C  	   Ninit				 - Initial numbers-at-age (total)
C      DepI 				 - Initial depletion
C	   SSB0					 - Virgin SSB
C      Steepness             - Steepness
C      Alpha,Beta1,Beta2     - SR parameters
C      Rerr                  - Use Recruitment uncertainty?
C      Rtype                 - What rec model to use (0=dome shaped; 1=flat top)

      COMMON / Dynamics / N,B1P,SSB,F,ObsF,Catch,RecRep,AgeN,AgeC
      REAL*8 N(0:MaxYr,0:MaxAge),B1p(MaxSim,0:MaxYr),SSB(MaxSim,0:MaxYr)
      REAL*8 F(MaxSim,0:MaxYr),ObsF(MaxSim,0:MaxYr)
      REAL*8 Catch(MaxSim,0:MaxYr),RecRep(MaxSim,0:MaxYr)
      REAL*8 AgeN(MaxSim,0:MaxYr),AgeC(MaxSim,0:MaxYr)
C      N(y,a)                - Numbers-at-age      
C      B1p(s,y)              - Biomass of 1+ animals for year y and simulation s
C      SSB(s,y)              - SSB for year y and simulation s
C      F(s,y)                - True fishing mortality
C      ObsF(s,y)             - Fishing mortality from RBC
C      Catch(s,y)            - Catch by year
C      RecRep(s,y)           - Recruitment by year
C      AgeN(s,y)             - Mean Age of the population
C      AgeC(s,y)             - Mean age of the catch
C      
	  
      COMMON / Assessment / Bhat,Emsymin,Emsymax,Cutoff,FracMin,FracMax,
     +       Maxcatch,OFL,Pstar,OFLmax,SigmaB,Berr,CutOpt,EpsB,T1psi,
     +       T1gamma,T1Max,T1Min,T2Max,T3Max,Tier,OFLopt,HCRopt,ISEED4
      REAL*8 Bhat,Emsymin,Emsymax,Cutoff,FracMin,FracMax,Maxcatch,OFL
      REAL*8 OFLmax,Pstar,SigmaB,Berr,CutOpt,T1psi,T1gamma,T1Max,T1Min
      REAL*8 T2Max,T3Max,EpsB(0:MaxYr)
      INTEGER Tier,HCRopt,OFLopt,ISEED4
C      Bhat                  - Estimated biomass
C      Emsymin               - Min value for Emsy
C      Emsymax               - Max value for Emsy
C      Cuttoff               - Cutoff parameter of the RBC
C      FracMin               - Min fraction parameter of the RBC
C      FracMax               - Max fraction parameter of the RBC
C      Maxcatch              - Maxcatch parameter of the RBC
C      OFL                   - OFL (overfishing limit)
C      Pstar                 - P star value for ABC calculation
C      SigmaB                - Variance for biomass estimate
C      Berr                  - Use recruitment error?
C      T1psi                 - Psi parameter for Tier 1
C      T1gamma               - Gamma parameter for Tier 1
C      T1Max                 - Max harvest rate for Tier 1
C      T2Max                 - Max catch for Tier 2
C      T3Max                 - Max catch for Tier 3
C      Tier                  - What tier to use?
C      HCRopt                - Which HCR should be used? 
C      OFLopt                - Which OFL should be used? 

      COMMON / EnvVar / Gtype,ISEED1,ISEED2,ISEED3,G,V,Eps,I,Vdev,
     +       EpsInit,ProwR,SigmaR,Vinit,Ybar,Period,Chi,Power,ProwV,
     +       SigmaV,Phi,MeanSST,SigmaI,Ierr,Verr,UseG,Amplitude,Vrep   
      INTEGER Gtype
      INTEGER ISEED1,ISEED2,ISEED3
      REAL*8 G(0:MaxYr),V(0:MaxYr),Eps(0:MaxYr),I(0:MaxYr),Vdev(0:MaxYr)
      REAL*8 EpsInit,ProwR,SigmaR
      REAL*8 Vinit,Ybar,Period,Chi,Power,ProwV,SigmaV,Phi,MeanSST
      REAL*8 SigmaI,Ierr,Verr,UseG,Amplitude
      REAL*8 Vrep(MaxSim,0:MaxYr)
C      Gtype                 - Type of variables      
C      G(y)                  - Underlying trends in the Environment variables
C      V(y)                  - Environmental variable
C      I(y)                  - Observed environmental variable
C      Eps(y)                - Rec_dev
C      EpsInit               - last Eps
C      ProwR                 - Prow for rec_dev
C      SigmaR                - Sigma for rec_dev
C      Vinit                 - last V
C      Ybar                  - Centering parameter
C      Period                - Period of the signal
C      Chi                   - Parameter for Type 1 relationship
C      Power                 - Parameter for Type 2 relationship
C      SigmaI                - Sigma for measurement error on the environmental index
C      Verr                  - Use environment uncertainty?
C
      COMMON / Sensitivities /  TVselx,TVweight,MultFleet,HyperSt,TVM,
     +       HPSa,Fmex,Fcan,WeightOut,SelxOut,UseAnn,DiffSelx,WhatAssmt,
     +       VarSurvFreq,DEPMfrq,SurvBias,SelxSum,WeightSum,SimpleSelex,
     +       SelxMC1,SelxMC2,SelxMCy,SelxPNW,Weight80s,Weight90s,Myear
      REAL*8 TVselx,TVweight,MultFleet,HyperSt,TVM,HPSa,Fmex,Fcan
      REAL*8 WeightOut,SelxOut,UseAnn,DiffSelx,WhatAssmt,VarSurvFreq
      REAL*8 DEPMfrq,SurvBias
      REAL*8 SelxSum(0:MaxYr,0:MaxAge),WeightSum(0:MaxYr,0:MaxAge)
      REAL*8 SimpleSelex(0:MaxAge),SelxMC1(0:MaxAge),SelxMC2(0:MaxAge)
      REAL*8 SelxMCy(0:MaxAge),Weight80s(0:MaxAge),Weight90s(0:MaxAge)
      REAL*8 Myear(0:MaxYr)
	  