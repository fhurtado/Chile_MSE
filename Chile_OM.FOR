C Codigo para simulaciones del EEM Chileno. Escrito por F. Hurtado-Ferro, A.E. Punt y R. Licandeo
C Por favor referencie a los autores si utiliza este codigo.
C
      IMPLICIT NONE
      INCLUDE "Chile_OM.INC"
C
C     Local variables 
      INTEGER Year,Isim
C      
C     Read in the specifications
      CALL ReadInBiol()
C
C     Read in the parameters of the control rule
C
      OPEN(UNIT=14,FILE="SUMMARY.OUT")
C	
C     Generate the underlying environmental signal (common for all simulations)  
C     NO NEED
      CALL CalcG()
C
      DO 10000 Isim = 1,Nsim
C
C      Initialize
       CALL Initialize(Isim)
C
       DO 11000 Year = 1,Nyear
        CALL Proj1Yr(Isim,Year)
11000  CONTINUE        
C
10000 CONTINUE   
      WRITE(*,*) "Simulation finished running"
      WRITE(*,*) "Writing summary output"  
      CALL SUMMARY()
      CLOSE(14)
C
      STOP
      END
C
C ===============================================================================================
C
      SUBROUTINE ReadInBiol()
C
      IMPLICIT NONE
      INCLUDE "Chile_OM.INC"
C
C     Local variables
      INTEGER Age
C
C     Which species is going to be used
      OPEN(UNIT=13,FILE="SPFspp.spec")
      READ(13,*)
      READ(13,*) Spp
      READ(13,*)
      READ(13,*) EorW
      CLOSE(13)
C
C     Open the right species file
      IF(Spp.EQ.1) THEN 
       WRITE(*,*) "Simulation for: REDBAIT"
       OPEN(UNIT=13,FILE="Chile_OM.Dat")
       IF(EorW.EQ.1) THEN 
        CLOSE(13)
        OPEN(UNIT=13,FILE="W.Redbait.Dat")
       ENDIF
      ENDIF
      IF(Spp.EQ.2) THEN
       WRITE(*,*) "Simulation for: JACK MACKEREL"
       OPEN(UNIT=13,FILE="JackMac.Dat")
      ENDIF
      IF(Spp.EQ.3) THEN 
      WRITE(*,*) "Simulation for: BLUE MACKEREL"
       OPEN(UNIT=13,FILE="E.BlueMac.Dat")
       IF(EorW.EQ.1) THEN 
        CLOSE(13)
        OPEN(UNIT=13,FILE="W.BlueMac.Dat")
       ENDIF
      ENDIF
      IF(Spp.EQ.4) THEN 
       WRITE(*,*) "Simulation for: SARDINE"
       OPEN(UNIT=13,FILE="E.Sardine.Dat")
       IF(EorW.EQ.1) THEN 
        CLOSE(13)
        OPEN(UNIT=13,FILE="W.Sardine.Dat")
       ENDIF
      ENDIF
C      OPEN(UNIT=13,FILE="Sardine.Dat")	  
C
C     And read from that file
      READ(13,*)
      READ(13,*) Nage
      WRITE(*,'(A5,1x,I3)') "Nage",Nage
      READ(13,*)
      READ(13,*) (Weight90s(Age),Age=0,Nage)
C      READ(13,*)
C      READ(13,*) (Weight80s(Age),Age=0,Nage)
      READ(13,*)
      READ(13,*) (Matu(Age),Age=0,Nage)
      READ(13,*)
      READ(13,*) (Fecu(Age),Age=0,Nage)
      READ(13,*)
      READ(13,*) (Niprop(Age),Age=0,Nage)
      READ(13,*)
      READ(13,*) Ninit, DeplI
      READ(13,*)
      READ(13,*) (Selex(Age),Age=0,Nage)
      READ(13,*)
      READ(13,*) (SimpleSelex(Age),Age=0,Nage)
C      READ(13,*)
C      READ(13,*) (SelxMC1(Age),Age=0,Nage)
C      READ(13,*)
C      READ(13,*) (SelxMC2(Age),Age=0,Nage)
C      READ(13,*)
C      READ(13,*) (SelxPNW(Age),Age=0,Nage)
      READ(13,*)
      READ(13,*) NatM
      READ(13,*)
      READ(13,*) Rstar, Steepness
C      READ(13,*)
C      READ(13,*) Alpha, Beta1, Beta2
      READ(13,*)
      READ(13,*) EpsInit,ProwR,SigmaR,Rerr,Rtype
      READ(13,*)
      READ(13,*) Vinit,Gtype,UseG,Verr
      READ(13,*)
      READ(13,*) Ybar,Period,ProwV,SigmaV,Phi,MeanSST,Amplitude
      READ(13,*)
      READ(13,*) Chi
      READ(13,*)
      READ(13,*) Power
      READ(13,*)
      READ(13,*) SigmaI,Ierr
      READ(13,*)
      READ(13,*) SigmaB, Berr
      READ(13,*)
      READ(13,*) T2Max
      CLOSE(13)
C
C     Read management specifications
      OPEN(UNIT=13,FILE="Chile_OM.Spec")
      READ(13,*)
      READ(13,*) Nyear
      READ(13,*)
      READ(13,*) Nsim
      READ(13,*) 
      READ(13,*) Emsymin, Emsymax	  
      READ(13,*) 
      READ(13,*) FracMin, FracMax
      READ(13,*) 
      READ(13,*) Cutoff
      READ(13,*) 
      READ(13,*) Maxcatch
      READ(13,*) 
      READ(13,*) T1psi, T1gamma, T1Max
      READ(13,*) 
      READ(13,*) T3Max
      READ(13,*) 
      READ(13,*) OFLopt
      READ(13,*) 
      READ(13,*) Tier
      READ(13,*) 
      READ(13,*) CutOpt
      READ(13,*)
      READ(13,*) Pstar,OFLmax
      READ(13,*)
      READ(13,*) ISEED1, ISEED2, ISEED3, ISEED4
      READ(13,*) 
      READ(13,*) WhatAssmt
      READ(13,*) 
      READ(13,*) DiffSelx
      READ(13,*) 
      READ(13,*) VarSurvFreq, DEPMfrq, SurvBias
      READ(13,*) 
      READ(13,*) TVselx
      READ(13,*) 
      READ(13,*) TVweight
      READ(13,*) 
      READ(13,*) MultFleet
      READ(13,*) 
      READ(13,*) Fmex,Fcan
      READ(13,*) 
      READ(13,*) HyperSt, HPSa
      READ(13,*) 
      READ(13,*) TVM
      READ(13,*) 
      READ(13,*) SelxOut, WeightOut
      READ(13,*) 
      READ(13,*) UseAnn	  
      CLOSE(13)
C      
      RETURN
      END
C
C ===============================================================================================
C
      SUBROUTINE Initialize(Isim)
C
C This subroutine sets up variables which are common across a run
C
      IMPLICIT NONE
      INCLUDE "Chile_OM.INC"
C
C     Global variables
      INTEGER Isim
C
C     Local variables
      INTEGER Iage,Year
      REAL*8 Neq(0:MaxAge),S0,EpsV,EpsR
      REAL*8 SSBV,Bio1P,Idev
      REAL*8 XNORM
      REAL*8 ANcount1, ANcount2
      REAL*8 Lselx,Jselx,Qweight,MaxSelx
      REAL*8 Gslope
      EXTERNAL XNORM
C
C     Initial N-vector
      ANcount1 = 0      
      ANcount2 = 0
      AgeN(Isim,0) = 0
      DO 10000 Iage = 0,Nage
       N(0,Iage) = Niprop(Iage)*Ninit*DeplI
       ANcount1 = ANcount1 + N(0,Iage)*Iage
       ANcount2 = ANcount2 + N(0,Iage)
10000 CONTINUE
      AgeN(Isim,0) = ANcount1/ANcount2
      RecRep(Isim,0) = N(0,0)       
C
C     By default, set the weight vector equal to that of 1991-2000
      DO 11000 Iage = 0,Nage
       Weight(Iage) = Weight90s(Iage)
11000 CONTINUE
C
C     Compute the SSB and 1+ biomass
      SSBV = 0
      Bio1P = 0
      DO 12000 Iage = 1,Nage
       Bio1P = Bio1P + Weight(Iage)*N(0,Iage)
12000 CONTINUE
      DO 13000 Iage = 0,Nage
       SSBV = SSBV + Weight(Iage)*N(0,Iage)*Matu(Iage)
13000 CONTINUE
      SSB(Isim,0) = SSBV
      B1p(Isim,0) = Bio1P
C
C     Compute Numbers at equilibrium and SSB0
      R0 = Rstar*(5*Steepness-1)/(4*Steepness)
      Neq(0) = R0
C      WRITE(*,'(A4,1x,F5.0)') " R0:",R0
      DO 14000 Iage = 1,Nage-1
       Neq(Iage) = Neq(Iage-1)*EXP(-NatM)
14000 CONTINUE
      Neq(Nage) = Neq(Nage-1)*EXP(-NatM)/(1-EXP(-NatM))
      SSB0 = 0 
      DO 15000 Iage = 0,Nage
       SSB0 = SSB0 + Weight(Iage)*Neq(Iage)*Matu(Iage)
15000 CONTINUE
C	  
C     Generate the environmental variables
C     NO NEED
C      Vdev(0) = Vinit
C      V(0) = Vdev(0) + MeanSST
C      Vrep(Isim,0) = V(0)
CC	  
CC     Update
C      DO 31000 Year = 1,Nyear
C       EpsV = XNORM(1,0.0d0,1.0d0,ISEED1)*SigmaV
C       EpsV = EpsV * Verr
C       Vdev(Year) = ProwV*Vdev(Year-1)+(1-ProwV)*G(Year)*Amplitude
C       Vdev(Year) = Vdev(Year) + SQRT(1-ProwV**2)*EpsV
C       V(Year) = Vdev(Year)*UseG + MeanSST
C       Vrep(Isim,Year) = V(Year)
C31000 CONTINUE
C
C     Generate the rec_dev
      Eps(0) = EpsInit
      DO 32000 Year = 1,Nyear
       EpsR = XNORM(2,0.0d0,1.0d0,ISEED2)*SigmaR
       EpsR = EpsR * Rerr
       Eps(Year) = prowR*Eps(Year-1)+SQRT(1-ProwR**2.0)*EpsR
32000 CONTINUE
C      
C     Add measurement error to the environmental variable
C     NO NEED
C      DO 33000 Year = 0,Nyear
C       Idev = XNORM(1,0.0d0,1.0d0,ISEED3)*SigmaI*Ierr
C       I(Year) = V(Year)+Idev
C33000 CONTINUE
C 
C     Add measurement error to the 1+ biomass (i.e. an "assessment")
      EpsB(0) = 1
      DO 34000 Year = 1,NYear   
       EpsB(Year) = (XNORM(1,0.0d0,1.0d0,ISEED4)+SurvBias)*SigmaB- 
     +                  (SigmaB**2)/2.0
       EpsB(Year) = EXP(EpsB(Year)*Berr)
34000 CONTINUE	   
C     
C     Set the default for weight and selex values
      DO 40000 Year = 0,Nyear
       DO 41000 Iage = 0,Nage
        SelxSum(Year,Iage) = Selex(Iage)
        WeightSum(Year,Iage) = Weight(Iage)
41000  CONTINUE
40000 CONTINUE
C
C     SENSITIVITIES
C     1. Time varying weight at age
      IF (TVweight.EQ.1) THEN
       DO 50000 Year=0,Nyear
C       Compute Q
        Qweight = -33.245934 + 2.115014*V(Year)
        Qweight = MIN(1.0,Qweight)
        Qweight = MAX(0.0,Qweight)
C       Compute the weight-at-age vector
        DO 51000 Iage = 0,Nage
         Weight(Iage) = Qweight*Weight80s(Iage)+
     +                  (1-Qweight)*Weight90s(Iage)
51000   CONTINUE
C       And save it
        DO 52000 Iage = 0,Nage
         WeightSum(Year,Iage) = Weight(Iage)
52000   CONTINUE
C
C       For debugging
C        WRITE(*,'(I3,1x,F5.3,1x,20(F5.3,1x))') Year,Qweight,
C     +                                         (Weight(Iage),Iage=0,Nage)
C        WRITE(*,'(I3,1x,F6.3,1x,F5.3)') Year,V(Year),Qweight
C        WRITE(*,'(I3,1x,F4.2)') Year, Myear(Year)
50000  CONTINUE
C
C      Compute the SSB and 1+ biomass
       SSBV = 0
       Bio1P = 0
       DO 53000 Iage = 1,Nage
        Bio1P = Bio1P + WeightSum(0,Iage)*N(0,Iage)
53000  CONTINUE
       DO 54000 Iage = 2,Nage
        SSBV = SSBV + WeightSum(0,Iage)*N(0,Iage)
54000  CONTINUE
       SSB(Isim,0) = SSBV
       B1p(Isim,0) = Bio1P
C
      ENDIF
C
C     SENSITIVITIES
C     2. Time varying selectivities
      IF (TVselx.EQ.1) THEN
       DO 60000 Year=0,Nyear
C       Compute L
        Lselx = -31.407868 + 1.998082*V(Year)
        Lselx = MIN(1.0,Lselx)
        Lselx = MAX(0.0,Lselx)
C       Compute J
        Jselx = -23.790353 + 1.518223*V(Year)
        Jselx = MIN(1.0,Jselx)
        Jselx = MAX(0.0,Jselx)
C       Compute the MexCal selectivity for the year
        DO 61000 Iage = 0,Nage
         SelxMCy(Iage) = Lselx*SelxMC1(Iage) + 
     +                  (1-Lselx)*SelxMC2(Iage)
61000   CONTINUE
C       Normalize (set the max selx to 1)
        MaxSelx = MAXVAL(SelxMCy)
        DO 62000 Iage = 0,Nage
         SelxMCy(Iage) = SelxMCy(Iage)/MaxSelx
62000   CONTINUE
C       Now, compute the total selectivity
        DO 63000 Iage = 0,Nage
         Selex(Iage) = Jselx*SelxMCy(Iage) +
     +                (1-Jselx)*SelxPNW(Iage)
63000   CONTINUE
C       And normalize again
        MaxSelx = MAXVAL(Selex)
        DO 64000 Iage = 0,Nage
         Selex(Iage) = Selex(Iage)/MaxSelx
64000   CONTINUE
C
C       Save the selex values
        DO 65000 Iage = 0,Nage
         SelxSum(Year,Iage) = Selex(Iage)
65000   CONTINUE
60000  CONTINUE
      ENDIF
C
C     SENSITIVITIES
C     4. Time varying natural mortality
      IF(TVM.EQ.1) THEN
       DO 70000 Year=0,Nyear
        Gslope = G(Year)-G(Year-5)
        IF(Year.LT.5) Gslope = 0.0
        Myear(Year) = 0.4
        IF(Gslope.LT.0.0) Myear(Year)=0.8
70000  CONTINUE
      ENDIF
C
C     SENSITIVITIES
C     5. Different selectivities without time varying 
      IF(DiffSelx.EQ.1) THEN
       DO 40001 Year = 0,Nyear
        DO 41001 Iage = 0,Nage
         SelxSum(Year,Iage) = SimpleSelex(Iage)
41001   CONTINUE
40001  CONTINUE
      ENDIF
      IF(DiffSelx.EQ.2) THEN
       DO 40002 Year = 0,Nyear
        DO 41002 Iage = 0,Nage
         SelxSum(Year,Iage) = Matu(Iage)
41002   CONTINUE
40002  CONTINUE
      ENDIF
C      IF(DiffSelx.EQ.1) THEN
C       DO 40001 Year = 0,Nyear
C        DO 41001 Iage = 0,Nage
C         SelxSum(Year,Iage) = SelxMC1(Iage)
C41001   CONTINUE
C40001  CONTINUE
C      ENDIF
C      IF(DiffSelx.EQ.2) THEN
C       DO 40002 Year = 0,Nyear
C        DO 41002 Iage = 0,Nage
C         SelxSum(Year,Iage) = SelxPNW(Iage)
C41002   CONTINUE
C40002  CONTINUE
C      ENDIF
C
C     For debugging only
C      DO 50000 Year = 1,Nyear
C      WRITE(*,600) Year,G(Year),V(Year),I(Year),Eps(Year)
C50000 CONTINUE       
C
      RETURN
600   FORMAT(1x,I4,1x,A1,1x,7(F7.3,1x))      
      END
C
C ===============================================================================================
C
      SUBROUTINE Proj1Yr(Isim,Iyr)
C
C This subroutine updates the basic dynamics
C
      IMPLICIT NONE
      INCLUDE "Chile_OM.INC"
C
C     Global variables
      INTEGER Iyr,Isim
      REAL*8 RBC
      REAL*8 RecBioCatch	  
      REAL*8 Eff, ObsEff
      REAL*8 Yrcount
C
C     Local variables
      INTEGER Age,II
      REAL*8 FF,Fmin,Fmax,Cnum,Cpred
      REAL*8 SSBV,Bio1P,RecR,SSBflat,ExplB
      REAL*8 ANcount1, ANcount2, Ncheck
      REAL*8 ACcount1, ACcount2
      REAL*8 qind
      REAL*8 Z(0:MaxAge)
      REAL*8 RickerA, RickerB, TBhat
      RickerA = -4.598
      RickerB = -3.457e-07
      IF(Rtype.EQ.4) THEN
       RickerA = -3.6978
       RickerB = -1.4680e-06	
      ENDIF
C
C     Do the "Assessment"
      TBhat = B1p(Isim,Iyr-1)*EpsB(Iyr-1)
C     SENSITIVITIES
C     Use SSB for the assessment (i.e. DEPM)
      IF (WhatAssmt.EQ.1.0) TBhat = SSB(Isim,Iyr-1)*EpsB(Iyr-1)
C     Use exploitable biomass - Not an assessment, just for explorations
      IF (WhatAssmt.EQ.2.0) THEN
        ExplB = 0
        DO 10100 Age = 0,Nage
          ExplB = ExplB + N(Iyr-1,Age)*SelxSum(Iyr-1,Age)*
     +            WeightSum(Iyr-1,Age)
10100   CONTINUE
        TBhat = ExplB
      ENDIF
	  
C     SENSITIVITIES
C     Variable frequency for DEPM surveys
      Yrcount = REAL(Iyr)
      IF (VarSurvFreq.EQ.0) THEN
        Bhat = TBhat
      ELSE
        IF(MOD((Yrcount-1.0),DEPMfrq).EQ.0.0) Bhat = TBhat
      ENDIF
C	  
C     SENSITIVITIES 
C     Hyperstability
      IF (HyperSt.EQ.1.0) THEN
        qind = 1
        qind = HPSa*1/SQRT(B1p(Isim,Iyr-1))
        qind = MAX(qind,1.0)
        Bhat = qind*B1p(Isim,Iyr-1)*EpsB(Iyr-1)
      ENDIF
C
C     SENSITIVITIES
C     Time varying M
      IF(TVM.EQ.1) NatM = Myear(Iyr-1)
C
C     Calculate the harvest guideline
      RBC = 0
      Cpred = 0
      Eff = 0
      ObsEff = 0
C
      RBC = RecBioCatch(Isim,Iyr-1)
      ObsEff = RBC/Bhat
C     Make sure catch does not go below 0
      RBC = MAX(RBC,0.0)
C
C     SENSITIVITIES
C     3. Multiple fleets - Multiply RBC by distribution
      IF(MultFleet.EQ.1.0) THEN
C       Don't forget to still set the minimum to the bait catch
        RBC = MAX(RBC,0.0)
C        RBC = MAX(RBC,2000.0)
C       Multiply by DISTRIBUTION (0.87)
        RBC = RBC * 0.87
      ENDIF
C
C     Make sure the catch does not exceed the biomass
      IF(RBC.GE.0.95*B1p(Isim,Iyr-1)) RBC = B1p(Isim,Iyr-1)*0.95
C
C     Save the F values for the report
      Eff = RBC/B1p(Isim,Iyr-1)
C     In case SSB was used for the assessment
      IF (WhatAssmt.EQ.1.0) Eff = RBC/SSB(Isim,Iyr-1)
      IF(OFLopt.EQ.1) THEN
       Eff = (Emsymax-Emsymin)/(Nsim-1)*(Isim-1) + Emsymin
      ENDIF
C	  
C     Start calculation of mean age of the catch
      AgeC(Isim,Iyr) = 0
C     Solve for the catch
      IF (RBC.LE.0) THEN
       FF = 0
       DO 10000 Age = 0,Nage
        Z(Age) = NatM
10000  CONTINUE
      ELSE
       Fmin = 0
       Fmax = 2.996
C      Use a bisection algorithm to guarantee that the RBC is taken
       DO 11000 II = 1,23
        FF = (Fmin+Fmax)/2.0
        CPred = 0
        ACcount1 = 0
        ACcount2 = 0
        DO 11100 Age = 0,Nage
         Z(Age) = SelxSum(Iyr-1,Age)*FF+NatM
C
C        SENSITIVITIES
C        3. Multiple fleet - Update dynamics
         IF (MultFleet.EQ.1) THEN
          Z(Age) = SelxSum(Iyr-1,Age)*FF + NatM + SelxMC2(Age)*Fmex +
     +             SelxPNW(age)*Fcan
         ENDIF
C
         Cnum = N(Iyr-1,Age)*SelxSum(Iyr-1,Age)*FF/Z(Age)*
     +          (1-exp(-Z(Age)))
         Cpred = Cpred  + WeightSum(Iyr-1,Age)*Cnum
         ACcount1 = ACcount1 + Cnum*Age
         ACcount2 = ACcount2 + Cnum
11100   CONTINUE 
        IF (Cpred.GT.RBC) THEN
         Fmax = FF
        ELSE
         Fmin = FF
        ENDIF
11000  CONTINUE        
      ENDIF 
      F(Isim,Iyr-1) = Eff
      AgeC(Isim,Iyr-1) = ACcount1/ACcount2
      IF(Cpred.LE.0.0001) THEN
        AgeC(Isim,Iyr-1) = 0.000
        F(Isim,Iyr-1) = 0.000
      ENDIF
      Catch(Isim,Iyr-1) = Cpred
C
C     Update the dynamics
      AgeN(Isim,Iyr) = 0
      N(Iyr,Nage) = N(Iyr-1,Nage)*exp(-Z(Nage))
      N(Iyr,Nage) = N(Iyr,Nage) + N(Iyr-1,Nage-1)*exp(-Z(Nage-1))
      ANcount1 = NAge*N(Iyr,Nage)
      ANcount2 = N(Iyr,Nage)
      DO 20000 Age = 1,Nage-1
       N(Iyr,Age) = N(Iyr-1,Age-1)*exp(-Z(Age-1))
       ANcount1 = ANcount1 + N(Iyr,Age)*Age
       ANcount2 = ANcount2 + N(Iyr,Age)
20000 CONTINUE
C      WRITE(*,'(1x,F5.3,1x,F5.3,1x,F5.3,1x,F5.3)') Z(0),Z(1),Z(2),NatM
C
C     Calculate mean age of the population
      AgeN(Isim,Iyr) = ANcount1/ANcount2
      IF(ANcount2.LE.1.0) AgeN(Isim,Iyr) = 0
      Ncheck = 0
C	  
      DO 20010 Age=0,Nage
       Ncheck = Ncheck + N(Iyr,Age)
20010 CONTINUE
C
      IF(Ncheck.LE.0.0001) AgeN(Isim,Iyr-1) = 0.000
C
      IF(Iyr.EQ.Nyear) THEN
       IF(Ncheck.LE.0.0001) AgeN(Isim,Iyr) = 0.000 
      ENDIF
C
C     Compute the SSB and 1+ biomass
      Bio1P = 0
      DO 21000 Age = 1,Nage
       Bio1P = Bio1P + WeightSum(Iyr,Age)*N(Iyr,Age)
21000 CONTINUE
C
      SSBV = 0
      DO 22000 Age = 2,Nage
       SSBV = SSBV + WeightSum(Iyr,Age)*N(Iyr,Age)
22000 CONTINUE
C
      SSB(Isim,Iyr) = SSBV
      B1p(Isim,Iyr) = Bio1P
C
C     Get the recruitment ((EQUATION 2))
C
C     Base case is a Beverton-Holt model,
      RecR = (4*Steepness*R0*SSBV/SSB0)/
     +       ((1-Steepness)+(5*Steepness-1)*SSBV/SSB0)
C	 
C      WRITE(*,'(1x,F11.1,1x,F5.1,1x,F6.1,1x,F6.1)') SSBV, RecR,  
C     +                                       N(Iyr,1),N(Iyr,2)
C      RecR = SSBV/1000 * EXP(AlogRPS + BlogRPS*SSBV/1000 + 
C     +       PhilogRPS*V(Iyr))
CC     the log(R) model, dome-shaped      
C      IF(Rtype.EQ.1) THEN
C       RecR = EXP(Alpha + Beta1*SSBV/1000 + Beta2*(SSBV/1000)**2 +
C     +        Phi*V(Iyr))
C      ENDIF
CC     The log(R) model, flat-top
CC     First, set the upper limit of the recruitment curve
C      SSBflat = -Beta1/(2*Beta2)
C      IF(Rtype.EQ.2) THEN
C       IF((SSBV/1000).GT.SSBflat) THEN
C        RecR = EXP(Alpha+Beta1*SSBflat+Beta2*(SSBflat)**2+
C     +         Phi*V(Iyr))
C       ENDIF
C      ENDIF
CC     A simple Ricker model with no temperature term
C      IF(Rtype.EQ.3) THEN
C       RecR = SSBV*EXP(RickerA + RickerB * SSBV)
C      ENDIF	   
C      IF(Rtype.EQ.4) THEN
C       RecR = SSBV*EXP(RickerA + RickerB * SSBV)
C      ENDIF	 
      RecR = RecR * EXP(Eps(Iyr)-(SigmaR**2)/2.0)
      N(Iyr,0) = RecR
      RecRep(Isim,Iyr) = RecR
C      
      RETURN
600   FORMAT(I4,1x,100(F12.3,1x))      
      END
C
C ===============================================================================================
C
      SUBROUTINE Summary()
C
      IMPLICIT NONE
      INCLUDE "Chile_OM.INC"
C      
C     Local variables
      INTEGER Isim,Year,Iage
C
      DO 10000 Year = 0,Nyear
       WRITE(14,600) Year,(B1p(Isim,Year),SSB(Isim,Year),
     +                     RecRep(Isim,Year),Catch(Isim,Year),
     +                     F(Isim,Year),ObsF(Isim,Year),Vrep(Isim,Year),
     +                     AgeN(Isim,Year),AgeC(Isim,Year),Isim= 1,Nsim)
10000 CONTINUE     
C
      IF(SelxOut.EQ.1.0) THEN
       DO 11000 Year = 0,Nyear
        WRITE(14,700) Year, (SelxSum(Year,Iage),Iage=0,MaxAge)
11000 CONTINUE      
      ENDIF
C
      IF(WeightOut.EQ.1.0) THEN
        WRITE(14,700) Year, (WeightSum(Year,Iage),Iage=0,MaxAge)       
      ENDIF
C
C     These are only for debugging. Comment them once the model works      
C      WRITE(14,*) F(1,2), F(1,3), F(1,4)
C      WRITE(*,*) F(1,2), F(1,3), F(1,4)
C      WRITE(14,*) Niprop
C      WRITE(14,*) N(1,0), N(1,1), N(1,2), N(1,3), N(1,4), N(1,5)
C      WRITE(14,*) SSB(1,0)
C      WRITE(14,*) V(1),V(2),V(3),V(4),V(5),V(6),V(7),V(8),V(9),V(10)
C
C      DO 11000 Year=0,Nyear
C       WRITE(14,*) Eps(Year)
C11000 CONTINUE
C	   
      RETURN
600   FORMAT(1x,I5,1x,1000(F14.2,1x,F14.2,1x,F14.2,1x,F12.2,1x,F6.4,1x,
     +       F6.4,1x,F7.4,1x,F6.4,1x,F6.4))
700   FORMAT(1x,I5,1x,20(F5.3,1x))    
      END
C
C ===============================================================================================
C
      SUBROUTINE CalcG()
C
C This subroutine sets the values for the underlying environmental signal
C
      IMPLICIT NONE
      INCLUDE "Chile_OM.INC"
C
      INTEGER Year
C
      REAL*8 Term1,Term2,Term3,Term4
      REAL*8 TwoPi
      TwoPi = 6.283185
C
      DO 10000 Year=0,Nyear
C      Create G
C       Sine wave (c)
       Term1 = Sin(TwoPi*(Year-Ybar)/Period)
C 
C       Basic square wave (a1)
       IF (MOD(Year,NINT(Period)).LT.Period/2.0) THEN
        Term2 = 1
        ELSE
         Term2 = -1
       ENDIF
C
C       Square wave with shorter good periods (a2)
       IF (MOD(Year,NINT(Period)).LT.Period/3.0) THEN
        Term3 = 1
        ELSE
         Term3 = -1
       ENDIF
C	   
C       Square wave with slowly descending shift from god to bad (b)
       IF (MOD(Year,NINT(Period)).LT.Period/2.0) THEN
        Term4 = 1
        ELSE
         IF (MOD(Year,NINT(Period)).GE.Period/2.0+10) THEN
           Term4 = -1
         ELSE
           Term4 = -0.2*(MOD(Year,NINT(Period))-Period/2.0)+1	   
         ENDIF
       ENDIF
C
       IF (Gtype.EQ.1) G(Year) = Term1
       IF (Gtype.EQ.2) G(Year) = Term2
       IF (Gtype.EQ.3) G(Year) = Term3
       IF (Gtype.EQ.4) G(Year) = Term4
C	   
10000 CONTINUE
C
      END
C
C ===============================================================================================
C
      REAL*8 FUNCTION RecBioCatch(Isim,Iyr)
C
      IMPLICIT NONE
      INCLUDE "Chile_OM.INC"
      
      REAL*8 RBCt, Emsy,IEmsy,Frac,Imovave,ABC
      REAL*8 Yrcount,SurvAge
      INTEGER Isim, Iyr

C     Calculate Emsy
C     By default (Tier=0), set Emsy=Emsymin, i.e. constant harvest rate. 
C     Note this is also Tier 1 option 2
      Emsy = Emsymin 
C
C     Calculate the OFL
      OFL = Emsy * Bhat
C	  
C     Calculate the base RCB
      RBCt = Emsy * Bhat
C	  
C     If Tier=1 set Emsy to a value depending on the age of the survey
      IF (Tier.EQ.1) THEN
        SurvAge = 1.0
        IF  (VarSurvFreq.EQ.1) THEN
          Yrcount = REAL(Iyr)
          SurvAge = MOD((Yrcount),DEPMfrq)
          IF (SurvAge.EQ.0.0) SurvAge = DEPMfrq
          IF (Iyr.EQ.0) SurvAge = 0.0
        ENDIF
        Emsy = - T1psi * SurvAge + T1gamma
        Emsy = MIN(Emsy,T1Max)
        RBCt = Bhat * Emsy
      ENDIF
C     If Tier=2, set the RBC equal to the species/stock max
      IF (Tier.EQ.2) RBCt = T2Max
	  
C     If Tier=3, set the RBC equal to Tier 3 max (by default, 500t)
      IF (Tier.EQ.3) RBCt = T3Max

C     If OFLopt=1, do a range of values for Emsy (Between Emsymin and Emsymax)
C      IF(OFLopt.EQ.1) THEN
C       Emsy = (Emsymax-Emsymin)/(Nsim-1)*(Isim-1) + Emsymin
C      ENDIF
C
C
C
C     If OFLopt=3, use the relationship between Emsy and SST
C      IF (OFLopt.EQ.2) THEN
C       Emsy = - 18.46452 + 3.25209*Imovave - 0.19723*Imovave**2 + 
C     +        0.0041863*Imovave**3
C       IEmsy = Emsy
C       Emsy = MAX(Emsy,0.0)
C       Emsy = MIN(Emsy,OFLmax)
C       OFL = Emsy*Bhat	  
C      ENDIF	  
C     
C     Calculate the ABC (p*=0.40)
      ABC = Pstar*OFL
C	  
C     Cutoff can be looped over 
C      IF (CutOpt.EQ.1) Cutoff = 50000*(Isim-1)
C
C     Calculate the RBC
C     ================
C     By default (Tier=0), set the RBC = Emsy
C      RBCt = OFL
      ObsF(Isim,Iyr) = Emsy
C
C     It can also be set with a cutoff and maxcatch and pre-specified Fraction=FracMin
C
C     Constant F but constrained by Maxcatch, Cutoff and limits to the HCR
C      IF (Tier.EQ.1) THEN
C       Frac = FracMin
C       RBCt = Frac*(Bhat-Cutoff)
C       RBCt = MIN(RBCt,Maxcatch)
C       RBCt = MIN(RBCt,ABC)
C       RBCt = MAX(RBCt,0.0)
C       ObsF(Isim,Iyr) = Frac
C      ENDIF
CC
CC     Or include a temperature term to determine Fraction
C      IF (Tier.EQ.2) THEN
C       Frac = Emsy
C       Frac = MIN(Frac,FracMax)
C       Frac = MAX(Frac,FracMin)
C       RBCt = Frac*(Bhat-Cutoff)
C       RBCt = MIN(RBCt,Maxcatch)
C       RBCt = MIN(RBCt,ABC)
C       RBCt = MAX(RBCt,0.0)
C       ObsF(Isim,Iyr) = Frac
C      ENDIF
C	  
      RecBioCatch = RBCt
C
      RETURN
      END
C
C ===============================================================================================
C
      INCLUDE "COMMON.FOR"