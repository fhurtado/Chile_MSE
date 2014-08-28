C     Australian SPF simulation code. 
C     This program reads the summary.out file, orders the data, and calculates means, medians and quantiles
C     Written by F. Hurtado-Ferro and A.E. Punt
C     Please reference the authors if you use this code.
C
      PROGRAM OutSort
      IMPLICIT NONE
	  
      INCLUDE "Sort.inc"
  
      INTEGER Yr 

      WRITE(*,*) "Start calculation of performance indices"
      CALL ReadData()	
	  
      OPEN(UNIT=14,FILE="PerfInd.OUT")
	  
      WRITE(*,*) "Dimensions"
      WRITE(*,*) "Nyear", Nyear, "Nsim", Nsim
	  
      WRITE(*,*) "End reading, start summarizing"
      CALL SummarizeData()	  
      WRITE(*,*) "End summarizing, start sorting"

C     Sort the big objects
      CALL SORT(Biomass,(Nyear+1)*Nsim)
      CALL SORT(SSB,(Nyear+1)*Nsim)
      CALL SORT(Catches,(Nyear+1)*Nsim)	  
C     Sort the small objects
      CALL SORT(Bio5,5*Nsim)
      CALL SORT(SSB5,5*Nsim)
      CALL SORT(Catch5,5*Nsim)	
      WRITE(*,*) "End sorting, write performance indices"

      CALL Summary()
      CLOSE(14)
	  
      END PROGRAM OutSort
C
C -------------------------------------------------------------------------
C    
C     Quick sort routine, written by Andre E. Punt 
      SUBROUTINE SORT(X,M)

C     USE A QUICK-SORT TO SORT ALL THE DATA

      IMPLICIT NONE
      INCLUDE "Sort.inc"

      REAL*8 X((Nyear+1)*Nsim),ST1((Nyear+1)*Nsim),MID
      INTEGER*4 M,LEFT((Nyear+1)*Nsim),RIGHT((Nyear+1)*Nsim),STKLEN	  
      INTEGER*4 LEFTS,RIGHTS,LS,RS,IC

C     Check for Daft call
      IF (M.LT.2) RETURN

C     Set up initial conditions
      LEFT(1) = 1
      RIGHT(1) = M
      STKLEN = 1

99    IF (STKLEN.EQ.0) GOTO 100

C     Set up the Pointers for this run
      MID = x(LEFT(STKLEN))
      LEFTS = LEFT(STKLEN)
      RIGHTS = RIGHT(STKLEN)
      LS = LEFT(STKLEN)
      RS = RIGHT(STKLEN)
                                      
C     Do a one-level sort
      DO 10 IC = LEFT(STKLEN)+1,RIGHT(STKLEN)

C      Check whether the current is less than the middle
       IF (X(IC).GT.MID) THEN
         ST1(RIGHTS) = X(IC)
         RIGHTS = RIGHTS - 1
       ELSE
         ST1(LEFTS) = X(IC)
         LEFTS = LEFTS + 1
       ENDIF
10    CONTINUE

C     Store the middle value
      ST1(LEFTS) = x(LEFT(STKLEN))

C     Replace the data
      DO 11 IC = LEFT(STKLEN),RIGHT(STKLEN)
       x(IC) = ST1(IC)
11    CONTINUE
      STKLEN = STKLEN - 1
        
C     update right pointer
      IF ((LEFTS-LS).GT.1) THEN
        STKLEN = STKLEN + 1
        LEFT(STKLEN) = LS
        RIGHT(STKLEN) = LEFTS - 1
      ENDIF
        
C     update left pointer
      IF ((RS-RIGHTS).GT.1) THEN
        STKLEN = STKLEN + 1
        LEFT(STKLEN) = RIGHTS + 1
        RIGHT(STKLEN) = RS
      ENDIF

      GOTO 99
100   CONTINUE

      RETURN
      END
C
C -------------------------------------------------------------------------
C 
      SUBROUTINE ReadData

      IMPLICIT NONE	
      INCLUDE "Sort.INC"	  
      INTEGER Yr,Col   
	  
C	  Get dimensions from the spec file
      OPEN (UNIT=9, FILE="Chile_OM.spec")
      READ(9,*)
      READ(9,*) Nyear
      READ(9,*)
      READ(9,*) Nsim
      CLOSE(9)
      ADim = Nsim*9

C     Read the summary and save it into a big thing	 
      WRITE(*,*) "Start reading summary.out" 
      OPEN (UNIT=11, FILE="SUMMARY.OUT")
      DO 10000 Yr=0,Nyear
        READ(11,*) (Tres(Col),Col=0,ADim)
        DO 11000 Col=1,Adim	  
          Resmat(Yr,Col) = Tres(Col)	
11000   CONTINUE	
10000 CONTINUE  
      CLOSE(11)
	  
      RETURN
      END
C
C -------------------------------------------------------------------------
C 
      SUBROUTINE SummarizeData	  
	  
      IMPLICIT NONE	
      INCLUDE "Sort.INC"	  
      INTEGER Yr,Col,Ind
C	  
C     Organize the full time series
      Ind=0
      DO 12000 Col=1,Nsim	  
       DO 12100 Yr=0,Nyear
        Biomass(Ind) = Resmat(Yr,(Col-1)*9+1)
        SSB(Ind) = Resmat(Yr,(Col-1)*9+2)
        Catches(Ind) = Resmat(Yr,(Col-1)*9+4)		
        Ind=Ind+1
12100  CONTINUE	  
12000 CONTINUE
C      WRITE(*,*) "Ind is ", Ind
C     
C     Organize the last 5 years of the time series
      Ind=0 
      DO 12300 Col=1,Nsim	  
       DO 12400 Yr=Nyear-6,Nyear-1
        Bio5(Ind) = Resmat(Yr,(Col-1)*9+1)
        SSB5(Ind) = Resmat(Yr,(Col-1)*9+2)
        Catch5(Ind) = Resmat(Yr,(Col-1)*9+4)		
        Ind=Ind+1
12400  CONTINUE	  
12300 CONTINUE
C      WRITE(*,*) "Ind is ", Ind
C      WRITE(*,*) "Finished"
      END

C
C -------------------------------------------------------------------------
C 
      SUBROUTINE Summary()

      IMPLICIT NONE
      INCLUDE "Sort.inc"
      INTEGER Med,Q05,Q25,Q75,Q95
      INTEGER Med5,Q505,Q525,Q575,Q595
      INTEGER Perc,tPerc,tPerc5
      Real*8 Mean,StdDev
      
      Med = (Nyear+1)*Nsim*0.50
      Q05 = (Nyear+1)*Nsim*0.05
      Q25 = (Nyear+1)*Nsim*0.25
      Q75 = (Nyear+1)*Nsim*0.75
      Q95 = (Nyear+1)*Nsim*0.95
	  
      Med5 = 5*Nsim*0.50
      Q505 = 5*Nsim*0.05
      Q525 = 5*Nsim*0.25
      Q575 = 5*Nsim*0.75
      Q595 = 5*Nsim*0.95

      WRITE(14,600) "Mean",Mean(Biomass,(Nyear+1)*Nsim),
     +              Mean(SSB,(Nyear+1)*Nsim),
     +              Mean(Catches,(Nyear+1)*Nsim),
     +              Mean(Bio5,5*Nsim),
     +              Mean(SSB5,5*Nsim),
     +              Mean(Catch5,5*Nsim)
      WRITE(14,600) "StdDev",StdDev(Biomass,(Nyear+1)*Nsim),
     +              StdDev(SSB,(Nyear+1)*Nsim),
     +              StdDev(Catches,(Nyear+1)*Nsim),
     +              StdDev(Bio5,5*Nsim),
     +              StdDev(SSB5,5*Nsim),
     +              StdDev(Catch5,5*Nsim)
      WRITE(14,600) "Q05",Biomass(Q05),SSB(Q05),Catches(Q05),
     +              Bio5(Q505),SSB5(Q505),Catch5(Q505)
      WRITE(14,600) "Q25",Biomass(Q25),SSB(Q25),Catches(Q25),
     +              Bio5(Q525),SSB5(Q525),Catch5(Q525)
      WRITE(14,600) "Median",Biomass(Med),SSB(Med),Catches(Med),
     +              Bio5(Med5),SSB5(Med5),Catch5(Med5)
      WRITE(14,600) "Q75",Biomass(Q75),SSB(Q75),Catches(Q75),
     +              Bio5(Q575),SSB5(Q575),Catch5(Q575)
      WRITE(14,600) "Q95",Biomass(Q95),SSB(Q95),Catches(Q95),
     +              Bio5(Q595),SSB5(Q595),Catch5(Q595)
      WRITE(*,*) "Stuff"	  
      DO 20000 Perc=1,100
        tPerc = (Nyear+1)*Nsim*REAL(Perc)/100
        tPerc5 = 5*Nsim*REAL(Perc)/100
        WRITE(14,700) Perc,Biomass(tPerc),SSB(tPerc),Catches(tPerc),
     +                Bio5(tPerc5),SSB5(tPerc5),Catch5(tPerc5)       
20000 CONTINUE	 
	   
      RETURN
600   FORMAT(1x,A6,1x,F14.2,1x,F14.2,1x,F14.2,1x,F14.2,1x,
     +       F14.2,1x,F14.2,1x)
700   FORMAT(1x,I6,1x,F14.2,1x,F14.2,1x,F14.2,1x,F14.2,1x,
     +       F14.2,1x,F14.2,1x)
      END
C
C -------------------------------------------------------------------------
C 
      REAL*8 FUNCTION Mean(X,M)	  
	  
      IMPLICIT NONE
      INCLUDE "Sort.inc"
      REAL*8 X((Nyear+1)*Nsim),Tsum
      INTEGER*4 M, Ind

      Tsum=0
      DO 30000 Ind=0,M
       Tsum=Tsum+X(Ind)
30000 CONTINUE

      Mean=Tsum/M
	  
      RETURN
      END
	  
C
C -------------------------------------------------------------------------
C 

      REAL*8 FUNCTION StdDev(X,M)
	  
      IMPLICIT NONE
      INCLUDE "Sort.inc"
      REAL*8 X((Nyear+1)*Nsim),Tsum,Tmean,SSQ
      INTEGER*4 M, Ind	  
	  
      Tsum=0
      DO 40000 Ind=0,M
       Tsum=Tsum+X(Ind)
40000 CONTINUE	  

      Tmean=Tsum/M

      SSQ = 0
      DO 40001 Ind=0,M
       SSQ=SSQ+(X(Ind)-Tmean)**2
40001 CONTINUE

      StdDev = SQRT(SSQ/(M-1))
	  
      RETURN
      END
	  