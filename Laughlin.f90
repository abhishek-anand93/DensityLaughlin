  !Calculates number density expectation over ground state Laughlin w.f. for a given filling factor
  !Plots a 2D polar profile of density over circular sample 
  !Uses Monte-Carlo method for integral calculation
  !Importance sampling using Metropolis algorithm

PROGRAM MONTE_CARLO
  
   IMPLICIT NONE
   
   REAL, PARAMETER:: pi = 3.1415927
   INTEGER, PARAMETER:: THERM = 1000
   INTEGER*8, PARAMETER:: ITER = 100000
   INTEGER, PARAMETER:: NUMBINS = 50
   
   !===================================================
   
   INTEGER:: N,filling,Nbin,stepsize1,trial
   INTEGER:: val(8),k
   INTEGER, DIMENSION(:), ALLOCATABLE :: seed
   INTEGER:: test,i,iacc,j,i1,i2
   INTEGER*8:: iter1
   REAL*8:: stepsize,radius,gaussian,d1,d2,ratio,Rmax,ran,r1
   REAL*8:: dtheta,arg
   
   !====================================================
   COMPLEX*16::additionalpiece,ap1,hole
   COMPLEX*16, DIMENSION(:), ALLOCATABLE:: zo,zn,swap
   REAL*8, DIMENSION(:), ALLOCATABLE::r,theta,onedim
   REAL*8, DIMENSION(:,:), ALLOCATABLE:: density,sdensity
   CHARACTER(1000)::file1,file2

   CALL DATE_AND_TIME(VALUES=val)
   CALL RANDOM_SEED(size=k)
   ALLOCATE(seed(k))
   seed(:)=val(8)
   CALL RANDOM_SEED(PUT=seed)

10 PRINT*,"Enter the value of N, filling, stepsize1 and trial :"
   READ*,N,filling,stepsize1,trial

   !========================================================================================

   ALLOCATE(density(NUMBINS,NUMBINS),STAT=test)
   ALLOCATE(sdensity(NUMBINS,NUMBINS),STAT=test)
   ALLOCATE(zo(N),STAT=test)
   ALLOCATE(zn(N),STAT=test)
   ALLOCATE(swap(N),STAT=test)
   ALLOCATE(onedim(NUMBINS),STAT=test)
   ALLOCATE(r(N),STAT=test)
   ALLOCATE(theta(N),STAT=test)

   WRITE(file1,50)"density_N",N,"_filling",filling,"_trial",trial
50 FORMAT(A9,I2,A8,I1,A6,I1)
   WRITE(file2,60)"init_",N,"_",filling,"_",trial
60 FORMAT(A5,I2,A1,I1,A1,I1)
   
   radius = sqrt(2.0*N*filling)
   Rmax = radius*2
  
   stepsize = stepsize1/10000.0

   CALL RANDOM_NUMBER(r)
   CALL RANDOM_NUMBER(theta)
  
  r = radius*r
  theta = 2*pi*theta
  
  r1 = radius
  dtheta = 2*pi*0.5
  
  hole = cmplx( r1*cos(dtheta), r1*sin(dtheta))
  
  DO i=1,N
     zo(i) = cmplx( r(i)*cos(theta(i)), r(i)*sin(theta(i)) )
  ENDDO

  iacc = 0

  DO iter1=1,THERM
     gaussian = 0
     
     CALL RANDOM_NUMBER(theta)
     theta = 2*pi*theta
     
     DO i=1,N
        zn(i) = zo(i) + cmplx(stepsize*cos(theta(i)),stepsize*sin(theta(i)))
        d1 = abs(zo(i))
        d2 = abs(zn(i))
        gaussian = gaussian + (d1*d1 - d2*d2)
     ENDDO
     
     ratio = 1
     additionalpiece = 1
     ap1 = 1
     
     DO i=1,N
        IF (i<N) THEN
           DO j=i+1,N
              ratio = ratio*( abs(zn(i)-zn(j)) / abs(zo(i)-zo(j)) )
           ENDDO
        ENDIF
        additionalpiece = additionalpiece*(zn(i)/zo(i))
        ap1 = ap1*(abs(zn(i)-hole)/abs(zo(i)-hole))
     ENDDO
     
     ratio = (ratio**(2*filling))*exp(gaussian/2.0)*(abs(additionalpiece)**2)*(ap1**2)

     CALL RANDOM_NUMBER(ran)
     IF ((ratio > 1) .OR. (ratio > ran)) THEN
        swap = zn
        zn = zo
        zo = swap
        iacc = iacc + 1
     ENDIF
  ENDDO
  
  PRINT*,"Acceptance ratio is :",real(iacc)/THERM

  r1 = Rmax*Rmax/NUMBINS
  dtheta = 2*pi/NUMBINS

  DO i=1,NUMBINS
     DO j=1,NUMBINS
        density(i,j) = 0.0
        sdensity(i,j) = 0.0
     ENDDO
  ENDDO
  
  
  DO i=1,N
     
     d1 = ABS(zo(i))
     i1 = INT( (d1*d1)/(r1) ) + 1

     arg = ATAN2(AIMAG(zo(i)),REAL(zo(i)))
     
     IF (arg < 0) THEN
        arg = 2*pi + arg
     ENDIF
     
     i2 = INT(arg/dtheta) + 1

     IF (i1<=NUMBINS .AND. i2<=NUMBINS) THEN
        sdensity(i1,i2) = sdensity(i1,i2) + 1
     ENDIF
  ENDDO
  
  OPEN(UNIT=100,FILE=file1,STATUS='NEW')
  OPEN(UNIT=200,FILE=file2,status='new')
  !================================================================
  
  DO iter1=1,ITER
     gaussian = 0
     CALL RANDOM_NUMBER(theta)
     theta = 2*pi*theta
     !------------------------------------------------------------------
     DO i=1,N
        zn(i) = zo(i) + cmplx(stepsize*cos(theta(i)),stepsize*sin(theta(i)))
        d1 = abs(zo(i))
        d2 = abs(zn(i))
        gaussian = gaussian + d1*d1 - d2*d2
     ENDDO
     !-------------------------------------------------------------------
     ratio = 1
     additionalpiece = 1
     ap1 = 1
     !-------------------------------------------------------------------
     DO i=1,N
        IF (i<N) THEN
           DO j=i+1,N
                ratio = ratio*abs(zn(i)-zn(j))/abs(zo(i)-zo(j))
             ENDDO
          ENDIF
          additionalpiece = additionalpiece*(zn(i)/zo(i))
          ap1 = ap1*(abs(zn(i)-hole)/abs(zo(i)-hole))
       ENDDO
       !------------------------------------------------------------------
       ratio = (ratio**(2*filling))*exp(gaussian/2.0)*(abs(additionalpiece)**2)*(ap1**2)
     
       
     CALL RANDOM_NUMBER(ran)
     IF ((ratio > 1) .OR. (ratio > ran)) THEN
        swap = zn
        zn = zo
        zo = swap
        
        DO i=1,NUMBINS
           DO j=1,NUMBINS
           sdensity(i,j)=0.0
           ENDDO
        ENDDO
        
        DO i=1,N
           d1 = ABS(zo(i))
           i1 = INT(d1*d1/r1) + 1

           arg = ATAN2(AIMAG(zo(i)),REAL(zo(i)))
           
           IF (arg < 0) THEN
              arg = 2*pi + arg
           ENDIF
           
           i2 = INT(arg/dtheta) + 1

           IF (i1<=NUMBINS .AND. i2<=NUMBINS) THEN
              sdensity(i1,i2) = sdensity(i1,i2) + 1
           ENDIF
        ENDDO
        
     ENDIF

     DO i=1,NUMBINS
        DO j=1,NUMBINS
           density(i,j) = density(i,j) + sdensity(i,j)
        ENDDO
     ENDDO
     
  ENDDO

  DO i=1,NUMBINS
     onedim(i) = 0
  ENDDO
  
  
70  WRITE(200,'(I3)')NUMBINS
  WRITE(200,'(F8.6)')r1
  WRITE(200,'(F8.6)')dtheta
  WRITE(200,'(F8.5)')Rmax
  
  DO i=1,NUMBINS
     DO j =1,NUMBINS
        onedim(j) = density(i,j)
        IF (j==NUMBINS) then
           write(100,'(E9.3,a2)')onedim(j)/ITER,' '
        else  
           write(100,'(E9.3,a2)',advance='no')onedim(j)/ITER,' '
        endif
     ENDDO
     !write(100,'(e9.3)')sum(onedim)/ITER
  ENDDO
  
  PRINT*,sum(density/ITER)
  
  
  CLOSE(100)
  close(200)
  !========================================================================

  
  END PROGRAM MONTE_CARLO
  
  
  
