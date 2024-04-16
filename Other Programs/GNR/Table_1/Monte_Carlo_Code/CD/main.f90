PROGRAM MAIN
  USE PROBABILITY
  USE INTERPOL
  USE RANDOM
  USE MINIMIZATION
  USE SIMPLEX
  USE INTEGRATION
  USE GLOBVAR
  USE MPI
  IMPLICIT NONE
  
  INTEGER :: i,j,s,t
  REAL(8) :: theta(1),dist,point(2),it,etat,epst,wtm1,prod,kt1,aux,lo,hi,iti,mti,yti
  REAL(8), POINTER :: fx_ch0(:,:),fx_ch_temp(:)
  INTEGER :: num1,num2,num3,num4
  INTEGER, POINTER :: mylo(:),myhi(:)
  INTEGER :: ierr,myid,numprocs
  
  myid=0
  numprocs=1
  
  CALL MPI_INIT( ierr ) ! Always call it at the begining of the program
  CALL MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr ) ! Always call it at the begining of the program, myid is going to be the id for the processor
  CALL MPI_COMM_SIZE( MPI_COMM_WORLD,numprocs,ierr) ! Always call it at the begining of the program, numprocs is the number of processors that you gave when you ran it     
  WRITE(*,*) "process",myid," of",numprocs," is alive"
  
  ALLOCATE(mylo(0:numprocs-1),myhi(0:numprocs-1))
  
!   IF (myid==0) THEN
!        WRITE(*,*) 'number'
!        READ(*,*) j
!   END IF
 
   d=(/0.05d0,0.075d0,0.1d0,0.125d0,0.15d0/)
  limits(:,1)=(/1.0d0,0.5d0/)
  limits(:,2)=(/6000.0d0,6.0d0/)
  DO j=1,N_d
     CALL MLINEAR_INIT(vt1(j),(/20,20/),limits,(/1,0/))
     CALL MLINEAR_INIT(i_pol(j),(/20,20/),limits,(/1,0/))
     CALL MLINEAR_INIT(mt_pol(j),(/20,20/),limits,(/1,0/))
  END DO
  ALLOCATE(fx_ch(N_d,vt1(1)%lt),fx_ch0(N_d,vt1(1)%lt),i_ch(N_d,vt1(1)%lt),mt_ch(N_d,vt1(1)%lt))
  ALLOCATE(fx_ch_temp(vt1(1)%lt))
  
  CALL GAUSS_HERMITE(ab,we)
  ab=ab*SQRT(2.0d0*vin)
  DO j=1,N_int
     we(j)=we(j)/SQRT(3.1415926535897932384626433832d0)
  END DO
  
  ! For MPI
  num1=vt1(1)%lt/numprocs
  num2=vt1(1)%lt-(num1*numprocs)
  DO i=1,numprocs
     num3=0
     IF (i<=num2) num3=1
     num4=num1+num3
     IF (i==1) THEN
        mylo(i-1)=1
     ELSE
        mylo(i-1)=myhi(i-2)+1
     END IF
     myhi(i-1)=mylo(i-1)+num4-1
  END DO
   
  ! Let's solve the model, for that I need to first guess the value function
  ! with that guess I can do the interpolation as an initial guess
  ! So let's do that
  DO h=1,N_d
      DO i=1,vt1(1)%lt
         fx_ch(h,i)=1.0d0 !Sample_Uniform(0.0d0,1.0d0)
      END DO
      i=0
      theta=20.0d0
      DO
         fx_ch0(h,:)=fx_ch(h,:)
         i=i+1
         
         DO j=mylo(myid),myhi(myid)
            kt=vt1(h)%v(1)%x(vt1(h)%i(j)%i(1))
            wt=vt1(h)%v(2)%x(vt1(h)%i(j)%i(2))
            mt=(LOG(am)+ak*LOG(kt)+LOG(wt)-LOG(Pm/P)) / (1.0d0-am)                 
            mt=EXP(mt)
            theta=20.0d0
            ! Ok so now we have a point in the grid        
!            IF (i<50) THEN
!              CALL NELDER_MEADE(theta,1.0d-7,OBJECTIVE,myid,100000)
!              CALL NELDER_MEADE(theta,1.0d-7,OBJECTIVE,myid,100000)
!              CALL NELDER_MEADE(theta,1.0d-7,OBJECTIVE,myid,100000)
!              CALL NELDER_MEADE(theta,1.0d-7,OBJECTIVE,myid,100000)
!              CALL NELDER_MEADE(theta,1.0d-7,OBJECTIVE,myid,100000)
!              CALL NELDER_MEADE(theta,1.0d-7,OBJECTIVE,myid,100000)
!              CALL NELDER_MEADE(theta,1.0d-7,OBJECTIVE,myid,100000)
!              CALL NELDER_MEADE(theta,1.0d-7,OBJECTIVE,myid,1000000)
!!               CALL BFGS(theta,1.0d-6,1.0d-6,OBJECTIVE,myid,.TRUE.)
!            ELSE IF(i<200) THEN
!              CALL NELDER_MEADE(theta,1.0d-9,OBJECTIVE,myid,100000)
!              CALL NELDER_MEADE(theta,1.0d-9,OBJECTIVE,myid,100000)
!              CALL NELDER_MEADE(theta,1.0d-9,OBJECTIVE,myid,100000)
!              CALL NELDER_MEADE(theta,1.0d-9,OBJECTIVE,myid,100000)
!              CALL NELDER_MEADE(theta,1.0d-7,OBJECTIVE,myid,100000)
!              CALL NELDER_MEADE(theta,1.0d-7,OBJECTIVE,myid,100000)
!              CALL NELDER_MEADE(theta,1.0d-9,OBJECTIVE,myid,100000)
!              CALL NELDER_MEADE(theta,1.0d-9,OBJECTIVE,myid,1000000)
!!               CALL BFGS(theta,1.0d-8,1.0d-8,OBJECTIVE,myid,.FALSE.)
!            ELSE
              CALL NELDER_MEADE(theta,1.0d-11,OBJECTIVE,myid,100000)
              CALL NELDER_MEADE(theta,1.0d-11,OBJECTIVE,myid,100000)
              CALL NELDER_MEADE(theta,1.0d-7,OBJECTIVE,myid,100000)
              CALL NELDER_MEADE(theta,1.0d-7,OBJECTIVE,myid,100000)
              CALL NELDER_MEADE(theta,1.0d-11,OBJECTIVE,myid,100000)
              CALL NELDER_MEADE(theta,1.0d-11,OBJECTIVE,myid,100000)
              CALL NELDER_MEADE(theta,1.0d-11,OBJECTIVE,myid,100000)
              CALL NELDER_MEADE(theta,1.0d-11,OBJECTIVE,myid,1000000)
!               CALL BFGS(theta,1.0d-10,1.0d-10,OBJECTIVE,myid,.FALSE.)
!            END IF
            fx_ch_temp(j)=-OBJECTIVE(theta)
            lo=MAX(0.0d0,limits(1,1)-(1.0d0-d(h))*kt)
            hi=limits(1,2)-(1.0d0-d(h))*kt
            hi=hi-lo-1.0d-10
            i_ch(h,j)=lo + (hi/(1.0d0+ (1.4d0**theta(1))))
            mt_ch(h,j)=mt
         END DO
         DO j=0,numprocs-1
            num1=myhi(j)-mylo(j)+1
            CALL MPI_BCAST(fx_ch_temp(mylo(j):myhi(j)),num1,MPI_REAL8,j,MPI_COMM_WORLD,ierr)
            CALL MPI_BCAST(i_ch(h,mylo(j):myhi(j)),num1,MPI_REAL8,j,MPI_COMM_WORLD,ierr)
            CALL MPI_BCAST(mt_ch(h,mylo(j):myhi(j)),num1,MPI_REAL8,j,MPI_COMM_WORLD,ierr)
         END DO
         fx_ch(h,:)=fx_ch_temp
         
         IF (i<50) THEN
            DO s=1,100
               DO j=1,vt1(h)%lt
                  kt=vt1(h)%v(1)%x(vt1(h)%i(j)%i(1))
                  wt=vt1(h)%v(2)%x(vt1(h)%i(j)%i(2))
                  point=(/kt,wt/)
                  it=MLINEAR_INTER(i_pol(h),point,i_ch(h,:))
                  mt=MLINEAR_INTER(mt_pol(h),point,mt_ch(h,:))
!                  mt=(LOG(am)+ak*LOG(kt)+LOG(wt)-LOG(Pm/P)) / (1.0d0-am)                 
 !                 mt=EXP(mt)
                  theta=(/it/)
                  fx_ch_temp(j)=-OBJECTIVE_LEVEL(theta)              
               END DO
               fx_ch(h,:)=fx_ch_temp
            END DO
         END IF
         
         dist=SUM((fx_ch(h,:)-fx_ch0(h,:))**2.0d0)
         CALL FLUSH(6)
         IF (myid==0) WRITE(*,*) h,i,dist,MAXVAL(abs(fx_ch(h,:)-fx_ch0(h,:))),0.00001d0*SQRT(SUM(fx_ch(h,:)*fx_ch(h,:)))
         
         IF (ANY(abs(fx_ch(h,:)-fx_ch0(h,:))>0.00001d0* SQRT(SUM(fx_ch(h,:)*fx_ch(h,:))) )) THEN
            j=1
         ELSE
            EXIT
         END IF
         
      END DO
   END DO
   
   IF (myid==0) THEN
      OPEN(523,file="data.out")
      ! Ok so now I have the solution to the Bellman Problem, let's create firms
      DO i=1,N_firms
         h=Sample_Uniform(1.0d0,DBLE(N_d+1))
         kt1=Sample_Uniform(11.0d0,400.0d0)
         wtm1=Sample_Uniform(1.0d0,3.0d0)
         DO t=1,N_time
            kt=kt1
            etat=Sample_Normal(0.0d0,vin)
            epst=Sample_Normal(0.0d0,0.04d0)
            wt=d0+d1*LOG(wtm1) + etat
            wt=EXP(wt)
            hi=limits(2,2)
            IF (wt>hi) wt=hi
            point=(/kt,wt/)
            iti=MLINEAR_INTER(i_pol(h),point,i_ch(h,:))
            mti=MLINEAR_INTER(mt_pol(h),point,mt_ch(h,:))
!            mt=(LOG(am)+ak*LOG(kt)+LOG(wt)-LOG(Pm/P)) / (1.0d0-am)
!            mt=EXP(mt)
            theta=20.0d0
!            WRITE(*,*) i,t,etat,epst,iti,mti,mt
            mt=mti
            CALL NELDER_MEADE(theta,1.0d-8,OBJECTIVE,0,10000)
            CALL NELDER_MEADE(theta,1.0d-8,OBJECTIVE,0,10000)
            !CALL BFGS(theta,1.0d-10,1.0d-10,OBJECTIVE,0,.TRUE.)
            lo=MAX(0.0d0,limits(1,1)-(1.0d0-d(h))*kt)
            hi=limits(1,2)-(1.0d0-d(h))*kt
            hi=hi-lo-1.0d-10
            it=lo + (hi/(1.0d0+ (1.4d0**theta(1))))
            prod=(kt**ak)*(mt**am)*wt*exp(epst)
            yti=(kt**ak)*(mti**am)*wt*exp(epst)
            ! Now prepare for next period
            kt1=(1.0d0-d(h))*kt + it
            wtm1=wt
            s=t-150
            IF (t>170) WRITE(523,'(3I20,10F32.16)') i,t,h,kt,mt,prod!,wt,etat,epst,it,iti,mti,yti            
         END DO
      END DO
   END IF
   
   CALL MPI_FINALIZE(ierr)
   
 END PROGRAM MAIN
 
 
 
