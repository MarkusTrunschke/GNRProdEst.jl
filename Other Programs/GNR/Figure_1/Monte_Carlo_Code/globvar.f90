MODULE GLOBVAR
  USE INTERPOL
  IMPLICIT NONE
  REAL(8), PARAMETER :: b=0.0d0
  INTEGER, PARAMETER :: N_time=200,N_firms=500
  REAL(8), PARAMETER :: EPm=1.0d0,P=1.0d0,Pi=8.0d0
  REAL(8), PARAMETER :: am=0.65d0,ak=0.25d0,disc=0.985d0
  REAL(8), PARAMETER :: d0=0.2d0,d1=0.8d0,vin=0.04d0
  INTEGER, PARAMETER :: N_int=8
  INTEGER, PARAMETER :: N_d=5
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!! THIS IS WHAT NEEDS TO BE CHANGED TO GENERATE THE FOUR LEVELS OF TIME SERIES VARIATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!! In particular var_pr should be set to 0.00005d0 for half_base, 0.0001d0 for base,     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  !!!!!!!!!!!!!!!!!! 0.0002d0 for two_base and 0.001d0 for ten_base. Change the name of the data in main.f90 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL(8), PARAMETER :: per_pr=0.6d0,var_pr=0.0001d0 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL(8) :: d(N_d)
  TYPE(MLINEARSTRUC) :: vt1(N_d),i_pol(N_d)
  REAL(8) :: kt,pt,wt,mtm1,lrpt,mt
  REAL(8) :: we(N_int),ab(N_int),ab_p(N_int)
  INTEGER :: h
  REAL(8), POINTER :: fx_ch(:,:),i_ch(:,:),mt_ch(:,:)
  REAL(8) :: limits(3,2)


CONTAINS

  REAL(8) FUNCTION OBJECTIVE(theta)
    USE INTERPOL
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: theta(:)
    REAL(8) :: kt1,wt1,point(3),it
    REAL(8) :: deltam,EV,mw,lo,hi,Pm,lpm
    INTEGER :: j,j2
    lo=MAX(0.0d0,limits(1,1)-(1.0d0-d(h))*kt)
    hi=limits(1,2)-(1.0d0-d(h))*kt
    hi=hi-lo - 1.0d-10
    it=lo + (hi/(1.0d0+ (1.4d0**theta(1))))
    kt1=(1.0d0-d(h))*kt+it
!    IF (kt1<=lo) kt1=lo+0.00001d0
!    IF (kt1>=hi) kt1=hi-0.00001d0
    point(1)=kt1
    mw=d0+d1*LOG(wt)    
    EV=0.0d0
    DO j=1,N_int
       wt1=EXP(mw+ab(j))
       lo=limits(3,1)
       IF (wt1<=lo) wt1=lo+1.0d-10
       hi=limits(3,2)
       IF (wt1>=hi) wt1=hi-1.0d-10
       point(3)=wt1
       DO j2=1,N_int
          lpm=LOG(Epm)+per_pr*lrpt+ab_p(j2)
          lo=limits(2,1)
          IF (lpm<=lo) lpm=lo+1.0d-10
          hi=limits(2,2)
          IF (lpm>=hi) lpm=hi-1.0d-10
          point(2)=lpm
          EV=EV+ MLINEAR_INTER(vt1(h),point,fx_ch(h,:))*we(j)*we(j2)
       END DO
    END DO
    Pm=EXP(lrpt)    
    mw=P*(kt**ak)*(mt**am)*wt
    mw=mw - 0.5d0*b*deltam - Pm*mt - Pi*it
    OBJECTIVE=-mw-disc*EV
  END FUNCTION OBJECTIVE

  REAL(8) FUNCTION OBJECTIVE_LEVEL(theta)
    USE INTERPOL
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: theta(:)
    REAL(8) :: kt1,wt1,point(3),it
    REAL(8) :: deltam,EV,mw,lo,hi,Pm,lpm
    INTEGER :: j,j2
    it=theta(1)
    lo=limits(1,1)
    hi=limits(1,2)
    kt1=(1.0d0-d(h))*kt+it
    IF (kt1<=lo) kt1=lo+1.0d-10
    IF (kt1>=hi) kt1=hi-1.0d-10
    point(1)=kt1
    mw=d0+d1*LOG(wt)    
    EV=0.0d0
    DO j=1,N_int
       wt1=EXP(mw+ab(j))
       lo=limits(3,1)
       IF (wt1<=lo) wt1=lo+1.0d-10
       hi=limits(3,2)
       IF (wt1>=hi) wt1=hi-1.0d-10
       point(3)=wt1
       DO j2=1,N_int
          lpm=LOG(Epm)+per_pr*lrpt+ab_p(j2)
          lo=limits(2,1)
          IF (lpm<=lo) lpm=lo+1.0d-10
          hi=limits(2,2)
          IF (lpm>=hi) lpm=hi-1.0d-10
          point(2)=lpm
          EV=EV+ MLINEAR_INTER(vt1(h),point,fx_ch(h,:))*we(j)*we(j2)
       END DO
    END DO
    mw=P*(kt**ak)*(mt**am)*wt
    Pm=EXP(lrpt)
    mw=mw - 0.5d0*b*deltam - Pm*mt - Pi*it
    OBJECTIVE_LEVEL=-mw-disc*EV
  END FUNCTION OBJECTIVE_LEVEL
  
END MODULE GLOBVAR
