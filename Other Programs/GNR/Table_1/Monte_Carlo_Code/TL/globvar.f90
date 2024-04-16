MODULE GLOBVAR
  USE INTERPOL
  IMPLICIT NONE
  INTEGER, PARAMETER :: N_time=200,N_firms=50000
  REAL(8), PARAMETER :: Pm=1.0d0,P=1.0d0,Pi=8.0d0,disc=0.985d0
  REAL(8), PARAMETER :: am=0.65d0,ak=0.25d0,amm=0.015d0,akk=0.015d0,akm=-0.032d0!amm=0.02d0,akk=0.02d0,akm=-0.032d0
  REAL(8), PARAMETER :: d0=0.2d0,d1=0.8d0,vin=0.04d0
  INTEGER, PARAMETER :: N_int=8
  INTEGER, PARAMETER :: N_d=5
  REAL(8) :: d(N_d)
  TYPE(MLINEARSTRUC) :: vt1(N_d),i_pol(N_d),mt_pol(N_d)
  REAL(8) :: kt,wt,mt
  REAL(8) :: we(N_int),ab(N_int)
  INTEGER :: h
  REAL(8), POINTER :: fx_ch(:,:),i_ch(:,:),mt_ch(:,:)
  REAL(8) :: limits(2,2)


CONTAINS

  REAL(8) FUNCTION OBJECTIVE(theta)
    USE INTERPOL
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: theta(:)
    REAL(8) :: kt1,wt1,point(2),it,lkt,lmt
    REAL(8) :: EV,mw,lo,hi,arg
    INTEGER :: j
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
!    wt1=EXP(mw)
       lo=limits(2,1)
       IF (wt1<=lo) wt1=lo+1.0d-10
       hi=limits(2,2)
       IF (wt1>=hi) wt1=hi-1.0d-10
       point(2)=wt1
       EV=EV+ MLINEAR_INTER(vt1(h),point,fx_ch(h,:))*we(j)
    END DO
    lkt=LOG(kt)
    lmt=LOG(mt)
    arg=(ak*lkt)+(am*lmt)+(akk*lkt*lkt)+(amm*lmt*lmt)+(akm*lkt*lmt)
    arg=EXP(arg)
    mw=P*arg*wt
    mw=mw - Pm*mt - Pi*it
    OBJECTIVE=-mw-disc*EV
  END FUNCTION OBJECTIVE

  REAL(8) FUNCTION OBJECTIVE_LEVEL(theta)
    USE INTERPOL
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: theta(:)
    REAL(8) :: kt1,wt1,point(2),it,lkt,lmt
    REAL(8) :: EV,mw,lo,hi,arg
    INTEGER :: j
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
!    wt1=EXP(mw)
       lo=limits(2,1)
       IF (wt1<=lo) wt1=lo+1.0d-10
       hi=limits(2,2)
       IF (wt1>=hi) wt1=hi-1.0d-10
       point(2)=wt1
       EV=EV+ MLINEAR_INTER(vt1(h),point,fx_ch(h,:))*we(j)
    END DO
    lkt=LOG(kt)
    lmt=LOG(mt)
    arg=(ak*lkt)+(am*lmt)+(akk*lkt*lkt)+(amm*lmt*lmt)+(akm*lkt*lmt)
    arg=EXP(arg)
    mw=P*arg*wt
    mw=mw - Pm*mt - Pi*it

    
    OBJECTIVE_LEVEL=-mw-disc*EV
  END FUNCTION OBJECTIVE_LEVEL
  
END MODULE GLOBVAR
