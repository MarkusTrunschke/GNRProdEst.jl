MODULE SIMPLEX
  IMPLICIT NONE

CONTAINS
  SUBROUTINE Nelder_Meade(x,ftol,func,id,itmax)
    USE nrutil, ONLY : assert_eq,imaxloc,iminloc,nrerror,swap
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: ftol
    INTEGER, INTENT(IN) :: id
    REAL(8), DIMENSION(:), INTENT(INOUT) :: x
    INTERFACE
       FUNCTION func(x)
         IMPLICIT NONE
         REAL(8), DIMENSION(:), INTENT(IN) :: x
         REAL(8) :: func
       END FUNCTION func
    END INTERFACE
    INTEGER, INTENT(IN) :: itmax
    !INTEGER, PARAMETER :: ITMAX=100000
    REAL(8), PARAMETER :: TINY=1.0d-12
    INTEGER :: ihi,ndim,i,j,iter
    REAL(8) :: psum(SIZE(x)),p(SIZE(x)+1,SIZE(x)),y(SIZE(x)+1),pb(SIZE(x)),ver
    !Initialize the simplex
    DO i=2,SIZE(x)+1
       DO j=1,SIZE(x)
          IF (i-1==j) THEN
             ver=ABS(x(j))*0.75d0
             IF (ver<1.0d-5) ver=0.75d0
!             ver=1.0d0
             pb(j) = x(j) + ver
             p(i,j) = x(j) + ver
          ELSE
             pb(j) = x(j)
             p(i,j) = x(j)                
          END IF
       END DO
       y(i)=FUNC(pb)
    END DO
    DO j=1,SIZE(x)
       pb(j)=x(j)
       p(1,j)=x(j)
    END DO
    y(1)=FUNC(pb)    
    call amoeba_private
  CONTAINS

    SUBROUTINE amoeba_private
      IMPLICIT NONE
      INTEGER :: j,i,ilo,inhi
      REAL(8) :: rtol,ysave,ytry,ytmp,tim
      INTEGER :: tin,tout,hz
      ndim=assert_eq(size(p,2),size(p,1)-1,size(y)-1,'amoeba')
      iter=0
      psum(:)=sum(p(:,:),dim=1)
      DO
         ilo=iminloc(y(:))
         x=p(ilo,:)
         ihi=imaxloc(y(:))
         ytmp=y(ihi)
         y(ihi)=y(ilo)
         inhi=imaxloc(y(:))
         y(ihi)=ytmp
         CALL SYSTEM_CLOCK(count_rate=hz) 
         CALL SYSTEM_CLOCK(count=tin)        
         rtol=2.0d0*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo))+TINY)
         if (rtol < ftol) then
            call swap(y(1),y(ilo))
            call swap(p(1,:),p(ilo,:))
            EXIT
         end if
         IF (iter >= ITMAX) EXIT
         ytry=amotry(-1.0d0)
         iter=iter+1
         if (ytry <= y(ilo)) then
            ytry=amotry(2.0d0)
            iter=iter+1
         else if (ytry >= y(inhi)) then
            ysave=y(ihi)
            ytry=amotry(0.5d0)
            iter=iter+1
            if (ytry >= ysave) then
               p(:,:)=0.5d0*(p(:,:)+spread(p(ilo,:),1,size(p,1)))
               do i=1,ndim+1
                  if (i /= ilo) y(i)=func(p(i,:))
               end do
               iter=iter+ndim
               psum(:)=sum(p(:,:),dim=1)
            end if
         end if
      end do
    END SUBROUTINE amoeba_private

    FUNCTION amotry(fac)
      IMPLICIT NONE
      REAL(8), INTENT(IN) :: fac
      REAL(8) :: amotry
      REAL(8) :: fac1,fac2,ytry
      REAL(8), DIMENSION(size(p,2)) :: ptry
      fac1=(1.0d0-fac)/ndim
      fac2=fac1-fac
      ptry(:)=psum(:)*fac1-p(ihi,:)*fac2
      ytry=func(ptry)
      if (ytry < y(ihi)) then
         y(ihi)=ytry
         psum(:)=psum(:)-p(ihi,:)+ptry(:)
         p(ihi,:)=ptry(:)
      end if
      amotry=ytry
    END FUNCTION amotry
  END SUBROUTINE Nelder_Meade
END MODULE SIMPLEX
