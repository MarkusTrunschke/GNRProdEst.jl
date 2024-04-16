MODULE MINIMIZATION
  USE NRUTIL, ONLY : vabs,assert_eq,nrerror,NR_EPS,MMUL,DOTP
  IMPLICIT NONE
  PRIVATE 
  PUBLIC BFGS!FRPR,BFGS_STATA
  
CONTAINS
    
  SUBROUTINE BFGS(theta,gtol,told,func,id,forward)
    IMPLICIT NONE
    REAL(8), INTENT(IN)    :: gtol,told
    INTEGER, INTENT(IN)    :: id        ! For MPI
    REAL(8), INTENT(INOUT) :: theta(:)	! parameter vector over which you which to maximize
    LOGICAL, INTENT(IN) :: forward
    INTERFACE	                        ! This lets the subroutine know that func is not a variable, it is a function
       REAL(8) FUNCTION func(theta)   ! This is not the likelihood itself, it just let the program know
         IMPLICIT NONE		  
         REAL(8), INTENT(IN) :: theta(:)
       END FUNCTION func
    END INTERFACE
    INTEGER, PARAMETER :: ITMAX=400	! Maximum number of iterations (in case it does not converge)
    REAL(8), PARAMETER :: STPMX=0.05d0,EPS=2.220446049250313D-9 !,TOLD=1.0d-12 ! Tolerance criterions
    INTEGER :: j,its ! indices
    LOGICAL :: check 
    REAL(8) :: den,fac,fad,fae,ftheta,sumdg,sumd ! working matrices
    REAL(8) :: dg(SIZE(theta)),g(SIZE(theta)),hdg(SIZE(theta)),thetanew(SIZE(theta)),d(SIZE(theta))
    REAL(8) :: B(SIZE(theta),SIZE(theta)),fret,stpmax,tim
    INTEGER :: tin,tout,hz,n
    CHARACTER(len=500) :: file_open
    n=SIZE(theta)
    ftheta=func(theta)        ! initial value of the function
!    IF (id==0) WRITE(*,'(A,F16.8)') 'First value of the likelihood ',ftheta
    CALL FLUSH(6)
    IF (forward) THEN
       g=gradient1(theta,ftheta,func)	! initial value of the gradient, if analytical gradient known here
    ELSE
       g=gradient2(theta,func)	! initial value of the gradient, if analytical gradient known here
    END IF
    B=0.0d0				! initialize the matrix to a negative definite matrix
    DO j = 1, n
       B(j,j)=1.0d0
!       IF (id==0) WRITE(*,*) j,g(j)
!       CALL FLUSH(6)
    END DO
    d=-MMUL(B,g)	! initial line direction
    stpmax=STPMX*MAX(vabs(theta),DBLE(n))
    DO its=1,ITMAX	! Main loop over iterations
       CALL SYSTEM_CLOCK(count_rate=hz) 
       CALL SYSTEM_CLOCK(count=tin) 
       CALL lnsrch(theta,ftheta,g,d,thetanew,fret,stpmax,check,func)
       ftheta=fret	! update the value of the fuction
       d=thetanew-theta	! update the direction of movement
       theta=thetanew	! update the point
       IF (MAXVAL(ABS(d)/MAX(ABS(theta),1.0D0)) < TOLD) THEN ! Check convergence in Delta-d
!          IF (id==0) THEN
!             WRITE(*,*) 'TOLD'
!             OPEN(6424,file="point.out")
!             OPEN(6425,file="hess.out")
!             DO j = 1, n
!                WRITE(6424,'(F32.16)') theta(j)
!                WRITE(6425,'(2000F32.16)') B(j,j)
!             END DO
!             CLOSE(6424)	
!             CLOSE(6425)	
!          END IF
          RETURN     
       END IF
       dg=g				! save the old gradient
       IF (forward) THEN
          g=gradient1(theta,ftheta,func)	! and get the new gradient
       ELSE
          g=gradient2(theta,func)	! and get the new gradient
       END IF
       den=MAX(ftheta,1.0d0)	
       IF (MAXVAL(ABS(g)*MAX(ABS(theta),1.0d0)/den) < gtol) THEN ! Check for convergence on zero gradient
!          IF (id==0) THEN
!             WRITE(*,*) 'GTOL'
!             OPEN(6424,file="point.out")
!             OPEN(6425,file="hess.out")
!             DO j = 1, n
!                WRITE(6424,'(F32.16)') theta(j)
!                WRITE(6425,'(2000F32.16)') B(j,j)
!             END DO
!             CLOSE(6424)	
!             CLOSE(6425)	
!          END IF
          RETURN
       END IF
       dg=g-dg			! compute difference in gradients
       hdg=MMUL(B,dg)	        ! and difference times current matrix
       fac=DOTP(dg,d)	! Dot products for denominators
       fae=DOTP(dg,hdg)
       sumdg=DOTP(dg,dg)
       sumd=DOTP(d,d)
       IF (fac > SQRT(EPS*sumdg*sumd)) THEN	! Skip update if fac not positive enough
          fac=1.0d0/fac
          fad=1.0d0/fae
          dg=fac*d-fad*hdg			! Vector that makes BFGS different from DFP
          ! BFGS updating formula
          B=B+fac*P_Outer_Product_Min(d,d)
          B=B-fad*P_Outer_Product_Min(hdg,hdg)
          B=B+fae*P_Outer_Product_Min(dg,dg)
       END IF
       d=-MMUL(B,g)	! Next direction to go
!       IF (id==0) THEN
!          OPEN(6424,file="point.out")
!          DO j = 1, n
!             WRITE(6424,'(F32.16)') theta(j)
!          END DO
!          CLOSE(6424)	
!       END IF
       CALL SYSTEM_CLOCK(count=tout) 
       tout=tout-tin
       tim=DBLE(tout)/DBLE(hz)
!       IF (id==0) THEN
!          write(*,'(A,I10,A,F32.16,A,F16.6)') 'iteration:',its,' value:',fret,' time:',tim
!          CALL FLUSH(6)
!       END IF
    END DO		! Go back for another iteration
!    CALL NRERROR('BFGS: too many iterations')
!    IF (id==0) THEN
!       WRITE(*,*) 'NO GTOL'
!       OPEN(6424,file="point.out")
!       OPEN(6425,file="hess.out")
!       DO j = 1, n
!          WRITE(6424,'(F32.16)') theta(j)
!          WRITE(6425,'(2000F32.16)') B(j,j)
!       END DO
!       CLOSE(6424)	
!       CLOSE(6425)	
!    END IF
!    WRITE(*,*) 'BFGS:Too many iterations'
  END SUBROUTINE BFGS
  
!!$  SUBROUTINE BFGS_STATA(theta,gtol,told,func,id,forward)
!!$    IMPLICIT NONE
!!$    REAL(8), INTENT(IN)    :: gtol,told
!!$    INTEGER, INTENT(IN)    :: id        ! For MPI
!!$    REAL(8), INTENT(INOUT) :: theta(:)	! parameter vector over which you which to maximize
!!$    LOGICAL, INTENT(IN) :: forward
!!$    INTERFACE	                        ! This lets the subroutine know that func is not a variable, it is a function
!!$       REAL(8) FUNCTION func(theta)   ! This is not the likelihood itself, it just let the program know
!!$         IMPLICIT NONE		  
!!$         REAL(8), INTENT(IN) :: theta(:)
!!$       END FUNCTION func
!!$    END INTERFACE
!!$    INTEGER, PARAMETER :: ITMAX=100	! Maximum number of iterations (in case it does not converge)
!!$    REAL(8), PARAMETER :: STPMX=0.01d0,EPS=2.220446049250313D-9 !,TOLD=1.0d-12 ! Tolerance criterions
!!$    INTEGER :: j,its ! indices
!!$    LOGICAL :: check 
!!$    REAL(8) :: den,fac,fad,fae,ftheta,sumdg,sumd ! working matrices
!!$    REAL(8) :: dg(SIZE(theta)),g(SIZE(theta)),hdg(SIZE(theta)),thetanew(SIZE(theta)),d(SIZE(theta))
!!$    REAL(8) :: B(SIZE(theta),SIZE(theta)),fret,stpmax,tim,SD(SIZE(theta))
!!$    INTEGER :: tin,tout,hz,n
!!$    CHARACTER(len=500) :: file_open
!!$    n=SIZE(theta)
!!$    ftheta=func(theta)        ! initial value of the function
!!$    IF (id==0) WRITE(*,'(A,F16.8)') 'First value of the likelihood ',ftheta
!!$    CALL FLUSH(6)
!!$    IF (forward) THEN
!!$       g=gradient1(theta,ftheta,func)	! initial value of the gradient, if analytical gradient known here
!!$    ELSE
!!$       g=gradient2(theta,func)	! initial value of the gradient, if analytical gradient known here
!!$    END IF
!!$    B=0.0d0				! initialize the matrix to a negative definite matrix
!!$    DO j = 1, n
!!$       B(j,j)=1.0d0
!!$       IF (id==0) WRITE(*,*) j,g(j)
!!$       CALL FLUSH(6)
!!$    END DO
!!$    d=-MMUL(B,g)	! initial line direction
!!$    stpmax=STPMX*MAX(vabs(theta),DBLE(n))
!!$    DO its=1,ITMAX	! Main loop over iterations
!!$       CALL SYSTEM_CLOCK(count_rate=hz) 
!!$       CALL SYSTEM_CLOCK(count=tin) 
!!$       CALL lnsrch(theta,ftheta,g,d,thetanew,fret,stpmax,check,func)
!!$       ftheta=fret	! update the value of the fuction
!!$       d=thetanew-theta	! update the direction of movement
!!$       theta=thetanew	! update the point
!!$       IF (MAXVAL(ABS(d)/MAX(ABS(theta),1.0D0)) < TOLD) THEN ! Check convergence in Delta-d
!!$          IF (id==0) THEN
!!$             WRITE(*,*) 'TOLD'
!!$             OPEN(6424,file="point.out")
!!$             DO j = 1, n
!!$                WRITE(6424,'(F32.16)') theta(j)
!!$             END DO
!!$             CLOSE(6424)	
!!$          END IF
!!$          RETURN     
!!$       END IF
!!$       dg=g				! save the old gradient
!!$       IF (forward) THEN
!!$          g=gradient1(theta,ftheta,func)	! and get the new gradient
!!$       ELSE
!!$          g=gradient2(theta,func)	! and get the new gradient
!!$       END IF
!!$       den=MAX(ftheta,1.0d0)	
!!$       IF (MAXVAL(ABS(g)*MAX(ABS(theta),1.0d0)/den) < gtol) THEN ! Check for convergence on zero gradient
!!$          IF (id==0) THEN
!!$             WRITE(*,*) 'GTOL'
!!$             OPEN(6424,file="point.out")
!!$             DO j = 1, n
!!$                WRITE(6424,'(F32.16)') theta(j)
!!$             END DO
!!$             CLOSE(6424)	
!!$          END IF
!!$          RETURN
!!$       END IF
!!$       dg=g-dg			! compute difference in gradients
!!$       hdg=MMUL(B,dg)	        ! and difference times current matrix
!!$       fac=DOTP(dg,d)	! Dot products for denominators
!!$       fae=DOTP(dg,hdg)
!!$       sumdg=DOTP(dg,dg)
!!$       sumd=DOTP(d,d)
!!$       IF (ABS(fac*fae) > EPS*sumdg*sumd) THEN	! Skip update if fac not positive enough
!!$          fac=1.0d0/fac
!!$          fad=1.0d0/fae
!!$          SD=d*fac - hdg*fad
!!$          B = B - fac*P_Outer_Product_Min(d,d)
!!$          B = B + fae*P_Outer_Product_Min(SD,SD)
!!$          B = B - fad*P_Outer_Product_Min(hdg,hdg)
!!$          B = 0.5d0*(B + TRANSPOSE(B))
!!$       END IF
!!$       d=-MMUL(B,g)	! Next direction to go
!!$       IF (id==0) THEN
!!$          OPEN(6424,file="point.out")
!!$          DO j = 1, n
!!$             WRITE(6424,'(F32.16)') theta(j)
!!$          END DO
!!$          CLOSE(6424)	
!!$       END IF
!!$       CALL SYSTEM_CLOCK(count=tout) 
!!$       tout=tout-tin
!!$       tim=DBLE(tout)/DBLE(hz)
!!$       IF (id==0) THEN
!!$          write(*,'(A,I10,A,F32.16,A,F16.6)') 'iteration:',its,' value:',fret,' time:',tim
!!$          CALL FLUSH(6)
!!$       END IF
!!$    END DO		! Go back for another iteration
!!$!    CALL NRERROR('BFGS: too many iterations')
!!$    WRITE(*,*) 'BFGS_STATA:Too many iterations'
!!$  END SUBROUTINE BFGS_STATA
!!$
!!$  SUBROUTINE FRPR(p,ftol,iter,fret,func,id,forward)
!!$    USE nrutil, ONLY : nrerror
!!$    IMPLICIT NONE
!!$    INTERFACE
!!$       FUNCTION func(p)
!!$         IMPLICIT NONE
!!$         REAL(8), DIMENSION(:), INTENT(IN) :: p
!!$         REAL(8) :: func
!!$       END FUNCTION func
!!$    END INTERFACE
!!$    LOGICAL, INTENT(IN) :: forward
!!$    INTEGER, INTENT(IN) :: id
!!$    INTEGER, INTENT(OUT) :: iter
!!$    REAL(8), INTENT(IN) :: ftol
!!$    REAL(8), INTENT(OUT) :: fret
!!$    REAL(8), DIMENSION(:), INTENT(INOUT) :: p
!!$    INTEGER, PARAMETER :: ITMAX=200
!!$    REAL(8), PARAMETER :: EPS=1.0d-10
!!$    INTEGER :: its,tin,tout,hz,j
!!$    REAL(8) :: dgg,fp,gam,gg,tim
!!$    REAL(8), DIMENSION(size(p)) :: g,h,xi
!!$    fp=func(p)
!!$    IF (id==0) WRITE(*,'(A,F16.8)') 'First value of the likelihood ',fp
!!$    CALL FLUSH(6)
!!$    IF (forward) THEN
!!$       xi=gradient1(p,fp,func)	! initial value of the gradient, if analytical gradient known here
!!$    ELSE
!!$       xi=gradient2(p,func)
!!$    END IF
!!$    DO j = 1,SIZE(p)
!!$       IF (id==0) WRITE(*,*) j,xi(j)
!!$       CALL FLUSH(6)
!!$    END DO
!!$    g=-xi
!!$    h=g
!!$    xi=h
!!$    DO its=1,ITMAX
!!$       CALL SYSTEM_CLOCK(count_rate=hz) 
!!$       CALL SYSTEM_CLOCK(count=tin) 
!!$       iter=its
!!$       CALL LINMIN(p,xi,fret,func)
!!$       IF (id==0) THEN
!!$          WRITE(*,*) 'Linmin'
!!$          CALL FLUSH(6)
!!$          OPEN(6424,file="point.out")
!!$          DO j = 1, SIZE(p)
!!$             WRITE(6424,'(F32.16)') p(j)
!!$          END DO
!!$          CLOSE(6424)	
!!$       END IF
!!$       IF (2.0d0*abs(fret-fp) <= ftol*(ABS(fret)+ABS(fp)+EPS)) RETURN
!!$       fp=fret
!!$       IF (forward) THEN
!!$          xi=gradient1(p,fp,func)	! initial value of the gradient, if analytical gradient known here
!!$       ELSE
!!$          xi=gradient2(p,func)
!!$       END IF
!!$       gg=DOTP(g,g)
!!$       !    dgg=dot_product(xi,xi)
!!$       dgg=DOTP(xi+g,xi)
!!$       IF (gg == 0.0d0) RETURN
!!$       gam=dgg/gg
!!$       g=-xi
!!$       h=g+gam*h
!!$       xi=h
!!$       CALL SYSTEM_CLOCK(count=tout) 
!!$       tout=tout-tin
!!$       tim=DBLE(tout)/DBLE(hz)
!!$       IF (id==0) THEN
!!$          write(*,'(A,I10,A,F32.16,A,F16.6)') 'iteration:',its,' value:',fret,' time:',tim
!!$          CALL FLUSH(6)
!!$       END IF
!!$    END DO
!!$    CALL nrerror('frpr: maximum iterations exceeded')
!!$  END SUBROUTINE FRPR
!!$  
!!$  SUBROUTINE linmin(p,xi,fret,func)
!!$    USE nrutil, ONLY : assert_eq
!!$    USE FRPR_PRIVATE_F1DIM
!!$    IMPLICIT NONE
!!$    INTERFACE
!!$       FUNCTION func(p)
!!$         IMPLICIT NONE
!!$         REAL(8), DIMENSION(:), INTENT(IN) :: p
!!$         REAL(8) :: func
!!$       END FUNCTION func
!!$    END INTERFACE
!!$    REAL(8), INTENT(OUT) :: fret
!!$    REAL(8), DIMENSION(:), TARGET, INTENT(INOUT) :: p,xi
!!$    REAL(8), PARAMETER :: TOL=1.0d-4
!!$    REAL(8) :: ax,bx,fa,fb,fx,xmin,xx
!!$    ncom=assert_eq(size(p),size(xi),'linmin')
!!$    pcom=>p
!!$    xicom=>xi
!!$    ax=0.0d0
!!$    xx=1.0d0
!!$    CALL mnbrak(ax,xx,bx,fa,fx,fb,f1dim,func)
!!$    WRITE(*,*) 'mnbrak'
!!$    fret=brent(ax,xx,bx,f1dim,TOL,xmin,func)
!!$    xi=xmin*xi
!!$    p=p+xi
!!$  END SUBROUTINE linmin
!!$  
!!$  FUNCTION brent(ax,bx,cx,feval,tol,xmin,func)
!!$    USE nrutil, ONLY : nrerror
!!$    IMPLICIT NONE
!!$    REAL(8), INTENT(IN) :: ax,bx,cx,tol
!!$    REAL(8), INTENT(OUT) :: xmin
!!$    REAL(8) :: brent
!!$    INTERFACE
!!$       FUNCTION func(p)
!!$         IMPLICIT NONE
!!$         REAL(8), DIMENSION(:), INTENT(IN) :: p
!!$         REAL(8) :: func
!!$       END FUNCTION func
!!$       FUNCTION feval(x,func)
!!$         IMPLICIT NONE
!!$         INTERFACE
!!$            FUNCTION func(p)
!!$              IMPLICIT NONE
!!$              REAL(8), DIMENSION(:), INTENT(IN) :: p
!!$              REAL(8) :: func
!!$            END FUNCTION func
!!$         END INTERFACE
!!$         REAL(8), INTENT(IN) :: x
!!$         REAL(8) :: feval
!!$       END FUNCTION feval
!!$    END INTERFACE
!!$    INTEGER, PARAMETER :: ITMAX=100
!!$    REAL(8), PARAMETER :: CGOLD=0.3819660d0,ZEPS=1.0d-3*epsilon(ax)
!!$    INTEGER :: iter
!!$    REAL(8) :: a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
!!$    a=MIN(ax,cx)
!!$    b=MAX(ax,cx)
!!$    v=bx
!!$    w=v
!!$    x=v
!!$    e=0.0d0
!!$    fx=feval(x,func)
!!$    fv=fx
!!$    fw=fx
!!$    DO iter=1,ITMAX
!!$       xm=0.5d0*(a+b)
!!$       tol1=tol*ABS(x)+ZEPS
!!$       tol2=2.0d0*tol1
!!$       IF (ABS(x-xm) <= (tol2-0.5d0*(b-a))) THEN
!!$          xmin=x
!!$          brent=fx
!!$          RETURN
!!$       END IF
!!$       IF (ABS(e) > tol1) THEN
!!$          r=(x-w)*(fx-fv)
!!$          q=(x-v)*(fx-fw)
!!$          p=(x-v)*q-(x-w)*r
!!$          q=2.0d0*(q-r)
!!$          IF (q > 0.0d0) p=-p
!!$          q=ABS(q)
!!$          etemp=e
!!$          e=d
!!$          IF (ABS(p) >= ABS(0.5d0*q*etemp) .OR. &
!!$               p <= q*(a-x) .OR. p >= q*(b-x)) THEN
!!$             e=MERGE(a-x,b-x, x >= xm )
!!$             d=CGOLD*e
!!$          ELSE
!!$             d=p/q
!!$             u=x+d
!!$             IF (u-a < tol2 .OR. b-u < tol2) d=SIGN(tol1,xm-x)
!!$          END IF
!!$       ELSE
!!$          e=MERGE(a-x,b-x, x >= xm )
!!$          d=CGOLD*e
!!$       END IF
!!$       u=MERGE(x+d,x+SIGN(tol1,d), ABS(d) >= tol1 )
!!$       fu=feval(u,func)
!!$       IF (fu <= fx) THEN
!!$          IF (u >= x) THEN
!!$             a=x
!!$          ELSE
!!$             b=x
!!$          END IF
!!$          CALL shft(v,w,x,u)
!!$          CALL shft(fv,fw,fx,fu)
!!$       ELSE
!!$          IF (u < x) THEN
!!$             a=u
!!$          ELSE
!!$             b=u
!!$          END IF
!!$          IF (fu <= fw .OR. w == x) THEN
!!$             v=w
!!$             fv=fw
!!$             w=u
!!$             fw=fu
!!$          ELSE IF (fu <= fv .OR. v == x .OR. v == w) THEN
!!$             v=u
!!$             fv=fu
!!$          END IF
!!$       END IF
!!$    END DO
!!$    CALL nrerror('brent: exceed maximum iterations')
!!$  CONTAINS
!!$    
!!$    SUBROUTINE shft(a,b,c,d)
!!$      REAL(8), INTENT(OUT) :: a
!!$      REAL(8), INTENT(INOUT) :: b,c
!!$      REAL(8), INTENT(IN) :: d
!!$      a=b
!!$      b=c
!!$      c=d
!!$    END SUBROUTINE shft
!!$  END FUNCTION brent
!!$  
!!$  SUBROUTINE mnbrak(ax,bx,cx,fa,fb,fc,feval,func)
!!$    USE nrutil, ONLY : swap
!!$    IMPLICIT NONE
!!$    REAL(8), INTENT(INOUT) :: ax,bx
!!$    REAL(8), INTENT(OUT) :: cx,fa,fb,fc
!!$    INTERFACE
!!$       FUNCTION func(p)
!!$         IMPLICIT NONE
!!$         REAL(8), DIMENSION(:), INTENT(IN) :: p
!!$         REAL(8) :: func
!!$       END FUNCTION func
!!$       FUNCTION feval(x,func)
!!$         IMPLICIT NONE
!!$         INTERFACE
!!$            FUNCTION func(p)
!!$              IMPLICIT NONE
!!$              REAL(8), DIMENSION(:), INTENT(IN) :: p
!!$              REAL(8) :: func
!!$            END FUNCTION func
!!$         END INTERFACE
!!$         REAL(8), INTENT(IN) :: x
!!$         REAL(8) :: feval
!!$       END FUNCTION feval
!!$    END INTERFACE
!!$    REAL(8), PARAMETER :: GOLD=1.618034d0,GLIMIT=100.0d0,TINY=1.0d-20
!!$    REAL(8) :: fu,q,r,u,ulim
!!$    fa=feval(ax,func)
!!$    fb=feval(bx,func)
!!$    IF (fb > fa) THEN
!!$       CALL swap(ax,bx)
!!$       CALL swap(fa,fb)
!!$    END IF
!!$    cx=bx+GOLD*(bx-ax)
!!$    fc=feval(cx,func)
!!$    DO
!!$       IF (fb < fc) RETURN
!!$       r=(bx-ax)*(fb-fc)
!!$       q=(bx-cx)*(fb-fa)
!!$       u=bx-((bx-cx)*q-(bx-ax)*r)/(2.0d0*sign(max(abs(q-r),TINY),q-r))
!!$       ulim=bx+GLIMIT*(cx-bx)
!!$       if ((bx-u)*(u-cx) > 0.0d0) then
!!$          fu=feval(u,func)
!!$          if (fu < fc) then
!!$             ax=bx
!!$             fa=fb
!!$             bx=u
!!$             fb=fu
!!$             RETURN
!!$          else if (fu > fb) then
!!$             cx=u
!!$             fc=fu
!!$             RETURN
!!$          end if
!!$          u=cx+GOLD*(cx-bx)
!!$          fu=feval(u,func)
!!$       else if ((cx-u)*(u-ulim) > 0.0d0) then
!!$          fu=feval(u,func)
!!$          if (fu < fc) then
!!$             bx=cx
!!$             cx=u
!!$             u=cx+GOLD*(cx-bx)
!!$             call shft(fb,fc,fu,feval(u,func))
!!$          end if
!!$       else if ((u-ulim)*(ulim-cx) >= 0.0d0) then
!!$          u=ulim
!!$          fu=feval(u,func)
!!$       else
!!$          u=cx+GOLD*(cx-bx)
!!$          fu=feval(u,func)
!!$       end if
!!$       call shft(ax,bx,cx,u)
!!$       call shft(fa,fb,fc,fu)
!!$    end do
!!$  CONTAINS
!!$    
!!$    SUBROUTINE shft(a,b,c,d)
!!$      REAL(8), INTENT(OUT) :: a
!!$      REAL(8), INTENT(INOUT) :: b,c
!!$      REAL(8), INTENT(IN) :: d
!!$      a=b
!!$      b=c
!!$      c=d
!!$    END SUBROUTINE shft
!!$  END SUBROUTINE mnbrak

  FUNCTION GRADIENT1(x0,f0,func)
    ! Routine to get the numerical gradient of a function using forward differences
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: x0(:),f0
    INTERFACE	                       ! This lets the subroutine know that func is not a variable, it is a function
       REAL(8) FUNCTION func(theta)  ! This is not the likelihood itself, it just let the program know
         IMPLICIT NONE		       ! the arguments and type of arguments it accepts
         REAL(8), INTENT(IN) :: theta(:)
       END FUNCTION func
    END INTERFACE
    REAL(8), PARAMETER :: EPS=1.49011611938477D-08 !SQRT(NR_EPS)
    REAL(8) :: GRADIENT1(SIZE(x0)),grdd(SIZE(x0)),dh(SIZE(x0)),xdh(SIZE(x0))
    REAL(8) :: arg(SIZE(x0),SIZE(x0))
    INTEGER :: i,k
    k=SIZE(x0)
    grdd = 0.0D0
    ! Computation of stepsize (dh) for gradient
    dh = EPS*ABS(x0)
    WHERE (dh<1.0d-13) dh=EPS
    xdh = x0+dh
    dh=xdh-x0
    arg = SPREAD(x0,DIM=2,NCOPIES=k)
    DO i = 1, k
       arg(i,i)=xdh(i)    
       grdd(i)=func(arg(:,i))
    END DO
    GRADIENT1 = (grdd-f0)/dh
  END FUNCTION GRADIENT1
  
  FUNCTION GRADIENT2(x0,func)
    ! Routine to get the numerical gradient of a function using central differences
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: x0(:)
    INTERFACE	                       ! This lets the subroutine know that func is not a variable, it is a function
       REAL(8) FUNCTION func(theta)  ! This is not the likelihood itself, it just let the program know
         IMPLICIT NONE		       ! the arguments and type of arguments it accepts
         REAL(8), INTENT(IN) :: theta(:)
       END FUNCTION func
    END INTERFACE
    REAL(8), PARAMETER :: EPS=6.05545445239334D-06 !NR_EPS**(1.0d0/3.0d0)
    REAL(8) :: GRADIENT2(SIZE(x0)),grdd(SIZE(x0)),dh(SIZE(x0)),xdh1(SIZE(x0)),xdh0(SIZE(x0)),dhh(SIZE(x0))
    REAL(8) :: arg(SIZE(x0),SIZE(x0))
    INTEGER :: i,k
    k=SIZE(x0)
    grdd = 0.0D0
    ! Computation of stepsize (dh) for gradient
    dh = EPS*ABS(x0)
    WHERE (dh==0.0D0) dh=EPS
    xdh1 = x0+dh
    arg = SPREAD(x0,DIM=2,NCOPIES=k)
    DO i = 1, k
       arg(i,i)=xdh1(i)    
       grdd(i)=func(arg(:,i))
    END DO
    xdh0 = x0-dh
    dhh = xdh1-xdh0
    arg = SPREAD(x0,DIM=2,NCOPIES=k)
    DO i = 1, k
       arg(i,i)=xdh0(i)    
       grdd(i)=grdd(i)-func(arg(:,i))
    END DO
    GRADIENT2 = grdd/dhh
  END FUNCTION GRADIENT2

  FUNCTION GRADIENT3(x0,func)
    ! Routine to get the numerical gradient of a function using a 5 point stencil
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: x0(:)
    INTERFACE	                       ! This lets the subroutine know that func is not a variable, it is a function
       REAL(8) FUNCTION func(theta)  ! This is not the likelihood itself, it just let the program know
         IMPLICIT NONE		       ! the arguments and type of arguments it accepts
         REAL(8), INTENT(IN) :: theta(:)
       END FUNCTION func
    END INTERFACE
    REAL(8), PARAMETER :: EPS=6.05545445239334D-06!NR_EPS**(1.0d0/3.0d0)
    REAL(8) :: GRADIENT3(SIZE(x0)),grdd(SIZE(x0)),dh(SIZE(x0)),xdh(SIZE(x0))
    REAL(8) :: arg(SIZE(x0),SIZE(x0))
    INTEGER :: i,k
    k=SIZE(x0)
    grdd = 0.0D0
    ! Computation of stepsize (dh) for gradient
    dh = EPS*ABS(x0)
    WHERE (dh==0.0D0) dh=EPS
    xdh = x0+dh
    dh = xdh-x0
    arg = SPREAD(x0,DIM=2,NCOPIES=k)
    DO i = 1, k
       arg(i,i)=xdh(i)    
       grdd(i)=func(arg(:,i))
    END DO
    xdh = x0-dh
    dh = x0-xdh
    arg = SPREAD(x0,DIM=2,NCOPIES=k)
    DO i = 1, k
       arg(i,i)=xdh(i)    
       grdd(i)=grdd(i)-func(arg(:,i))
    END DO
    GRADIENT3 = grdd/(2.0d0*dh)
  END FUNCTION GRADIENT3

  SUBROUTINE lnsrch(xold,fold,g,p,x,f,stpmax,check,func)
    ! Given an N-dimensional point xold, the value of the function and gradient there, fold
    ! and g, and a direction p, finds a new point x along the direction p from xold where the
    ! function func has decreased sufficiently. xold, g, p, and x are all arrays of length N.
    ! The new function value is returned in f. stpmax is an input quantity that limits the length
    ! of the steps so that you do not try to evaluate the function in regions where it is undefined
    ! or subject to overflow. p is usually the Newton direction. The output quantity check is
    ! false on a normal exit. It is true when x is too close to xold. In a minimization algorithm,
    ! this usually signals convergence and can be ignored. However, in a zero-finding algorithm
    ! the calling program should check whether the convergence is spurious.
    ! Parameters: ALF ensures su.cient decrease in function value; TOLX is the convergence
    ! criterion on .x.
    IMPLICIT NONE
    REAL(8), DIMENSION(:), INTENT(IN) :: xold,g
    REAL(8), DIMENSION(:), INTENT(INOUT) :: p(SIZE(xold))
    REAL(8), INTENT(IN) :: fold,stpmax
    REAL(8), INTENT(OUT) :: x(SIZE(xold))
    REAL(8), INTENT(OUT) :: f
    LOGICAL, INTENT(OUT) :: check
    INTERFACE	! This lets the subroutine know that func is not a variable, it is a function
       REAL(8) FUNCTION func(theta)  ! This is not the likelihood itself, it just let the program know
         IMPLICIT NONE		
         REAL(8), INTENT(IN) :: theta(:)
       END FUNCTION func
    END INTERFACE
    REAL(8), PARAMETER :: ALF=1.0d-4,TOLX=1.0d-9
    INTEGER :: ndum
    REAL(8) :: a,alam,alam2,alamin,b,disc,f2,pabs,rhs1,rhs2,slope,tmplam,zero
    zero=0.0D0
    ndum=assert_eq(SIZE(g),SIZE(xold),'Dimensions do not agree: LNSRCH')
    check=.false.
    pabs=vabs(p(:))
    IF (pabs > stpmax) p(:)=p(:)*stpmax/pabs ! Scale if attempted step is too big.
    slope=DOTP(g,p)
    IF (slope >= zero) CALL NRERROR('roundoff problem in lnsrch')
    alamin=TOLX/MAXVAL(DABS(p(:))/MAX(DABS(xold(:)),1.0D0)) !Compute .min.
    alam=1.0D0 ! Always try full Newton step first.
    DO ! Start of iteration loop.
       x(:)=xold(:)+alam*p(:)
       F=FUNC(x)
       IF (alam < alamin) THEN                 ! Convergence on .x. For zero finding,
          !the calling program should verify the convergence.
          x(:)=xold(:)
          check=.true.
          RETURN
       ELSE IF (f <= fold+ALF*alam*slope) THEN ! Sufficient function decrease.
          RETURN
       ELSE !Backtrack.
          IF (alam == 1.0D0) THEN !First time.
             tmplam=-slope/(2.0D0*(f-fold-slope))
          ELSE !Subsequent backtracks.
             rhs1=f-fold-alam*slope
             rhs2=f2-fold-alam2*slope
             a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
             b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/&
                  (alam-alam2)
             IF (a == 0.0D0) THEN
                tmplam=-slope/(2.0D0*b)
             ELSE
                disc=b*b-3.0*a*slope
                IF (disc < 0.0D0) THEN
                   tmplam=0.5D0*alam
                ELSE IF (b <= 0.0D0) THEN
                   tmplam=(-b+SQRT(disc))/(3.0D0*a)
                ELSE
                   tmplam=-slope/(b+SQRT(disc))
                END IF
             END IF
             IF (tmplam > 0.5D0*alam) tmplam=0.5D0*alam
          END IF
       END IF
       alam2=alam
       f2=f
       alam=MAX(tmplam,0.1D0*alam) 
    END DO !Try again.
  END SUBROUTINE lnsrch
  
  FUNCTION P_Outer_Product_Min(a,b)
    IMPLICIT NONE
    REAL(8), DIMENSION(:), INTENT(IN) :: a,b
    REAL(8), DIMENSION(size(a),size(b)) :: P_Outer_Product_Min
    P_Outer_Product_Min = spread(a,dim=2,ncopies=size(b)) * &
         spread(b,dim=1,ncopies=size(a))
  END FUNCTION P_Outer_Product_Min
  
END MODULE MINIMIZATION
