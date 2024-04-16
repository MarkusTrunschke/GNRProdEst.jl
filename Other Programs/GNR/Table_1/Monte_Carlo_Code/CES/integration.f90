MODULE INTEGRATION
  USE NRUTIL, ONLY : nrerror,arth,iminloc,NR_PI,NR_EPS,assert_eq
  IMPLICIT NONE
  
  PRIVATE :: P_POLINT_INT,P_TRAPZD_INT
  
CONTAINS
  
  REAL(8) FUNCTION Simpson(func,a,b,eps,i)
    ! Returns the integral of the function func from a to b. The parameter EPS should be set to
    ! the desired fractional accuracy and JMAX so that 2 to the power JMAX-1 is the maximum
    ! allowed number of steps. Integration is performed by Simpson's rule.
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i
    REAL(8), INTENT(IN) :: a,b,eps
    INTERFACE
       REAL(8) FUNCTION func(x,i)
         INTEGER, INTENT(IN) :: i
         REAL(8), INTENT(IN) :: x
       END FUNCTION func
    END INTERFACE
    INTEGER, PARAMETER :: JMAX=50
    INTEGER :: j
    REAL(8) :: os,ost,st
    ost=0.0d0
    os= 0.0d0
    DO j=1,JMAX
       CALL P_TRAPZD_INT(func,a,b,st,j,i)
       Simpson=(4.0d0*st-ost)/3.0d0
       IF (j > 5) THEN
          IF (ABS(Simpson-os) < eps*ABS(os) .OR. (Simpson == 0.0d0 .AND. os == 0.0d0)) RETURN
       END IF
       os=Simpson
       ost=st
    END DO
    CALL NRERROR ('too many steps: Simpson')
  END FUNCTION Simpson
  
  REAL(8) FUNCTION Romberg(func,a,b,eps,i)
    ! Returns the integral of the function func from a to b. Integration is performed by Romberg's
    ! method of order 2K, where, e.g., K=2 is Simpson's rule.
    ! Parameters: EPS is the fractional accuracy desired, as determined by the extrapolation error
    ! estimate; JMAX limits the total number of steps; K is the number of points used in the
    ! extrapolation.
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: a,b,eps
    INTERFACE
       REAL(8) FUNCTION func(x,i)
         INTEGER, INTENT(IN) :: i
         REAL(8), INTENT(IN) :: x
       END FUNCTION func
    END INTERFACE
    INTEGER, PARAMETER :: JMAX=40,JMAXP=JMAX+1,K=15,KM=K-1
    REAL(8), DIMENSION(JMAXP) :: h,s
    INTEGER, INTENT(IN) :: i
    REAL(8) :: dqromb
    INTEGER :: j
    h(1)=1.0d0
    DO j=1,JMAX
       CALL P_TRAPZD_INT(func,a,b,s(j),j,i)
       IF (j >= K) THEN
          CALL P_POLINT_INT(h(j-KM:j),s(j-KM:j),0.0d0,Romberg,dqromb)
          IF (ABS(dqromb) <= eps*ABS(Romberg)) RETURN
       END IF
       s(j+1)=s(j)
       h(j+1)=0.25d0*h(j)
    END DO
    CALL NRERROR('too many steps: Romberg')
  END FUNCTION Romberg

  REAL(8) FUNCTION Gauss_16(func,a,b,i)
    INTEGER, INTENT(IN) :: i
    REAL(8), INTENT(IN) :: a,b
    INTERFACE
       REAL(8) FUNCTION func(x,i)
         INTEGER, INTENT(IN) :: i
         REAL(8), INTENT(IN) :: x
       END FUNCTION func
    END INTERFACE
    REAL(8) :: xm,xr,ss
    REAL(8), DIMENSION(8) :: dx,w = (/0.18945061045507d0,0.18260341504492d0,0.16915651939500d0, &
         0.14959598881658d0,0.12462897125554d0,0.09515851168249d0,0.06225352393865d0,0.02715245941175d0/),&
         x = (/0.09501250983764d0,0.28160355077926d0,0.45801677765723d0,0.61787624440264d0, &
         0.75540440835500d0,0.86563120238783d0,0.94457502307323d0,0.98940093499165d0/)
    INTEGER :: j
    xm=0.5d0*(b+a)
    xr=0.5d0*(b-a)
    ss=0.0d0
    dx=xr*x
    DO j = 1, 8
       ss=ss+w(j)*(func(xm+dx(j),i)+func(xm-dx(j),i))
    END DO
    ss=xr*ss
    Gauss_16=ss
  END FUNCTION Gauss_16

  REAL(8) FUNCTION Gauss_8(func,a,b,i)
    INTEGER, INTENT(IN) :: i
    REAL(8), INTENT(IN) :: a,b
    INTERFACE
       REAL(8) FUNCTION func(x,i)
         INTEGER, INTENT(IN) :: i
         REAL(8), INTENT(IN) :: x
       END FUNCTION func
    END INTERFACE
    REAL(8) :: xm,xr,ss
    REAL(8), DIMENSION(4) :: dx,w = (/0.36268378337836d0,0.31370664587789d0,0.22238103445337d0, &
         0.10122853629038d0/),&
         x = (/0.18343464249565d0,0.52553240991633d0,0.79666647741363d0,0.96028985649754d0/)
    INTEGER :: j
    xm=0.5d0*(b+a)
    xr=0.5d0*(b-a)
    ss=0.0d0
    dx=xr*x
    DO j = 1, 4
       ss=ss+w(j)*(func(xm+dx(j),i)+func(xm-dx(j),i))
    END DO
    ss=xr*ss
    Gauss_8=ss
  END FUNCTION Gauss_8

  RECURSIVE FUNCTION Adaptive_Quadrature_16_8(func,lo,hi,eps,i) RESULT (res)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: lo,hi,eps
    INTEGER, INTENT(IN) :: i
    INTERFACE
       REAL(8) FUNCTION func(x,i)
         INTEGER, INTENT(IN) :: i
         REAL(8), INTENT(IN) :: x
       END FUNCTION func
    END INTERFACE
    REAL(8) :: res
    REAL(8) :: g8,g16,middle_point
    REAL(8) :: left_area,right_area
    REAL(8) :: diff
    middle_point=0.5d0*(hi+lo)
    g8 = Gauss_8(func,lo,hi,i)
    g16 = Gauss_16(func,lo,hi,i)
    diff = ABS(g16-g8)
    IF ((diff < eps*ABS(g8)).OR.((g16==0.0).AND.g8==0.0).OR.(ABS(lo-middle_point)<=0.01d0)) THEN
       res = g16
    ELSE
       left_area = Adaptive_Quadrature_16_8(func,lo,middle_point,eps,i)
       right_area = Adaptive_Quadrature_16_8(func,middle_point,hi,eps,i)
       res = left_area + right_area
    END IF
  END FUNCTION Adaptive_Quadrature_16_8
  
  SUBROUTINE Gauss_Legendre(x,w,x1,x2)
    !Given the lower and upper limits of integration x1 and x2, this routine returns arrays x and w
    !of length N containing the abscissas and weights of the Gauss-Legendre N-point quadrature
    !formula. The parameter EPS is the relative precision.
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: x1,x2
    REAL(8), DIMENSION(:), INTENT(OUT) :: x,w
    INTEGER :: its,j,m,n
    INTEGER, PARAMETER :: MAXIT=20
    REAL(8) :: xl,xm
    REAL(8), DIMENSION((size(x)+1)/2) :: p1,p2,p3,pp,z,z1
    LOGICAL, DIMENSION((size(x)+1)/2) :: unfinished
    n=SIZE(x)
    m=(n+1)/2
    xm=0.5d0*(x2+x1)
    xl=0.5d0*(x2-x1)
    z=COS(NR_PI*(arth(1,1,m)-0.25d0)/(n+0.5d0))
    unfinished=.TRUE.
    DO its=1,MAXIT
       WHERE (unfinished)
          p1=1.0d0
          p2=0.0d0
       END WHERE
       DO j=1,n
          WHERE (unfinished)
             p3=p2
             p2=p1
             p1=((2.0d0*j-1.0d0)*z*p2-(j-1.0d0)*p3)/j
          END WHERE
       END DO
       WHERE (unfinished)
          pp=n*(z*p1-p2)/(z*z-1.0d0)
          z1=z
          z=z1-p1/pp
          unfinished=(ABS(z-z1) > NR_EPS)
       END WHERE
       IF (.NOT. ANY(unfinished)) EXIT
    END DO
    IF (its == MAXIT+1) CALL NRERROR('too many iterations: Gauss Legendre')
    x(1:m)=xm-xl*z
    x(n:n-m+1:-1)=xm+xl*z
    w(1:m)=2.0d0*xl/((1.0d0-z**2)*pp**2)
    w(n:n-m+1:-1)=w(1:m)
  END SUBROUTINE Gauss_Legendre

  SUBROUTINE Gauss_Laguerre(x,w,alf)
    ! Given alf, the parameter á of the Laguerre polynomials, this routine returns arrays x and w
    ! of length N containing the abscissas and weights of the N-point Gauss-Laguerre quadrature
    ! formula. The abscissas are returned in ascending order. The parameter EPS is the relative
    ! precision.
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: alf
    REAL(8), DIMENSION(:), INTENT(OUT) :: x,w
    INTEGER :: its,j,n
    INTEGER, PARAMETER :: MAXIT=20
    REAL(8) :: anu
    REAL(8), PARAMETER :: C1=9.084064d-01,C2=5.214976d-02,&
         C3=2.579930d-03,C4=3.986126d-03
    REAL(8), DIMENSION(size(x)) :: rhs,r2,r3,theta
    REAL(8), DIMENSION(size(x)) :: p1,p2,p3,pp,z,z1
    LOGICAL, DIMENSION(size(x)) :: unfinished
    n=SIZE(x)
    anu=4.0d0*n+2.0d0*alf+2.0d0
    rhs=arth(4*n-1,-4,n)*NR_PI/anu
    r3=rhs**(1.0d0/3.0d0)
    r2=r3**2
    theta=r3*(C1+r2*(C2+r2*(C3+r2*C4)))
    z=anu*COS(theta)**2
    unfinished=.TRUE.
    DO its=1,MAXIT
       WHERE (unfinished)
          p1=1.0d0
          p2=0.0d0
       END WHERE
       DO j=1,n
          WHERE (unfinished)
             p3=p2
             p2=p1
             p1=((2.0d0*j-1.0d0+alf-z)*p2-(j-1.0d0+alf)*p3)/j
          END WHERE
       END DO
       WHERE (unfinished)
          pp=(n*p1-(n+alf)*p2)/z
          z1=z
          z=z1-p1/pp
          unfinished=(ABS(z-z1) > NR_EPS*z)
       END WHERE
       IF (.NOT. ANY(unfinished)) EXIT
    END DO
    IF (its == MAXIT+1) CALL NRERROR('too many iterations: Gauss Laguerre')
    x=z
    w=-EXP(LOG_GAMMA(alf+n)-LOG_GAMMA(REAL(n,8)))/(pp*n*p2)
  END SUBROUTINE Gauss_Laguerre
  
  SUBROUTINE Gauss_Hermite(x,w)
    !This routine returns arrays x and w of length N containing the abscissas and weights of
    !the N-point Gauss-Hermite quadrature formula. The abscissas are returned in descending
    !order. Note that internal computations are done in double precision.
    !Parameters: EPS is the relative precision.
    IMPLICIT NONE
    REAL(8), DIMENSION(:), INTENT(OUT) :: x,w
    REAL(8), PARAMETER :: PIM4=0.7511255444649425d0
    INTEGER :: its,j,m,n
    INTEGER, PARAMETER :: MAXIT=500
    REAL(8) :: anu
    REAL(8), PARAMETER :: C1=9.084064d-01,C2=5.214976d-02,&
         C3=2.579930d-03,C4=3.986126d-03
    REAL(8), DIMENSION((size(x)+1)/2) :: rhs,r2,r3,theta
    REAL(8), DIMENSION((size(x)+1)/2) :: p1,p2,p3,pp,z,z1
    LOGICAL, DIMENSION((size(x)+1)/2) :: unfinished
    LOGICAL :: fin
    n=SIZE(x)
    m=(n+1)/2
    anu=2.0d0*n+1.0d0
    rhs=arth(3,4,m)*NR_PI/anu
    r3=rhs**(1.0d0/3.0d0)
    r2=r3**2
    theta=r3*(C1+r2*(C2+r2*(C3+r2*C4)))
    z=SQRT(anu)*COS(theta)
    unfinished=.TRUE.
    fin=.FALSE.
    DO its=1,MAXIT
       WHERE (unfinished)
          p1=PIM4
          p2=0.0d0
       END WHERE
       DO j=1,n
          WHERE (unfinished)
             p3=p2
             p2=p1
             p1=z*SQRT(2.0d0/j)*p2-SQRT(REAL(j-1,8)/REAL(j,8))*p3
          END WHERE
       END DO
       WHERE (unfinished)
          pp=SQRT(2.0d0*n)*p2
          z1=z
          z=z1-p1/pp
          unfinished=(ABS(z-z1) > 1.0d-12)! NR_EPS)
       END WHERE
       IF (ALL(unfinished.EQV..FALSE.)) fin=.TRUE.
       IF (fin) EXIT
    END DO
    IF (its == MAXIT+1) CALL NRERROR('too many iterations: Gauss_Hermite')
    x(1:m)=z
    x(n:n-m+1:-1)=-z
    w(1:m)=2.0d0/pp**2
    w(n:n-m+1:-1)=w(1:m)
  END SUBROUTINE Gauss_Hermite

  SUBROUTINE Gauss_Jacobi(x,w,alf,bet)
    !Given alf and bet, the parameters á and â of the Jacobi polynomials, this routine returns
    !arrays x and w of length N containing the abscissas and weights of the N-point Gauss-
    !Jacobi quadrature formula. The abscissas are returned in descending order. The parameter
    !EPS is the relative precision.
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: alf,bet
    REAL(8), DIMENSION(:), INTENT(OUT) :: x,w
    INTEGER :: its,j,n
    INTEGER, PARAMETER :: MAXIT=20
    REAL(8) :: alfbet,a,c,temp
    REAL(8), DIMENSION(size(x)) :: b,p1,p2,p3,pp,z,z1
    LOGICAL, DIMENSION(size(x)) :: unfinished
    n=SIZE(x)
    alfbet=alf+bet
    z=COS(NR_PI*(arth(1,1,n)-0.25d0+0.5d0*alf)/(n+0.5d0*(alfbet+1.0d0)))
    unfinished=.TRUE.
    DO its=1,MAXIT
       temp=2.0d0+alfbet
       WHERE (unfinished)
          p1=(alf-bet+temp*z)/2.0d0
          p2=1.0d0
       END WHERE
       DO j=2,n
          a=2*j*(j+alfbet)*temp
          temp=temp+2.0d0
          c=2.0d0*(j-1.0d0+alf)*(j-1.0d0+bet)*temp
          WHERE (unfinished)
             p3=p2
             p2=p1
             b=(temp-1.0d0)*(alf*alf-bet*bet+temp*(temp-2.0d0)*z)
             p1=(b*p2-c*p3)/a
          END WHERE
       END DO
       WHERE (unfinished)
          pp=(n*(alf-bet-temp*z)*p1+2.0d0*(n+alf)*&
               (n+bet)*p2)/(temp*(1.0d0-z*z))
          z1=z
          z=z1-p1/pp
          unfinished=(ABS(z-z1) > NR_EPS)
       END WHERE
       IF (.NOT. ANY(unfinished)) EXIT
    END DO
    IF (its == MAXIT+1) CALL NRERROR('too many iterations: Gauss_Jacobi')
    x=z
    w=EXP(LOG_GAMMA(alf+n)+LOG_GAMMA(bet+n)-LOG_GAMMA(n+1.0d0)-&
         LOG_GAMMA(n+alf+bet+1.0d0))*temp*2.0d0**alfbet/(pp*p2)
  END SUBROUTINE Gauss_Jacobi
  
  SUBROUTINE Gauss_Chebyshev(x,w)
    INTEGER :: i,n
    REAL(8), INTENT(OUT) :: x(:),w(:)
    REAL(8) :: nd
    n=SIZE(x)
    nd=DBLE(n)
    DO i=1,n
       x(i)=COS( NR_PI*(DBLE(i)-0.5)/nd);
       w(i)=NR_PI/nd;
    END DO
  END SUBROUTINE Gauss_Chebyshev

  SUBROUTINE Lobatto_Gauss_Legendre(z,w,z1,z2)
    !Given the lower and upper limits of integration x1 and x2, this routine returns arrays x and w
    !of length N containing the abscissas and weights of the Lobatto-Gauss-Legendre N-point quadrature
    !formula. The parameter EPS is the relative precision.
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: z1,z2
    REAL(8), DIMENSION(:), INTENT(OUT) :: z,w
    INTEGER :: its,k,n,n1
    INTEGER, PARAMETER :: MAXIT=20
    REAL(8) :: zl,zm,xold(SIZE(z))
    REAL(8) :: x(SIZE(z)),p(SIZE(z),SIZE(z))
    n=SIZE(z)-1
    ! Truncation + 1
    n1=n+1
    ! Use the Chebyshev-Gauss-Lobatto nodes as the first guess
    x=ARTH(0,1,n1)
    x=NR_PI*x/DBLE(n)
    x=COS(x)
    ! The Legendre Vandermonde Matrix
    p=0.0d0
    ! Compute P_(N) using the recursion relation
    ! Compute its first and second derivatives and 
    ! update x using the Newton-Raphson method.
    xold=2.0d0
    its=0
    DO WHILE (MAXVAL(SQRT((x-xold)**2))>NR_EPS)
       xold=x
       P(:,1)=1.0d0
       P(:,2)=x
       DO k=2,n
          P(:,k+1)=( (DBLE(2*k-1)*x*P(:,k))-(DBLE(k-1)*P(:,k-1)) )/DBLE(k)
       END DO
       x=xold-( (x*P(:,n1))-P(:,n) ) / ( n1*P(:,n1) )
       its=its+1
    END DO
    IF (its == MAXIT) CALL NRERROR('too many iterations: Lobatto Gauss Legendre')
    zm=0.5d0*(z2+z1)
    zl=0.5d0*(z2-z1)
    w=2.0d0*zl/(DBLE(n)*DBLE(n1)*(P(:,n1)**2))
    z=zm+zl*x
  END SUBROUTINE Lobatto_Gauss_Legendre

!!!!!!!!!!!!HELPER FUNCTIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE P_TRAPZD_INT(func,a,b,s,n,i)
    ! This routine computes the nth stage of refinement of anrutil.f90 anneal.f90 integration.f90 matrix.f90 minimization.f90 nonparametric.f90 probability.f90 random.f90 statistics.f90n extended trapezoidal rule. func is
    ! input as the name of the function to be integrated between limits a and b, also input. When
    ! called with n=1, the routine returns as s the crudest estimate of the integral. Subsequent
    ! calls with n=2,3,... (in that sequential order) will improve the accuracy of s by adding 2n-2
    ! additional interior points. s should not be modified between sequential calls.
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i
    REAL(8), INTENT(IN) :: a,b
    REAL(8), INTENT(INOUT) :: s
    INTEGER, INTENT(IN) :: n
    INTERFACE
       REAL(8) FUNCTION FUNC(x,i)
         INTEGER, INTENT(IN) :: i
         REAL(8), INTENT(IN) :: x
       END FUNCTION FUNC
    END INTERFACE
    REAL(8) :: del,fsum
    REAL(8), ALLOCATABLE :: ar(:)
    INTEGER :: it,j
    IF (n==1) THEN
       fsum=FUNC(a,i)
       fsum=fsum+FUNC(b,i)
       s=0.5d0*(b-a)*fsum
    ELSE
       it=2**(n-2)
       del=(b-a)/it
       ALLOCATE(ar(it))
       ar=ARTH(a+0.5d0*del,del,it)
       fsum=FUNC(ar(1),i)
       DO j = 2, it
          fsum=fsum+FUNC(ar(j),i)
       END DO
       DEALLOCATE(ar)
       s=0.5d0*(s+del*fsum)
    END IF
  END SUBROUTINE P_TRAPZD_INT
  
  SUBROUTINE P_POLINT_INT(xa,ya,x,y,dy)
    IMPLICIT NONE
    REAL(8), DIMENSION(:), INTENT(IN) :: xa,ya
    REAL(8), INTENT(IN) :: x
    REAL(8), INTENT(OUT) :: y,dy
    INTEGER :: m,n,ns
    REAL, DIMENSION(SIZE(xa)) :: c,d,den,ho
    c=ya
    d=ya
    ho=xa-x
    ns=IMINLOC(ABS(x-xa))
    y=ya(ns)
    ns=ns-1
    n=SIZE(xa)
    DO m=1,n-1
       den(1:n-m)=ho(1:n-m)-ho(1+m:n)
       IF (ANY(den(1:n-m) == 0.0)) CALL NRERROR('calculation failure: P_POLINT_INT')
       den(1:n-m)=(c(2:n-m+1)-d(1:n-m))/den(1:n-m)
       d(1:n-m)=ho(1+m:n)*den(1:n-m)
       c(1:n-m)=ho(1:n-m)*den(1:n-m)
       IF (2*ns < n-m) THEN
          dy=c(ns+1)
       ELSE
          dy=d(ns)
          ns=ns-1
       END IF
       y=y+dy
    END DO
  END SUBROUTINE P_POLINT_INT
 
END MODULE INTEGRATION
