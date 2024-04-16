MODULE RANDOM
  USE NRUTIL, ONLY : NRERROR,NR_SMALL,NR_BIG,NR_EPS,MMUL
  USE LAPACK95, ONLY : POTRF
  IMPLICIT NONE

  PRIVATE :: P_Cholesky_Ran,P_CDF_Normal_Ran,P_CDF_Normal_Inverse_Ran

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                         PARAMETERS FOR:
  ! Marsaglia & Tsang generator for random normals & random exponentials.
  ! Translated from C by Alan Miller (amiller@bigpond.net.au)
  ! Marsaglia, G. & Tsang, W.W. (2000) `The ziggurat method for generating
  ! random variables', J. Statist. Software, v5(8).
  ! This is an electronic journal which can be downloaded from:
  ! http://www.jstatsoft.org/v05/i08
  ! N.B. It is assumed that all integers are 32-bit.
  ! N.B. The value of zigm2 has been halved to compensate for the lack of
  !      unsigned integers in Fortran.
  ! Latest version - February 3 2004
  ! Dec 18, 2007 changed the uniform random number generator to use KISS instead
  
  INTEGER(4), SAVE :: kiss1=123456789,kiss2=362436069,kiss3=521288629,kiss4=916191069
  REAL(8), PARAMETER  :: zigm1=2147483648.0d0,zigm2=2147483648.0d0,half=0.5d0
  REAL(8) :: zigdn=3.442619855899d0,zigtn=3.442619855899d0,zigvn=0.00991256303526217d0, &
       zigq,zigde=7.697117470131487d0,zigte=7.697117470131487d0,zigve=0.003949659822581572d0
  INTEGER(4), SAVE :: zigiz,zigjz,seed,zigkn(0:127),zigke(0:255),zighz
  REAL(8), SAVE :: zigwn(0:127),zigfn(0:127),zigwe(0:255),zigfe(0:255)
  LOGICAL, SAVE :: initialized=.FALSE.

CONTAINS

  SUBROUTINE Set_Seed(k1,k2,k3,k4)
    IMPLICIT NONE
    INTEGER(4), INTENT(IN), OPTIONAL :: k1,k2,k3,k4
    INTEGER :: i
    IF (PRESENT(k1)) THEN
       kiss1=k1
       kiss2=k2
       kiss3=k3
       kiss4=k4
    ELSE
       CALL SYSTEM_CLOCK(seed)
       kiss1=SHR3()
       kiss2=SHR3()
       kiss3=SHR3()
       kiss4=SHR3()
    END IF
    zigdn=3.442619855899d0
    zigtn=3.442619855899d0
    zigvn=0.00991256303526217d0
    zigde=7.697117470131487d0
    zigte=7.697117470131487d0
    zigve=0.003949659822581572d0

    !  Tables for NORMAL
    zigq = zigvn*EXP(half*zigdn*zigdn)
    zigkn(0) = (zigdn/zigq)*zigm1
    zigkn(1) = 0
    zigwn(0) = zigq/zigm1
    zigwn(127) = zigdn/zigm1
    zigfn(0) = 1.0d0
    zigfn(127) = EXP( -half*zigdn*zigdn )
    DO  i = 126, 1, -1
       zigdn = SQRT( -2.0d0 * LOG( zigvn/zigdn + EXP( -half*zigdn*zigdn ) ) )
       zigkn(i+1) = (zigdn/zigtn)*zigm1
       zigtn = zigdn
       zigfn(i) = EXP(-half*zigdn*zigdn)
       zigwn(i) = zigdn/zigm1
    END DO
    !  Tables for EXPONENTIAL
    zigq = zigve*EXP( zigde )
    zigke(0) = (zigde/zigq)*zigm2
    zigke(1) = 0
    zigwe(0) = zigq/zigm2
    zigwe(255) = zigde/zigm2
    zigfe(0) = 1.0d0
    zigfe(255) = EXP( -zigde )
    DO  i = 254, 1, -1
       zigde = -LOG( zigve/zigde + EXP( -zigde ) )
       zigke(i+1) = zigm2 * (zigde/zigte)
       zigte = zigde
       zigfe(i) = EXP( -zigde )
       zigwe(i) = zigde/zigm2
    END DO
    initialized = .TRUE.
  END SUBROUTINE Set_Seed

  ! Generate random 32-bit integers
  INTEGER(4) FUNCTION shr3()
    zigjz = seed
    seed = IEOR(seed,ISHFT(seed,13))
    seed = IEOR(seed,ISHFT(seed,-17))
    seed = IEOR(seed,ISHFT(seed,5))
    shr3 = zigjz + seed
  END FUNCTION shr3

  INTEGER(4) FUNCTION KISS ()
    ! The  KISS (Keep It Simple Stupid) random number generator. Combines:
    ! (1) The congruential generator x(n)=69069*x(n-1)+1327217885, period 2^32.
    ! (2) A 3-shift shift-register generator, period 2^32-1,
    ! (3) Two 16-bit multiply-with-carry generators, period 597273182964842497>2^59
    !  Overall period>2^123;
    IF (.NOT.INITIALIZED) CALL SET_SEED()
    kiss1 = 69069 * kiss1 + 1327217885
    kiss2 = m (m (m (kiss2, 13), - 17), 5)
    kiss3 = 18000 * iand (kiss3, 65535) + ishft (kiss3, - 16)
    kiss4 = 30903 * iand (kiss4, 65535) + ishft (kiss4, - 16)
    kiss = kiss1 + kiss2 + ishft (kiss3, 16) + kiss4
  CONTAINS
    INTEGER FUNCTION m(k, n)
      INTEGER :: k, n
      m = ieor (k, ishft (k, n) )
    END FUNCTION m
  END FUNCTION KISS
  
  REAL(8) FUNCTION Sample_Uniform(a,b)
    !INTEGER(4) ::  ival,zigjz
    REAL(8), INTENT(IN) :: a,b
    IF (a>b) CALL NRERROR('upper limit lower than lower limit: Sample_Uniform')
    Sample_Uniform = (DBLE(KISS())+2147483649.0d0)/4294967297.0d0
    Sample_Uniform = Sample_Uniform*(b-a) + a
  END FUNCTION Sample_Uniform

  !  Generate random normals
  REAL(8) FUNCTION Sample_Normal(mu,var)
    REAL(8), INTENT(IN) :: mu,var
    REAL(8), PARAMETER ::  r = 3.442620d0
    REAL(8) :: x,y
    IF (var<0.0d0) CALL NRERROR('variance has to be positive: Sample_Normal')
    zighz = KISS()
    zigiz = IAND(zighz,127)
    IF(ABS(zighz) < zigkn(zigiz)) THEN
       Sample_Normal = zighz * zigwn(zigiz)
       Sample_Normal = Sample_Normal*SQRT(var) + mu
    ELSE
       DO
          IF(zigiz==0) THEN
             DO
                x = -0.2904764d0*LOG(Sample_Uniform(0.0d0,1.0d0))
                y = -LOG(Sample_Uniform(0.0d0,1.0d0))
                IF( y+y >= x*x ) EXIT
             END DO
             Sample_Normal = r+x			
             IF (zighz<=0) Sample_Normal = -Sample_Normal
             Sample_Normal = Sample_Normal*SQRT(var) + mu
             RETURN
          END IF
          x = zighz * zigwn(zigiz)
          IF (zigfn(zigiz) + Sample_Uniform(0.0d0,1.0d0)* &
               (zigfn(zigiz-1)-zigfn(zigiz))< EXP(-half*x*x)) THEN
             Sample_Normal = x
             Sample_Normal = Sample_Normal*SQRT(var) + mu
             RETURN
          END IF
          zighz = KISS()
          zigiz = IAND(zighz,127)
          IF(ABS(zighz) < zigkn(zigiz)) THEN
             Sample_Normal = zighz * zigwn(zigiz)
             Sample_Normal = Sample_Normal*SQRT(var) + mu
             RETURN
          END IF
       END DO
    END IF
  END FUNCTION Sample_Normal

  !  Generate random exponentials
  REAL(8) FUNCTION Sample_Exponential(a)
    REAL(8), INTENT(IN) :: a
    REAL(8)  ::  x
    IF (a<0.0d0) CALL NRERROR('a has to be positive: Sample_Exponential')
    zigjz = KISS()
    zigiz = IAND(zigjz,255)
    IF (ABS(zigjz) < zigke(zigiz)) THEN
       Sample_Exponential = (ABS(zigjz) * zigwe(zigiz))*a
       RETURN
    END IF
    DO
       IF (zigiz==0) THEN
          Sample_Exponential = (7.69711 - LOG(Sample_Uniform(0.0d0,1.0d0)))*a
          RETURN
       END IF
       x = ABS(zigjz) * zigwe(zigiz)
       IF (zigfe(zigiz) + Sample_Uniform(0.0d0,1.0d0)*(zigfe(zigiz-1) - zigfe(zigiz)) < EXP(-x)) THEN
          Sample_Exponential = x*a
          RETURN
       END IF
       zigjz = KISS()
       zigiz = IAND(zigjz,255)
       IF (ABS(zigjz) < zigke(zigiz)) THEN
          Sample_Exponential = (ABS(zigjz)*zigwe(zigiz))*a
          RETURN
       END IF
    END DO
  END FUNCTION Sample_Exponential

  REAL(8) FUNCTION Sample_Gamma(ain,b)
    ! Returns a number distributed as a gamma distribution with shape parameter a and inverse scale b.
    ! such that mean=a/b and var=a/(b*b)
    ! That is P(x) = [(b^a)/Gamma(a)] * [x^(a-1)] * exp(-x*b) for x>0 and a>=1.
    ! Uses the algorithm in
    ! Marsaglia, G. and Tsang, W.W. (2000) `A simple method for generating
    ! gamma variables', Trans. om Math. Software (TOMS), vol.26(3), pp.363-372.
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: ain,b
    REAL(8) :: c,d,u,v,x,correction,a
    correction=1.0d0
    a=ain
    IF (ain<0.0d0) CALL NRERROR('a has to be positive: Sample_Gamma')
    IF (b<0.0d0) CALL NRERROR('b has to be positive: Sample_Gamma')
    IF (a<1.0d0) THEN
       correction=Sample_Uniform(0.0d0,1.0d0)**(1.0d0/ain)
       a = ain + 1.0
    END IF
    d = a - 1.0d0/3.0d0
    c = 1.0d0/SQRT(9.0d0*d)
    ! Start of main loop
    DO
       ! Generate v = (1+cx)^3 where x is random normal; repeat if v <= 0.
       DO
          x = Sample_Normal(0.0d0,1.0d0)
          v = (1.0d0 + c*x)**3
          IF (v > 0.0d0) EXIT
       END DO
       ! Generate uniform variable U
       u = Sample_Uniform(0.0d0,1.0d0)
       IF (u < 1.0d0 - 0.0331d0*x**4.0d0) THEN
          Sample_Gamma = correction*d*v
          EXIT
       ELSE IF (LOG(u) < 0.5d0*(x**2.0d0) + d*(1.0d0 - v + LOG(v))) THEN
          Sample_Gamma = correction*d*v
          EXIT
       END IF
    END DO
    Sample_Gamma=Sample_Gamma/b
  END FUNCTION Sample_Gamma

  REAL(8) FUNCTION Sample_Chi2(df)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: df
    Sample_Chi2=2.0d0*Sample_Gamma(0.5d0*df,1.0d0)
  END FUNCTION Sample_Chi2
       
  REAL(8) FUNCTION Sample_Mixture_Gamma(a,b,p)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: a(:),b(:),p(:)
    REAL(8) :: u,cdf,x
    INTEGER :: j
    cdf=0.0d0
    u=Sample_Uniform(0.0d0,1.0d0)
    DO j=1,SIZE(p)
       cdf=cdf+p(j)
       x=Sample_Gamma(a(j),b(j))
       IF (u<=cdf) THEN
          Sample_Mixture_GAMMA=x
          RETURN
       END IF
    END DO
  END FUNCTION Sample_Mixture_Gamma
  
  REAL(8) FUNCTION Sample_Mixture_Normal(mu,var,p)
    ! Returns 1 draws from a mixture of normals with nmix mixture components
    ! p vector of weights, mu vector of means and sigma vector of variances
    IMPLICIT NONE
    REAL(8), INTENT(in) :: p(:),mu(:),var(:)
    INTEGER :: i, j,nmix
    REAL(8) :: cdf(SIZE(p)+1),nor,u
    IF (ANY(var<0.0d0)) CALL NRERROR('variance has to be positive: Sample_Mixture_Normal')
    IF ((SUM(p)>1.0000001d0).OR.(SUM(p)<0.9999999d0)) CALL NRERROR('weights have to add to one: Sample_Mixture_Normal')
    IF (ANY(p<0.0d0)) CALL NRERROR('weights have to be positive: Sample_Mixture_Normal')
    nmix=SIZE(p)
    u = Sample_Uniform(0.0d0,1.0d0)
    nor = Sample_Normal(0.0d0,1.0d0)
    cdf=0.0D0
    DO i = 1, nmix
       cdf(i+1) = cdf(i) + p(i)
       IF ((u>cdf(i)).AND.(u<=cdf(i+1))) THEN
          Sample_Mixture_Normal=mu(i)+SQRT(var(i))*nor
          RETURN
       END IF
    END DO
  END FUNCTION Sample_Mixture_Normal

  REAL(8) FUNCTION Sample_Truncated_Mixture_Normal(mu,var,p,a,lb)
    ! Returns 1 draws from a mixture of normals with nmix mixture components
    ! p vector of weights, mu vector of means and sigma vector of variances
    ! that is truncated to be between a and infinity (lb=true) or -infinity and a (lb=false)
    ! !!!!!!!!!!!!!!!!!!!!!!!!!NOT the same as sampling from a !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!mixture of truncated normals!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! the weights need to be updated
    IMPLICIT NONE
    REAL(8), INTENT(in) :: p(:),mu(:),var(:),a
    LOGICAL, INTENT(IN) :: lb
    INTEGER :: i, j,nmix
    REAL(8) :: cdf(SIZE(p)+1),u,w(SIZE(p))
    IF (ANY(var<0.0d0)) CALL NRERROR('variance has to be positive: Sample_Truncated_Mixture_Normal')
    IF ((SUM(p)>1.0000001d0).OR.(SUM(p)<0.9999999d0)) CALL NRERROR('weights have to add to one: Sample_Truncated_Mixture_Normal')
    IF (ANY(p<0.0d0)) CALL NRERROR('weights have to be positive: Sample_Truncated_Mixture_Normal')
    nmix=SIZE(p)
    IF (lb) THEN
       DO i = 1, nmix
          w(i) = (1.0d0 - P_CDF_Normal_Ran(a,mu(i),var(i)))*p(i)
       END DO
    ELSE
       DO i = 1, nmix
          w(i) = P_CDF_Normal_Ran(a,mu(i),var(i))*p(i)
       END DO
    END IF
    w=w/SUM(w)
    u = Sample_Uniform(0.0d0,1.0d0)
    cdf=0.0D0
    DO i = 1, nmix
       cdf(i+1) = cdf(i) + w(i)
       IF ((u>cdf(i)).AND.(u<=cdf(i+1))) THEN
          Sample_Truncated_Mixture_Normal=Sample_Truncated_Normal_Geweke(mu(i),var(i),a,lb)
          RETURN
       END IF
    END DO
  END FUNCTION Sample_Truncated_Mixture_Normal

  REAL(8) FUNCTION Sample_Truncated_Normal_Geweke(mu,var,a,lb)
    ! Returns one draw from a truncated normal with underlying
    ! mean=mu and VARIANCE=var with truncation
    !           (a,+infty)    IF lb=TRUE
    !           (-infty,a)    IF lb=FALSE
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: mu,var,a
    LOGICAL, INTENT(IN) :: lb
    REAL(8), PARAMETER :: t4=0.45D0
    REAL(8) :: u,z,phi_z,c,temp
    IF (var<0.0d0) CALL NRERROR('variance has to be positive: Sample_Truncated_Normal')
    c=((a-mu)/(SQRT(var)))
    IF (.NOT.lb) THEN
       c=-c
    END IF
    IF (c < t4) THEN
       ! normal rejection sampling
       DO
          u=Sample_Normal(0.0d0,1.0d0)
          IF (u > c) EXIT
       ENDDO
       temp=u
    ELSE
       ! exponential rejection sampling
       DO
          u = Sample_Uniform(0.0d0,1.0d0)
          z = Sample_Exponential(1.0d0/c)
          phi_z=EXP(-.5D0*(z*z))
          IF (u < phi_z) EXIT
       ENDDO
       temp=c + z
    ENDIF
    IF (.not.(lb)) THEN
       Sample_Truncated_Normal_Geweke = mu - (temp*SQRT(var))
    ELSE
       Sample_Truncated_Normal_Geweke = mu + (temp*SQRT(var))
    ENDIF
  END FUNCTION Sample_Truncated_Normal_Geweke

  REAL(8) FUNCTION Sample_Truncated_Normal(mu,var,a,lb)
    ! Returns one draw from a truncated normal with underlying
    ! mean=mu and VARIANCE=var with truncation
    !           (a,+infty)    IF lb=TRUE
    !           (-infty,a)    IF lb=FALSE
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: mu,var,a
    LOGICAL, INTENT(IN) :: lb
    REAL(8), PARAMETER :: t4=0.45D0
    REAL(8) :: u,cdf,aux
    IF (var<0.0d0) CALL NRERROR('variance has to be positive: Sample_Truncated_Normal')
    u = Sample_Uniform(0.0d0,1.0d0) 
    cdf=P_CDF_Normal_Ran(a,mu,var)
    IF (lb) THEN
       aux = (u*(1.0-cdf)) + cdf       
    ELSE
       aux = u*cdf
    END IF
    Sample_Truncated_Normal = P_CDF_Normal_Inverse_Ran(aux,mu,var)
  END FUNCTION Sample_Truncated_Normal
  
  REAL(8) FUNCTION Sample_Double_Truncated_Mixture_Normal(mu,var,p,a,b)
    ! Returns 1 draws from a mixture of normals with nmix mixture components
    ! p vector of weights, mu vector of means and sigma vector of variances
    ! that is truncated to be between a and b
    ! !!!!!!!!!!!!!!!!!!!!!!!!!NOT the same as sampling from a !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!mixture of truncated normals!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! the weights need to be updated
    IMPLICIT NONE
    REAL(8), INTENT(in) :: p(:),mu(:),var(:),a,b
    INTEGER :: i, j,nmix
    REAL(8) :: cdf(SIZE(p)+1),u,w(SIZE(p))
    IF (ANY(var<0.0d0)) CALL NRERROR('variance has to be positive: Sample_Double_Truncated_Mixture_Normal')
    IF ((SUM(p)>1.0000001d0).OR.(SUM(p)<0.9999999d0)) & 
         CALL NRERROR('weights have to add to one: Sample_Double_Truncated_Mixture_Normal')
    IF (ANY(p<0.0d0)) CALL NRERROR('weights have to be positive: Sample_Double_Truncated_Mixture_Normal')
    nmix=SIZE(p)
    DO i = 1, nmix
       w(i) = (P_CDF_Normal_Ran(b,mu(i),var(i)) - P_CDF_Normal_Ran(a,mu(i),var(i)))*p(i)
    END DO
    w=w/SUM(w)
    u = Sample_Uniform(0.0d0,1.0d0)
    cdf=0.0D0
    DO i = 1, nmix
       cdf(i+1) = cdf(i) + w(i)
       IF ((u>cdf(i)).AND.(u<=cdf(i+1))) THEN
          Sample_Double_Truncated_Mixture_Normal=Sample_Double_Truncated_Normal(mu(i),var(i),a,b)
          RETURN
       END IF
    END DO
  END FUNCTION Sample_Double_Truncated_Mixture_Normal

  REAL(8) FUNCTION Sample_Double_Truncated_Normal_Geweke(mu,var,a,b)
    ! Generates one draw from truncated standard normal distribution (mu,sigma) on (a,b)
    IMPLICIT NONE
    REAL(8),INTENT(in) :: a,b,mu,var
    REAL(8) :: c,c1,c2,u(2),x,cdel,f1,f2,z,az,bz,eps
    REAL(8),PARAMETER :: t1=0.375D0,t2=2.18D0,t3=0.725D0,t4=0.45D0
    LOGICAL :: lflip
    REAL :: aaa
    INTEGER :: j
    eps=2.220446049250313D-016
    az=(a-mu)/SQRT(var)
    bz=(b-mu)/SQRT(var)
    c1=az
    c2=bz
    lflip=.false.
    IF (c1*c2<0.0D0) THEN
       IF ((f(c1)>t1) .and. (f(c2)>t1)) THEN
          cdel=c2-c1
          DO
             DO j = 1, 2
                u(j) = Sample_Uniform(0.0d0,1.0d0)
             END DO
             x=c1+cdel*u(1)
             IF (u(2)<f(x)) EXIT
          END DO
       ELSE
          DO
             x=Sample_Normal(0.0d0,1.0d0)
             IF ((x>c1) .and. (x<c2)) EXIT
          END DO
       END IF
    ELSE
       IF (c1<0.0D0) THEN
          c=c1
          c1=-c2
          c2=-c
          lflip=.true.
       END IF
       f1=f(c1)
       f2=f(c2)
       IF ((f2<eps) .or. (f1/f2>t2)) THEN
          IF (c1>t3) THEN
             !exponential rejection sampling
             c=c2-c1
             DO
                u(1) = Sample_Uniform(0.0d0,1.0d0)
                z = Sample_Exponential(1.0d0/c1)
                IF ((z<c) .and. (u(1)<f(z))) EXIT
             END DO
             x=c1+z
          ELSE
             !half-normal rejection sampling
             DO
                x=Sample_Normal(0.0d0,1.0d0)
                x=abs(x)
                IF ((x>c1) .and. (x<c2)) EXIT
             END DO
          END IF
       ELSE
          !uniform rejection sampling
          cdel=c2-c1
          DO
             DO j = 1, 2
                u(j) = Sample_Uniform(0.0d0,1.0d0)
             END DO
             x=c1+cdel*u(1)
             IF (u(2)<(f(x)/f1)) EXIT
          END DO
       END IF
    END IF
    IF (lflip) THEN
       Sample_Double_Truncated_Normal_Geweke=mu-(SQRT(var)*x)
    ELSE
       Sample_Double_Truncated_Normal_Geweke=mu+(SQRT(var)*x)
    END IF
  CONTAINS
    REAL(8) FUNCTION f(x)
      REAL(8) :: x
      f=dexp(-.5D0*(x*x))
    END FUNCTION f
  END FUNCTION Sample_Double_Truncated_Normal_Geweke

  REAL(8) FUNCTION Sample_Double_Truncated_Normal(mu,var,a,b)
    ! Returns one draw from a truncated normal with underlying
    ! mean=mu and VARIANCE=var with truncation
    !           (a,+infty)    IF lb=TRUE
    !           (-infty,a)    IF lb=FALSE
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: mu,var,a,b
    REAL(8) :: u,cdfa,cdfb,aux
    IF (var<0.0d0) CALL NRERROR('variance has to be positive: Sample_Double_Truncated_Normal')
    u = Sample_Uniform(0.0d0,1.0d0) 
    cdfa=P_CDF_Normal_Ran(a,mu,var)
    cdfb=P_CDF_Normal_Ran(b,mu,var)
    aux = (u*(cdfb-cdfa))+cdfa
    Sample_Double_Truncated_Normal = P_CDF_Normal_Inverse_Ran(aux,mu,var)
  END FUNCTION Sample_Double_Truncated_Normal

  FUNCTION Sample_Dirichlet(k,a)
    ! DIRICHLET
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: k
    REAL(8), DIMENSION(:), INTENT(IN) :: a
    REAL(8) :: Sample_Dirichlet(k)
    REAL(8) :: rg(k),sg
    INTEGER :: i
    DO i = 1, k
       rg(i) = Sample_Gamma(a(i),1.0d0)
    ENDDO
    sg=SUM(rg)
    Sample_Dirichlet = rg/sg
  END FUNCTION Sample_Dirichlet

  INTEGER FUNCTION Sample_Multinomial(p)
    ! INPUT  :   p  is kx1
    ! Returns a random variable sampled from a multinomial distribution with
    ! k categories having probabilities p
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: p(:)
    REAL(8) :: u,cdf
    INTEGER :: j,k
!    IF (ABS(SUM(p)-1.0d0)<1.0d-8) THEN
!       WRITE(*,*) p
!       CALL NRERROR("p does not sum to 1")
!    END IF
    k=SIZE(p)
    u = Sample_Uniform(0.0d0,1.0d0)
    cdf=0.0d0
    DO j=1,k-1
       IF ((u>cdf).AND.(u<=cdf+p(j))) THEN
          Sample_Multinomial=j
          RETURN
       END IF
       cdf=cdf+p(j)
    END DO
    IF (u>cdf) Sample_Multinomial=k
  END FUNCTION Sample_Multinomial

  FUNCTION Sample_Multivariate_Normal(mu,var)
    ! Returns one draw from a multivariate normal with mean mu and varcovar var
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: mu(:),var(:,:)
    REAL(8) :: cvar(SIZE(var,1),SIZE(var,2)),Sample_Multivariate_Normal(SIZE(mu))
    INTEGER :: j
    cvar=P_Cholesky_Ran(var)
    DO j = 1, SIZE(mu)
       Sample_Multivariate_Normal(j) = Sample_Normal(0.0d0,1.0d0)
    END DO
    Sample_Multivariate_Normal = mu + MMUL(cvar,Sample_Multivariate_Normal)
  END FUNCTION Sample_Multivariate_Normal

  REAL(8) FUNCTION Sample_Logit()
    IMPLICIT NONE
    REAL(8) :: u
    u=Sample_Uniform(0.0d0,1.0d0)
    Sample_Logit=LOG((1.0d0-u)/u)
  END FUNCTION Sample_Logit

  REAL(8) FUNCTION Sample_EV1()
    REAL(8) :: u
    u=Sample_Uniform(0.0d0,1.0d0)
    Sample_EV1=-LOG(-LOG(u))
  END FUNCTION Sample_EV1

  FUNCTION P_Cholesky_RAN(AIN)
    IMPLICIT NONE
    REAL(8), DIMENSION(:,:), INTENT(IN) :: AIN
    REAL(8) :: A(SIZE(AIN,1),SIZE(AIN,2)),P_Cholesky_RAN(SIZE(AIN,1),SIZE(AIN,2))
    INTEGER :: i,j,n
    n=SIZE(ain,1)
    A=AIN
    CALL POTRF(A,'L')
    P_Cholesky_RAN=0.0d0
    FORALL(i=1:n,j=1:n,j<=i) P_Cholesky_RAN(i,j)=A(i,j)
  END FUNCTION P_Cholesky_RAN

  REAL(8) FUNCTION P_CDF_Normal_Ran(zin,mu,var)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: zin
    REAL(8), INTENT(IN) :: mu,var
    REAL(8) :: zabs,p,arg,logpdf,z,std
    REAL(8), PARAMETER :: p0=220.2068679123761D0,p1=221.2135961699311D0,p2=112.0792914978709D0, &
         p3 = 33.91286607838300D0,p4 = 6.373962203531650D0,p5 = .7003830644436881D0, &
         p6 = .3526249659989109D-01,q0 = 440.4137358247522D0,q1 = 793.8265125199484D0, &
         q2 = 637.3336333788311D0,q3 = 296.5642487796737D0,q4 = 86.78073220294608D0, &
         q5=16.06417757920695D0,q6=1.755667163182642D0,q7=.8838834764831844D-1,cutoff = 7.071D0, &
         logroot2pi = 0.918938533204672780563271317078D0
    IF (var<0.0d0) CALL NRERROR('variance has to be positive: P_CDF_Normal_Ran')
    z=zin-mu
    zabs=ABS(z)  
    IF (zabs<NR_Small) THEN
       P_CDF_Normal_Ran=0.5d0
       RETURN
    END IF
    std = SQRT(var)
    IF (std<NR_Small) THEN
       IF (zin-mu>0.0d0) THEN
          P_CDF_Normal_Ran = 1.0d0
       ELSE IF	(zin-mu<0.0d0) THEN
          P_CDF_Normal_Ran = 0.0d0
       END IF
    END IF
    zabs=zabs/std
    IF (z > 37.0D0) THEN
       P_CDF_Normal_Ran = 1.0D0
       RETURN
    ELSE IF (z < -37.0D0) THEN
       P_CDF_Normal_Ran = 0.0D0
       RETURN
    END IF
    arg = -0.5D0*zabs*zabs
    logpdf = -logroot2pi - LOG(std) + arg
    IF (zabs < cutoff) THEN
       p = arg + LOG(((((((p6*zabs + p5)*zabs + p4)*zabs + p3)*zabs + &
            p2)*zabs + p1)*zabs + p0)) - LOG((((((((q7*zabs + q6)*zabs + &
            q5)*zabs + q4)*zabs + q3)*zabs + q2)*zabs + q1)*zabs + &
            q0))
    ELSE
       p = logpdf - LOG((zabs + 1.0D0/(zabs + 2.0D0/(zabs + 3.0D0/(zabs + 4.0D0/ &
            (zabs + 0.65D0))))))
    END IF
    p = EXP(p)
    IF (z < 0.0D0) THEN
       P_CDF_Normal_Ran=p
       RETURN
    ELSE
       P_CDF_Normal_Ran = 1.0D0 - p
       RETURN
    END IF
    RETURN
  END FUNCTION P_CDF_Normal_Ran
  
  REAL(8) FUNCTION P_CDF_Normal_Inverse_Ran(P,muin,varin)
    !	Produces the normal deviate Z corresponding to a given lower
    !	tail area of P; Z is accurate to about 1 part in 10**16.
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: P
    REAL(8), OPTIONAL, INTENT(IN) :: muin,varin
    REAL(8) :: Q,R,mu,var
    REAL(8), PARAMETER :: ZERO=0.D0, ONE = 1.D0, HALF = 0.5D0,SPLIT1 = 0.425D0, SPLIT2 = 5.D0, &
         CONST1 = 0.180625D0, CONST2 = 1.6D0
    !	Coefficients for P close to 0.5
    REAL(8), PARAMETER :: A0 = 3.3871328727963666080D0,A1=1.3314166789178437745D+2, &
         A2=1.9715909503065514427D+3,A3=1.3731693765509461125D+4,A4=4.5921953931549871457D+4, &
         A5=6.7265770927008700853D+4,A6=3.3430575583588128105D+4,A7=2.5090809287301226727D+3, &
         B1=4.2313330701600911252D+1,B2=6.8718700749205790830D+2,B3=5.3941960214247511077D+3, &
         B4=2.1213794301586595867D+4,B5=3.9307895800092710610D+4,B6=2.8729085735721942674D+4, &
         B7=5.2264952788528545610D+3
    !	Coefficients for P not close to 0, 0.5 or 1.
    REAL(8), PARAMETER :: C0=1.42343711074968357734D0,C1=4.63033784615654529590D0, &
         C2=5.76949722146069140550D0,C3=3.64784832476320460504D0,C4=1.27045825245236838258D0, &
         C5=2.41780725177450611770D-1,C6=2.27238449892691845833D-2,C7=7.74545014278341407640D-4, &
         D1=2.05319162663775882187D0,D2=1.67638483018380384940D0,D3=6.89767334985100004550D-1, &
         D4=1.48103976427480074590D-1,D5=1.51986665636164571966D-2,D6=5.47593808499534494600D-4, &
         D7=1.05075007164441684324D-9
    !	CoefficientsforP near 0 or 1.
    REAL(8), PARAMETER :: E0=6.65790464350110377720D0,E1=5.46378491116411436990D0, &
         E2=1.78482653991729133580D0,E3=2.96560571828504891230D-1,E4=2.65321895265761230930D-2, &
         E5=1.24266094738807843860D-3,E6=2.71155556874348757815D-5,E7=2.01033439929228813265D-7, &
         F1=5.99832206555887937690D-1,F2=1.36929880922735805310D-1,F3=1.48753612908506148525D-2, &
         F4=7.86869131145613259100D-4,F5=1.84631831751005468180D-5,F6=1.42151175831644588870D-7, &
         F7=2.04426310338993978564D-15
    IF ((P>1.0D0).OR.(P<0.0D0)) CALL NRERROR('P has to be between zero and one: P_CDF_Normal_Inverse_Ran')
    mu=0.0d0
    IF (PRESENT(muin)) mu=muin
    var=1.0d0
    IF (PRESENT(varin)) THEN
       IF (varin<0.0d0) CALL NRERROR('variance has to be positive: P_CDF_Normal_Inverse_Ran')
       var=varin
    END IF
    Q = P - HALF
    IF (ABS(Q) .LE. SPLIT1) THEN
       R = CONST1 - Q * Q
       P_CDF_Normal_Inverse_Ran=Q*(((((((A7*R+A6)*R+A5)*R+A4)*R+A3)*R+A2)*R+A1)*R+ A0) / &
            (((((((B7 * R + B6) * R + B5) * R + B4) * R + B3)* R + B2) * R + B1) * R + ONE)
       P_CDF_Normal_Inverse_Ran=(SQRT(var)*P_CDF_Normal_Inverse_Ran)+mu
       RETURN
    ELSE
       IF (Q .LT. ZERO) THEN
          R = P
       ELSE
          R = ONE - P
       END IF
       R = SQRT(-LOG(R))
       IF (R .LE. SPLIT2) THEN
          R = R - CONST2
          P_CDF_Normal_Inverse_Ran=(((((((C7*R + C6) * R + C5) * R + C4) * R + C3)*R+C2)*R+ C1) * R + C0) / &
               (((((((D7 * R + D6) * R + D5) * R + D4) * R + D3)* R + D2) * R + D1) * R + ONE)
       ELSE
          R = R - SPLIT2
          P_CDF_Normal_Inverse_Ran=(((((((E7*R + E6) * R + E5) * R + E4) * R + E3)*R+E2)*R + E1) * R + E0) / &
               (((((((F7 * R + F6) * R + F5) * R + F4) * R + F3)* R + F2) * R + F1) * R + ONE)
       END IF
       IF (Q .LT. ZERO) P_CDF_Normal_Inverse_Ran = - P_CDF_Normal_Inverse_Ran
       P_CDF_Normal_Inverse_Ran=(SQRT(var)*P_CDF_Normal_Inverse_Ran)+mu
       RETURN
    END IF
  END FUNCTION P_CDF_Normal_Inverse_Ran
  
END MODULE RANDOM
