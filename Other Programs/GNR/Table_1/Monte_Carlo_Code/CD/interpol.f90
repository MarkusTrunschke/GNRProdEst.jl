MODULE INTERPOL
  USE NRUTIL, ONLY : ARRAY1D
  IMPLICIT NONE
  
  !PRIVATE :: PRIV_INTA1,PRIV_CHEB_VAR,PRIV_CHEB_COEF,
  PRIVATE :: halfpi,P_Sort_Interp,P_Inssort_Interp ! Helper function, same as in matrix module

  REAL(8), PARAMETER :: halfpi=1.57079632679489d0

  INTERFACE CHEB_INIT
     MODULE PROCEDURE INIT_CHEB_1,INIT_CHEB_2
  END INTERFACE
  INTERFACE CHEB_INTER
     MODULE PROCEDURE CHEB_INTER_1,CHEB_INTER_2
  END INTERFACE
  INTERFACE PLINEAR_INIT
     MODULE PROCEDURE PLINEAR_INIT_1,PLINEAR_INIT_2,PLINEAR_INIT_3,PLINEAR_INIT_4
  END INTERFACE
  INTERFACE MLINEAR_INIT
     MODULE PROCEDURE MLINEAR_INIT_1,MLINEAR_INIT_2
  END INTERFACE
  INTERFACE PLINEAR_INTER
     MODULE PROCEDURE PLINEAR_INTER_1,PLINEAR_INTER_2
  END INTERFACE
  
  TYPE PRIV_INTA1
     INTEGER, POINTER :: i(:)
  END TYPE PRIV_INTA1
  
  TYPE PRIV_PLINEAR_VAR
     ! Type for a variable in partially linear structure
     INTEGER :: n             ! Number of points in this variable
     REAL(8), POINTER :: x(:) ! Actual variable
  END TYPE PRIV_PLINEAR_VAR

  TYPE PRIV_CHEB_VAR 
     ! Type for a variable in the chebyshev structure
     INTEGER :: o             ! Number of points in this variable, equivalent to the order since we are doing interpolation not regression
     INTEGER :: d             ! Degree of the polynomial in this variable, equivalent to ord-1
     REAL(8), POINTER :: g(:) ! the np zeros of the polynomial of degree=order
     REAL(8) :: lim(2)        ! Lower and upper limit values that the variable can take
     REAL(8), POINTER :: x(:) ! Adjusted according to lim zeros of the polynomial
     REAL(8), POINTER :: t(:,:) ! Contains the value of the polynomial of degree i at value g(k) i=0,...,deg
  END TYPE PRIV_CHEB_VAR
  
  TYPE PRIV_CHEB_COEF
     REAL(8) :: c               ! The coefficient itself I WILL NOT USE THIS BUT RATHER THE a vector for mpi reasons
     REAL(8) :: den           ! The denominator for the coefficient since it does not change (i.e. does not depend on y)
     INTEGER, POINTER :: i(:) ! Index that describes to what combination this coefficient corresponds to 
  END TYPE PRIV_CHEB_COEF
  
  TYPE CHEBSTRUC
     TYPE(PRIV_CHEB_VAR), POINTER :: v(:) ! Variable type  
     TYPE(PRIV_CHEB_COEF),POINTER :: c(:) ! Coefficient type
     INTEGER :: nv,nc,cpd,pn ! Number of variables, number of coefficients,complete polynomial degree, product of n (order1*order2*...*order_nv)
     LOGICAL :: tensor ! Tensor product approximation or complete polynomial?
     TYPE(PRIV_INTA1), POINTER :: si(:) ! Will contain the sumation indices to do the tensor product for prediction
     ! TO make my life simpler in MPI I will define an outside vector containing the coefficients in the c(:) type
     REAL(8), POINTER :: a(:)
  END TYPE CHEBSTRUC
  
  TYPE PLINEARSTRUC
     TYPE(PRIV_PLINEAR_VAR), POINTER :: v(:) ! variable type
     INTEGER :: nv,lt           ! Number of variables, length of vectorized version of matrix
     TYPE(PRIV_INTA1), POINTER :: i(:) ! Will contain the indices to the vraiables?
  END TYPE PLINEARSTRUC

  TYPE MLINEARSTRUC
     TYPE(PRIV_PLINEAR_VAR), POINTER :: v(:) ! variable type
     INTEGER :: nv,lt           ! Number of variables, length of vectorized version of matrix
     TYPE(PRIV_INTA1), POINTER :: i(:) ! Will contain the indices to the vraiables?
  END TYPE MLINEARSTRUC
     
CONTAINS
  
  SUBROUTINE PLINEAR_INIT_1(ps,n,lim,typ)
    ! typ=0 is standard grid,1 is logarithmic grid
    IMPLICIT NONE
    TYPE(PLINEARSTRUC), INTENT(INOUT) :: ps
    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(IN), OPTIONAL :: typ
    REAL(8), INTENT(IN) :: lim(2)
    REAL(8) :: rn,step,x,llim(2),lo
    INTEGER :: j
    ps%nv=1
    ps%lt=n
    ALLOCATE(ps%v(ps%nv),ps%i(n))
    ps%v(1)%n=n
    ALLOCATE(ps%v(1)%x(n))
    rn=DBLE(n)
    IF (PRESENT(typ)) THEN
       IF (typ==1) THEN
          lo=MAX(0.0d0,-lim(1)+0.1d0)
          llim(1)=LOG(lim(1)+lo)
          llim(2)=LOG(lim(2)+lo)
          step=(llim(2)-llim(1))/(rn-1.0d0)
          ALLOCATE(ps%i(1)%i(1))
          ps%v(1)%x(1)=lim(1)
          ps%i(1)%i(1)=1
          x=llim(1)
          DO j=2,n
             ALLOCATE(ps%i(j)%i(1))
             x=x+step
             ps%v(1)%x(j)=EXP(x)-lo
             ps%i(j)%i(1)=j
          END DO
          ps%v(1)%x(n)=lim(2)
       ELSE 
          step=(lim(2)-lim(1))/(rn-1.0d0)
          ALLOCATE(ps%i(1)%i(1))
          ps%v(1)%x(1)=lim(1)
          ps%i(1)%i(1)=1
          DO j=2,n
             ALLOCATE(ps%i(j)%i(1))
             ps%v(1)%x(j)=ps%v(1)%x(j-1)+step
             ps%i(j)%i(1)=j
          END DO
       END IF
    ELSE
       step=(lim(2)-lim(1))/(rn-1.0d0)
       ALLOCATE(ps%i(1)%i(1))
       ps%v(1)%x(1)=lim(1)
       ps%i(1)%i(1)=1
       DO j=2,n
          ALLOCATE(ps%i(j)%i(1))
          ps%v(1)%x(j)=ps%v(1)%x(j-1)+step
          ps%i(j)%i(1)=j
       END DO
    END IF
  END SUBROUTINE PLINEAR_INIT_1
  
  SUBROUTINE PLINEAR_INIT_3(ps,n,lim,elim,typ)
    ! typ=0 is standard grid,1 is logarithmic grid
    ! lim(4), (lower,lower_grid,upper_grid,upper)    
    IMPLICIT NONE
    TYPE(PLINEARSTRUC), INTENT(INOUT) :: ps
    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(IN), OPTIONAL :: typ
    REAL(8), INTENT(IN) :: lim(2),elim(2)
    REAL(8) :: rn,step,x,llim(2),lo
    INTEGER :: j
    ps%nv=1
    ps%lt=n
    ALLOCATE(ps%v(ps%nv),ps%i(n))
    ps%v(1)%n=n
    ALLOCATE(ps%v(1)%x(n))
    rn=DBLE(n)
    IF (PRESENT(typ)) THEN
       IF (typ==1) THEN
          lo=MAX(0.0d0,-lim(1)+0.1d0)
          llim(1)=LOG(lim(1)+lo)
          llim(2)=LOG(lim(2)+lo)
          step=(llim(2)-llim(1))/(rn-3.0d0)
          ALLOCATE(ps%i(1)%i(1),ps%i(2)%i(1))
          ps%v(1)%x(1)=elim(1)
          ps%i(1)%i(1)=1
          ps%v(1)%x(2)=lim(1)
          ps%i(2)%i(1)=2
          x=llim(1)
          DO j=3,n-1
             ALLOCATE(ps%i(j)%i(1))
             x=x+step
             ps%v(1)%x(j)=EXP(x)-lo
             ps%i(j)%i(1)=j
          END DO
          ps%v(1)%x(n)=elim(2)
          ALLOCATE(ps%i(n)%i(1))
          ps%i(n)%i(1)=n
       ELSE
          step=(lim(2)-lim(1))/(rn-3.0d0)
          ALLOCATE(ps%i(1)%i(1),ps%i(2)%i(1))
          ps%v(1)%x(1)=elim(1)
          ps%i(1)%i(1)=1
          ps%v(1)%x(2)=lim(1)
          ps%i(2)%i(1)=2
          DO j=3,n-1
             ALLOCATE(ps%i(j)%i(1))
             ps%v(1)%x(j)=ps%v(1)%x(j-1)+step
             ps%i(j)%i(1)=j
          END DO
          ps%v(1)%x(n)=elim(2)
          ALLOCATE(ps%i(n)%i(1))
          ps%i(n)%i(1)=n
       END IF
    ELSE
       step=(lim(2)-lim(1))/(rn-3.0d0)
       ALLOCATE(ps%i(1)%i(1),ps%i(2)%i(1),ps%i(n)%i(1))
       ps%v(1)%x(1)=elim(1)
       ps%i(1)%i(1)=1
       ps%v(1)%x(2)=lim(1)
       ps%i(2)%i(1)=2
       DO j=3,n-1
          ALLOCATE(ps%i(j)%i(1))
          ps%v(1)%x(j)=ps%v(1)%x(j-1)+step
          ps%i(j)%i(1)=j
       END DO
       ps%v(1)%x(n)=elim(2)
       ps%i(n)%i(1)=n
    END IF
  END SUBROUTINE PLINEAR_INIT_3

  SUBROUTINE PLINEAR_INIT_2(ps,n,lim,typ)
    IMPLICIT NONE
    TYPE(PLINEARSTRUC), INTENT(INOUT) :: ps
    INTEGER, INTENT(IN) :: n(:)
    INTEGER, INTENT(IN), OPTIONAL :: typ(:)
    REAL(8), INTENT(IN) :: lim(:,:)
    REAL(8) :: x,rn,step,llim(2),lo
    INTEGER :: m,j,v,cumprod(0:SIZE(n))
    ps%nv=SIZE(n)
    ps%lt=PRODUCT(n)
    ALLOCATE(ps%v(ps%nv),ps%i(ps%lt))
    cumprod(0)=1.0d0
    DO v=1,ps%nv
       rn=DBLE(n(v))
       ps%v(v)%n=n(v)
       ALLOCATE(ps%v(v)%x(n(v)))
       cumprod(v)=cumprod(v-1)*ps%v(v)%n
       IF (PRESENT(typ)) THEN
          IF (typ(v)==1) THEN
             lo=MAX(0.0d0,-lim(v,1)+0.1d0)
             llim(1)=LOG(lim(v,1)+lo)
             llim(2)=LOG(lim(v,2)+lo)
             step=(llim(2)-llim(1))/(rn-1.0d0)
             ps%v(v)%x(1)=lim(v,1)
             x=llim(1)
             DO j=2,ps%v(v)%n
                x=x+step
                ps%v(v)%x(j)=EXP(x)-lo
             END DO
             ps%v(v)%x(ps%v(v)%n)=lim(v,2)
          ELSE
             step=(lim(v,2)-lim(v,1))/(rn-1.0d0)
             ps%v(v)%x(1)=lim(v,1)
             DO j=2,ps%v(v)%n
                ps%v(v)%x(j)=ps%v(v)%x(j-1) + step
             END DO
          END IF
       ELSE
          step=(lim(v,2)-lim(v,1))/(rn-1.0d0)
          ps%v(v)%x(1)=lim(v,1)
          DO j=2,ps%v(v)%n
             ps%v(v)%x(j)=ps%v(v)%x(j-1) + step
          END DO
       END IF
    END DO
    DO m=1,ps%lt
       ALLOCATE(ps%i(m)%i(ps%nv))
       DO j=1,ps%nv
          v=((m-1)/cumprod(j-1)) - ((m-1)/cumprod(j))*ps%v(j)%n + 1
          ps%i(m)%i(j)=v
       END DO
    END DO
  END SUBROUTINE PLINEAR_INIT_2

  SUBROUTINE PLINEAR_INIT_4(ps,n,lim,elim,typ)
    IMPLICIT NONE
    TYPE(PLINEARSTRUC), INTENT(INOUT) :: ps
    INTEGER, INTENT(IN) :: n(:)
    INTEGER, INTENT(IN), OPTIONAL :: typ(:)
    REAL(8), INTENT(IN) :: lim(:,:),elim(:,:)
    REAL(8) :: x,rn,step,llim(2),lo
    INTEGER :: m,j,v,cumprod(0:SIZE(n))
    ps%nv=SIZE(n)
    ps%lt=PRODUCT(n)
    ALLOCATE(ps%v(ps%nv),ps%i(ps%lt))
    cumprod(0)=1.0d0
    DO v=1,ps%nv
       rn=DBLE(n(v))
       ps%v(v)%n=n(v)
       ALLOCATE(ps%v(v)%x(n(v)))
       cumprod(v)=cumprod(v-1)*ps%v(v)%n
       IF (PRESENT(typ)) THEN
          IF (typ(v)==1) THEN
             lo=MAX(0.0d0,-lim(v,1)+0.1d0)
             llim(1)=LOG(lim(v,1)+lo)
             llim(2)=LOG(lim(v,2)+lo)
             step=(llim(2)-llim(1))/(rn-3.0d0)
             ps%v(v)%x(1)=elim(v,1)
             ps%v(v)%x(2)=lim(v,1)
             x=llim(1)
             DO j=3,ps%v(v)%n-1
                x=x+step
                ps%v(v)%x(j)=EXP(x)-lo
             END DO
             ps%v(v)%x(ps%v(v)%n)=elim(v,2)
          ELSE
             step=(lim(v,2)-lim(v,1))/(rn-3.0d0)
             ps%v(v)%x(1)=elim(v,1)
             ps%v(v)%x(2)=lim(v,1)
             DO j=3,ps%v(v)%n-1
                ps%v(v)%x(j)=ps%v(v)%x(j-1) + step
             END DO
             ps%v(v)%x(ps%v(v)%n)=elim(v,2)
          END IF
       ELSE
          step=(lim(v,2)-lim(v,1))/(rn-3.0d0)
          ps%v(v)%x(1)=elim(v,1)
          ps%v(v)%x(2)=lim(v,1)
          DO j=3,ps%v(v)%n-1
             ps%v(v)%x(j)=ps%v(v)%x(j-1) + step
          END DO
          ps%v(v)%x(ps%v(v)%n)=elim(v,2)
       END IF
    END DO
    DO m=1,ps%lt
       ALLOCATE(ps%i(m)%i(ps%nv))
       DO j=1,ps%nv
          v=((m-1)/cumprod(j-1)) - ((m-1)/cumprod(j))*ps%v(j)%n + 1
          ps%i(m)%i(j)=v
       END DO
    END DO
  END SUBROUTINE PLINEAR_INIT_4

  SUBROUTINE MLINEAR_INIT_1(ps,n,lim,typ)
    IMPLICIT NONE
    TYPE(MLINEARSTRUC), INTENT(INOUT) :: ps
    INTEGER, INTENT(IN) :: n(:)
    INTEGER, INTENT(IN), OPTIONAL :: typ(:)
    REAL(8), INTENT(IN) :: lim(:,:)
    REAL(8) :: x,rn,step,llim(2),lo
    INTEGER :: m,j,v,cumprod(0:SIZE(n))
    ps%nv=SIZE(n)
    ps%lt=PRODUCT(n)
    ALLOCATE(ps%v(ps%nv),ps%i(ps%lt))
    cumprod(0)=1.0d0
    DO v=1,ps%nv
       rn=DBLE(n(v))
       ps%v(v)%n=n(v)
       ALLOCATE(ps%v(v)%x(n(v)))
       cumprod(v)=cumprod(v-1)*ps%v(v)%n
       IF (present(typ)) THEN
          IF (typ(v)==1) THEN
             lo=MAX(0.0d0,-lim(v,1)+0.1d0)
             llim(1)=LOG(lim(v,1)+lo)
             llim(2)=LOG(lim(v,2)+lo)
             step=(llim(2)-llim(1))/(rn-1.0d0)
             ps%v(v)%x(1)=lim(v,1)
             x=llim(1)
             DO j=2,ps%v(v)%n
                x=x+step
                ps%v(v)%x(j)=EXP(x)-lo
             END DO
             ps%v(v)%x(ps%v(v)%n)=lim(v,2)
          ELSE
             step=(lim(v,2)-lim(v,1))/(rn-1.0d0)
             ps%v(v)%x(1)=lim(v,1)
             DO j=2,ps%v(v)%n
                ps%v(v)%x(j)=ps%v(v)%x(j-1) + step
             END DO
          END IF
       ELSE
          step=(lim(v,2)-lim(v,1))/(rn-1.0d0)
          ps%v(v)%x(1)=lim(v,1)
          DO j=2,ps%v(v)%n
             ps%v(v)%x(j)=ps%v(v)%x(j-1) + step
          END DO
       END IF
    END DO
    DO m=1,ps%lt
       ALLOCATE(ps%i(m)%i(ps%nv))
       DO j=1,ps%nv
          v=((m-1)/cumprod(j-1)) - ((m-1)/cumprod(j))*ps%v(j)%n + 1
          ps%i(m)%i(j)=v
       END DO
    END DO
  END SUBROUTINE MLINEAR_INIT_1

  SUBROUTINE MLINEAR_INIT_2(ps,n,lim,elim,typ)
    IMPLICIT NONE
    TYPE(MLINEARSTRUC), INTENT(INOUT) :: ps
    INTEGER, INTENT(IN) :: n(:)
    INTEGER, INTENT(IN), OPTIONAL :: typ(:)
    REAL(8), INTENT(IN) :: lim(:,:),elim(:,:)
    REAL(8) :: x,rn,step,llim(2),lo
    INTEGER :: m,j,v,cumprod(0:SIZE(n))
    ps%nv=SIZE(n)
    ps%lt=PRODUCT(n)
    ALLOCATE(ps%v(ps%nv),ps%i(ps%lt))
    cumprod(0)=1.0d0
    DO v=1,ps%nv
       rn=DBLE(n(v))
       ps%v(v)%n=n(v)
       ALLOCATE(ps%v(v)%x(n(v)))
       cumprod(v)=cumprod(v-1)*ps%v(v)%n
       IF (PRESENT(typ)) THEN
          IF (typ(v)==1) THEN
             lo=MAX(0.0d0,-lim(v,1)+0.1d0)
             llim(1)=LOG(lim(v,1)+lo)
             llim(2)=LOG(lim(v,2)+lo)
             step=(llim(2)-llim(1))/(rn-3.0d0)
             ps%v(v)%x(1)=elim(v,1)
             ps%v(v)%x(2)=lim(v,1)
             x=llim(1)
             DO j=3,ps%v(v)%n-1
                x=x+step
                ps%v(v)%x(j)=EXP(x)-lo
             END DO
             ps%v(v)%x(ps%v(v)%n)=elim(v,2)
          ELSE
             step=(lim(v,2)-lim(v,1))/(rn-3.0d0)
             ps%v(v)%x(1)=elim(v,1)
             ps%v(v)%x(2)=lim(v,1)
             DO j=3,ps%v(v)%n-1
                ps%v(v)%x(j)=ps%v(v)%x(j-1) + step
             END DO
             ps%v(v)%x(ps%v(v)%n)=elim(v,2)
          END IF
       ELSE
          step=(lim(v,2)-lim(v,1))/(rn-3.0d0)
          ps%v(v)%x(1)=elim(v,1)
          ps%v(v)%x(2)=lim(v,1)
          DO j=3,ps%v(v)%n-1
             ps%v(v)%x(j)=ps%v(v)%x(j-1) + step
          END DO
          ps%v(v)%x(ps%v(v)%n)=elim(v,2)
       END IF
    END DO
    DO m=1,ps%lt
       ALLOCATE(ps%i(m)%i(ps%nv))
       DO j=1,ps%nv
          v=((m-1)/cumprod(j-1)) - ((m-1)/cumprod(j))*ps%v(j)%n + 1
          ps%i(m)%i(j)=v
       END DO
    END DO
  END SUBROUTINE MLINEAR_INIT_2

  REAL(8) FUNCTION MLINEAR_INTER(ps,x,y)
    IMPLICIT NONE
    TYPE(MLINEARSTRUC), INTENT(IN) :: ps
    REAL(8), INTENT(IN) :: x(:),y(:)
    INTEGER :: coord(ps%nv),v,cumprod(0:ps%nv),ind
    INTEGER :: new_coord,n,j,coun(ps%nv),pow(ps%nv),add(ps%nv)
    REAL(8) :: scal,den,coef,f
    ! Let's first find the coordinates of the hypercube where my point is located
    ! and put the values of the function that you need in the order needed
    ! To do that first we need to transform the j index into the appropriate "i" coordinates
    ! The formula is ind=cum_prod(0)*(coord(1)-1) + cum_prod(1)*(coord(2)-1) + ... + cum_prod(d-1)*(coord(d)-1) + 1
    n=2**ps%nv
    f=0.0d0

    ind=1
    cumprod(0)=1
    coef=1.0d0
    DO v=1,ps%nv
       coord(v)=LOCATE(ps%v(v)%x,x(v))
       ind=ind + cumprod(v-1)*(coord(v)-1)
       cumprod(v)=cumprod(v-1)*ps%v(v)%n
       den=ps%v(v)%x(coord(v)+1)-ps%v(v)%x(coord(v))       
       scal=1.0d0 - (ABS(x(v)-ps%v(v)%x(coord(v)))/den)
       coef=coef*scal
       pow(v)=2**(v-1)
    END DO
    f=f+y(ind)*coef

    coun=0
    add=0
    DO j=2,n
       ind=1
       cumprod(0)=1
       coef=1.0d0
       DO v=1,ps%nv
          coun(v)=coun(v)+1
          IF (coun(v)==pow(v)) THEN
             coun(v)=0
             IF (add(v)==0) THEN
                add(v)=1
             ELSE IF (add(v)==1) THEN
                add(v)=0
             END IF
          END IF
          new_coord=coord(v)+add(v)
          ind=ind + cumprod(v-1)*(new_coord-1)
          cumprod(v)=cumprod(v-1)*ps%v(v)%n
          den=ps%v(v)%x(coord(v)+1)-ps%v(v)%x(coord(v))       
          scal=1.0d0 - (ABS(x(v)-ps%v(v)%x(new_coord))/den)
          coef=coef*scal
       END DO
       f=f+y(ind)*coef
    END DO
    MLINEAR_INTER=f
  END FUNCTION MLINEAR_INTER

  REAL(8) FUNCTION PLINEAR_INTER_1(ps,x,y)
    IMPLICIT NONE 
    TYPE(PLINEARSTRUC), INTENT(IN) :: ps
    REAL(8), INTENT(IN) :: x,y(:)
    INTEGER :: coord
    REAL(8) :: den,scaled_x,f(0:1)
    coord=LOCATE(ps%v(1)%x,x)
    den=ps%v(1)%x(coord+1) - ps%v(1)%x(coord)
    scaled_x=(x-ps%v(1)%x(coord))/den
    f(0)=y(coord)
    PLINEAR_INTER_1=f(0)
    f(1)=y(coord+1)
    PLINEAR_INTER_1=PLINEAR_INTER_1 + scaled_x*(f(1)-f(0))
  END FUNCTION PLINEAR_INTER_1

  REAL(8) FUNCTION PLINEAR_INTER_2(ps,x,y)
    IMPLICIT NONE
    TYPE(PLINEARSTRUC), INTENT(IN) :: ps
    REAL(8), INTENT(IN) :: x(:),y(:)
    INTEGER :: coord(ps%nv),j,v,sortind(ps%nv),cumprod(0:ps%nv),ind,k
    INTEGER :: cumprod_p1(0:ps%nv),ind_p1(ps%nv),h,hh
    REAL(8) :: scaled_x(ps%nv),f(0:ps%nv),den
    ! Let's first find the coordinates of the hypercube where my point is located
    ! and put the values of the function that you need in the order needed
    ! To do that first we need to transform the j index into the appropriate "i" coordinates
    ! The formula is ind=cum_prod(0)*(coord(1)-1) + cum_prod(1)*(coord(2)-1) + ... + cum_prod(d-1)*(coord(d)-1) + 1
    ind=1
    cumprod(0)=1
     DO v=1,ps%nv
       coord(v)=LOCATE(ps%v(v)%x,x(v))
       den=ps%v(v)%x(coord(v)+1)-ps%v(v)%x(coord(v))
       scaled_x(v)=(x(v)-ps%v(v)%x(coord(v)))/den
       sortind(v)=v
       ind=ind + cumprod(v-1)*(coord(v)-1)
       cumprod(v)=cumprod(v-1)*ps%v(v)%n
    END DO
    ! sortind gives me the sorting index of the scaled point y in descending form
    CALL P_Sort_Interp(-scaled_x,f(1:ps%nv),sortind)
    DO j=1,ps%nv
       ! ind_p1(k) gives me the index when I move in the kth direction
       ! _p1 means plus 1 (in the coordinate sense of the unit hypercube)
       k=sortind(j)
       cumprod_p1(0)=1
       ind_p1(k)=1 
       coord(k)=coord(k)+1
       DO h=1,ps%nv
          ind_p1(k)=ind_p1(k) + cumprod_p1(h-1)*(coord(h)-1)
          cumprod_p1(h)=cumprod_p1(h-1)*ps%v(h)%n
       END DO
    END DO
    f(0)=y(ind)
    PLINEAR_INTER_2=f(0)
    DO j=1,ps%nv
       k=sortind(j)
       f(j)=y(ind_p1(k))
       PLINEAR_INTER_2=PLINEAR_INTER_2 + scaled_x(k)*(f(j)-f(j-1))
    END DO
  END FUNCTION PLINEAR_INTER_2

  SUBROUTINE INIT_CHEB_1(cs,n,limit)
    IMPLICIT NONE
    TYPE(CHEBSTRUC), INTENT(INOUT) :: cs
    REAL(8), INTENT(IN) :: limit(2)
    INTEGER, INTENT(IN) :: n
    INTEGER :: i,j,k
    REAL(8) :: rn
    cs%nv=1 ! Number of variables
    cs%pn=n ! cumulative product of orders of poynomials (trivial in this case)
    cs%tensor=.TRUE. ! Trivial definition of this case as a tensor product
    cs%nc=n ! trivial number of coefficients n
    ALLOCATE(cs%v(cs%nv),cs%c(cs%nc),cs%si(cs%pn),cs%a(cs%nc))
    rn=DBLE(n)
    DO i=1,cs%nv ! Trivial since there is only one variable
       cs%v(i)%o=n
       cs%v(i)%d=n-1
       cs%v(i)%lim=limit
       ALLOCATE(cs%v(i)%g(n),cs%v(i)%x(n),cs%v(i)%t(0:n-1,n))
       DO j=1,cs%v(i)%o
          cs%v(i)%g(j)=-COS( ( (2.0d0*DBLE(j))-1.0d0 ) * halfpi/rn )
          cs%v(i)%x(j)=( (cs%v(i)%g(j)+1.0d0)*0.5d0*(cs%v(i)%lim(2)-cs%v(i)%lim(1)) ) + cs%v(i)%lim(1)
          cs%v(i)%t(0,j)=1.0d0
          cs%v(i)%t(1,j)=cs%v(i)%g(j)
          DO k=2,cs%v(i)%d
             cs%v(i)%t(k,j) = ( 2.0d0*cs%v(i)%g(j)*cs%v(i)%t(k-1,j) ) - cs%v(i)%t(k-2,j)
          END DO
       END DO
    END DO
    ! In this case nc and pn are the same and so are the indices so I'll do it in one step
    DO j=1,cs%nc
       ALLOCATE(cs%c(j)%i(cs%nv))
       cs%c(j)%i(1)=j-1
       ! There is just one variable
       rn=0.0d0
       DO k=1,cs%v(1)%o
          rn=rn+cs%v(1)%t(cs%c(j)%i(1),k)**2
       END DO
       cs%c(j)%den=rn       
    END DO
    DO i=1,cs%pn
       ALLOCATE(cs%si(i)%i(cs%nv))
       cs%si(i)%i=i
    END DO
  END SUBROUTINE INIT_CHEB_1
  
  SUBROUTINE INIT_CHEB_2(cs,n,limit,tensor)
    IMPLICIT NONE
    TYPE(CHEBSTRUC), INTENT(INOUT) :: cs
    REAL(8), INTENT(IN) :: limit(:,:)
    INTEGER, INTENT(IN) :: n(:)
    LOGICAL, INTENT(IN) :: tensor
    INTEGER :: i,j,k,m,cumprod(0:SIZE(n)),coord(PRODUCT(n),SIZE(n))
    REAL(8) :: rn
    cs%nv=SIZE(n) ! Number of variables
    cs%pn=PRODUCT(n) ! cumulative product of orders of poynomials
    ALLOCATE(cs%v(cs%nv),cs%si(cs%pn))
    cumprod(0)=1.0d0
    cs%cpd=0
    DO i=1,cs%nv 
       rn=DBLE(n(i))
       cs%v(i)%o=n(i)
       cs%v(i)%d=n(i)-1
       cs%v(i)%lim=limit(i,:)
       ALLOCATE(cs%v(i)%g(n(i)),cs%v(i)%x(n(i)),cs%v(i)%t(0:n(i)-1,n(i)))
       DO j=1,cs%v(i)%o
          cs%v(i)%g(j)=-COS( ( (2.0d0*DBLE(j))-1.0d0 ) * halfpi/rn )
          cs%v(i)%x(j)=( (cs%v(i)%g(j)+1.0d0)*0.5d0*(cs%v(i)%lim(2)-cs%v(i)%lim(1)) ) + cs%v(i)%lim(1)
          cs%v(i)%t(0,j)=1.0d0
          cs%v(i)%t(1,j)=cs%v(i)%g(j)
          DO k=2,cs%v(i)%d
             cs%v(i)%t(k,j) = ( 2.0d0*cs%v(i)%g(j)*cs%v(i)%t(k-1,j) ) - cs%v(i)%t(k-2,j)
          END DO
       END DO
       cumprod(i)=cumprod(i-1)*cs%v(i)%o
       IF (cs%cpd<cs%v(i)%o) cs%cpd=cs%v(i)%o ! In case we do a complete polynomial basis instead of tensor product
    END DO
    cs%tensor=tensor ! Tensor product or complete polynomial?
    IF (cs%tensor) THEN
       cs%nc=cs%pn
       ALLOCATE(cs%c(cs%nc),cs%a(cs%nc))
       DO m=1,cs%nc
          ! We are already here so we might as well take advantage and do the outer index
          ALLOCATE(cs%c(m)%i(cs%nv),cs%si(m)%i(cs%nv))
          cs%c(m)%den=1.0d0
          DO j=1,cs%nv
             k=((m-1)/cumprod(j-1)) - ((m-1)/cumprod(j))*cs%v(j)%o
             cs%c(m)%i(j)=k
             cs%si(m)%i(j)=k+1
             rn=0.0d0
             DO k=1,cs%v(j)%o
                rn=rn+cs%v(j)%t(cs%c(m)%i(j),k)**2
             END DO
             cs%c(m)%den=cs%c(m)%den*rn
          END DO
       END DO
    ELSE
       i=0
       DO m=1,cs%pn
          ! We are already doing this so we might as well do it here for the outer index
          ALLOCATE(cs%si(m)%i(cs%nv))
          k=0
          DO j=1,cs%nv
             cs%si(m)%i(j)=((m-1)/cumprod(j-1)) - ((m-1)/cumprod(j))*cs%v(j)%o
             k=k+cs%si(m)%i(j)
          END DO
          IF (k<=cs%cpd) THEN
             i=i+1
             coord(i,:)=cs%si(m)%i
          END IF
          cs%si(m)%i=cs%si(m)%i+1
       END DO
       cs%nc=i ! This is the number of coefficients that satisy the cpd restriction
       ALLOCATE(cs%c(cs%nc),cs%a(cs%nc))
       DO m=1,cs%nc
          ALLOCATE(cs%c(m)%i(cs%nv))
          cs%c(m)%i=coord(m,:)
          cs%c(m)%den=1.0d0
          DO j=1,cs%nv
             rn=0.0d0
             DO k=1,cs%v(j)%o
                rn=rn+cs%v(j)%t(cs%c(m)%i(j),k)**2
             END DO
             cs%c(m)%den=cs%c(m)%den*rn
          END DO
       END DO
    END IF
  END SUBROUTINE INIT_CHEB_2
  
  SUBROUTINE CHEB_COEF(cs,y)
    IMPLICIT NONE
    TYPE(CHEBSTRUC), INTENT(INOUT) :: cs
    REAL(8), INTENT(IN) :: y(:)
    INTEGER :: i,k,v,j
    REAL(8) :: num,n1
    DO i=1,cs%nc
       num=0.0d0
       DO k=1,cs%pn
          n1=y(k)
          DO v=1,cs%nv
             n1=n1*cs%v(v)%t(cs%c(i)%i(v),cs%si(k)%i(v))
          END DO
          num=num+n1
       END DO
       cs%c(i)%c=num/cs%c(i)%den
       cs%a(i)=cs%c(i)%c
    END DO
  END SUBROUTINE CHEB_COEF
  
  REAL(8) FUNCTION CHEB_INTER_1(cs,x)
    IMPLICIT NONE
    TYPE(CHEBSTRUC), INTENT(IN) :: cs
    REAL(8), INTENT(IN) :: x
    REAL(8) :: out,xadj,t(0:cs%v(1)%d)
    INTEGER :: k
    xadj=(2.0d0*(x-cs%v(1)%lim(1))/(cs%v(1)%lim(2)-cs%v(1)%lim(1))) - 1.0d0
    t(0)=1.0d0
    t(1)=xadj
    DO k=2,cs%v(1)%d
       t(k)= ( 2.0d0*xadj*t(k-1) ) - t(k-2)
    END DO
    out=0.0d0
    DO k=1,cs%nc
       out=out+cs%a(k)*t(cs%c(k)%i(1))
    END DO
    CHEB_INTER_1=out
  END FUNCTION CHEB_INTER_1
  
  REAL(8) FUNCTION CHEB_INTER_2(cs,x)
    IMPLICIT NONE
    TYPE(CHEBSTRUC), INTENT(IN) :: cs
    REAL(8), INTENT(IN) :: x(:)
    REAL(8) :: xadj,out,aux
    TYPE(ARRAY1D) :: t(cs%nv)
    INTEGER :: j,m,k
    DO j=1,cs%nv
       xadj=(2.0d0*(x(j)-cs%v(j)%lim(1))/(cs%v(j)%lim(2)-cs%v(j)%lim(1))) - 1.0d0
       ALLOCATE(t(j)%a(0:cs%v(j)%d))
       t(j)%a(0)=1.0d0
       t(j)%a(1)=xadj
       DO k=2,cs%v(j)%d
          t(j)%a(k)= ( 2.0d0*xadj*t(j)%a(k-1) ) - t(j)%a(k-2)
       END DO
    END DO
    out=0.0d0
    DO m=1,cs%nc
       aux=cs%a(m)
       DO j=1,cs%nv
          aux=aux*t(j)%a(cs%c(m)%i(j))
       END DO
       out=out+aux
    END DO
    DO j=1,cs%nv
       DEALLOCATE(t(j)%a)
    END DO
    CHEB_INTER_2=out
  END FUNCTION CHEB_INTER_2
  
  FUNCTION LOCATE(xx,x)
    !Given an array xx(1:n), and given a value x, returns a value j such that x is between
    !xx(j) and xx(j+1). xx(1:n) must be monotonic, either increasing or decreasing. j=0
    !or j=n is returned to indicate that x is out of range.
    IMPLICIT NONE
    REAL(8), DIMENSION(:), INTENT(IN) :: xx
    REAL(8), INTENT(IN) :: x
    INTEGER :: locate
    INTEGER :: n,jl,jm,ju
    LOGICAL :: ascnd
    n=SIZE(xx)
    IF (ABS(xx(n)-x)<=2.220446049250313E-015) THEN
       locate=n-1
       RETURN
    END IF
    ascnd = (xx(n) >= xx(1))
    IF (ascnd.EQV..FALSE.) THEN
       WRITE(*,*) "Locate: data not in ascending order"
       STOP
    END IF
    jl=0
    ju=n+1
    DO
       IF (ju-jl <= 1) EXIT
       jm=(ju+jl)/2
       IF (ascnd .EQV. (x >= xx(jm))) THEN
          jl=jm
       ELSE
          ju=jm
       END IF
    END DO
    locate=jl
  END FUNCTION LOCATE
  
  
!!!!!!!!!!!!!!!!!!!!!! HELPER FUNCTIONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE P_Sort_Interp(Xin,xout,ind)
    !  Sorts xout into ascENDing order - _quicksort_stat
    REAL (8), DIMENSION (:), INTENT (in) :: Xin
    INTEGER, INTENT(INOUT) :: ind(:)
    REAL (8), DIMENSION (SIZE(ind)), INTENT (Out) :: xout
    xout=XIN
    Call P_SUBSOR_Interp(xout,1,Size(xout),ind)
    Call P_Inssort_Interp(xout,ind)
  CONTAINS
    RECURSIVE SUBROUTINE P_SUBSOR_Interp (xout, IDEB1, IFIN1,ind)
      !  Sorts xout from IDEB1 to IFIN1
      REAL(8), DIMENSION (:), INTENT (InOut) :: xout
      INTEGER, INTENT (In) :: IDEB1, IFIN1
      INTEGER, INTENT(INOUT) :: ind(:)
      INTEGER, PARAMETER :: NINS = 16 ! Max for insertion sort
      INTEGER :: ICRS, IDEB, IDCR, IFIN, IMIL,IWRK,IPIV
      REAL(8) :: XPIV, XWRK
      IDEB = IDEB1
      IFIN = IFIN1
      !  IF we DOn't have enough values to make it worth while, we leave
      !  them unsorted, and the final insertion sort will take care of them
      IF ((IFIN - IDEB) > NINS) THEN
         IMIL = (IDEB+IFIN) / 2
         !  One chooses a pivot, median of 1st, last, and middle values
         IF (xout(IMIL) < xout(IDEB)) THEN
            XWRK = xout (IDEB)
            IWRK = IND (IDEB)
            xout (IDEB) = xout (IMIL)
            IND (IDEB) = IND (IMIL)
            xout (IMIL) = XWRK
            IND (IMIL) = IWRK
         END IF
         IF (xout(IMIL) > xout(IFIN)) Then
            XWRK = xout (IFIN)
            IWRK = IND(IFIN)
            xout (IFIN) = xout (IMIL)
            IND (IFIN) = IND (IMIL)
            xout (IMIL) = XWRK
            IND (IMIL) = IWRK
            IF (xout(IMIL) < xout(IDEB)) Then
               XWRK = xout (IDEB)
               IWRK = IND (IDEB)
               xout (IDEB) = xout (IMIL)
               IND (IDEB) = IND (IMIL)
               xout (IMIL) = XWRK
               IND (IMIL) = IWRK
            END IF
         END IF
         XPIV = xout (IMIL)
         IPIV = IND (IMIL)
         !  One exchanges values to put those > pivot in the END and
         !  those <= pivot at the beginning
         ICRS = IDEB
         IDCR = IFIN
         ECH2: DO
            DO
               ICRS = ICRS + 1
               IF (ICRS >= IDCR) Then
                  !  the first  >  pivot is IDCR
                  !  the last   <= pivot is ICRS-1
                  !  Note: IF one arrives here on the first iteration, then
                  !  the pivot is the maximum of the set, the last value is equal
                  !  to it, and one can reduce by one the size of the set to process,
                  !  as IF xout (IFIN) > XPIV
                  Exit ECH2
               END IF
               IF (xout(ICRS) > XPIV) Exit
            END DO
            DO
               IF (xout(IDCR) <= XPIV) Exit
               IDCR = IDCR - 1
               IF (ICRS >= IDCR) Then
                  !  The last value < pivot is always ICRS-1
                  Exit ECH2
               END IF
            END DO
            XWRK = xout (IDCR)
            IWRK = IND (IDCR)
            xout (IDCR) = xout (ICRS)
            IND (IDCR) = IND (ICRS)
            xout (ICRS) = XWRK
            IND (ICRS) = IWRK
         END DO ECH2
         !  One now sorts each of the two sub-intervals
         Call P_SUBSOR_Interp (xout, IDEB1, ICRS-1,IND)
         Call P_SUBSOR_Interp (xout, IDCR, IFIN1,IND)
      END IF
      RETURN
    END SUBROUTINE P_SUBSOR_INTERP
  END SUBROUTINE P_Sort_Interp

  SUBROUTINE P_Inssort_Interp (xout,IND)
    !  Sorts xout into increasing order (Insertion sort)
    REAL(8), DIMENSION (:), INTENT (InOut) :: xout
    INTEGER, INTENT(INOUT) :: ind(:)
    INTEGER :: ICRS, IDCR,IWRK
    REAL(8) :: XWRK
    DO ICRS = 2, Size (xout)
       XWRK = xout (ICRS)
       IWRK = IND (ICRS)
       IF (XWRK >= xout(ICRS-1)) Cycle
       xout (ICRS) = xout (ICRS-1)
       IND (ICRS) = IND (ICRS-1)
       DO IDCR = ICRS - 2, 1, - 1
          IF (XWRK >= xout(IDCR)) Exit
          xout (IDCR+1) = xout (IDCR)
          IND (IDCR+1) = IND (IDCR)
       END DO
       xout (IDCR+1) = XWRK
       IND (IDCR+1) = IWRK
    END DO
    RETURN
  END SUBROUTINE P_Inssort_Interp

END MODULE INTERPOL
