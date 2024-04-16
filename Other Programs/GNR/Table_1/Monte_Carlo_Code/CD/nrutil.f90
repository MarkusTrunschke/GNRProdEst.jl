MODULE NRUTIL
  USE BLAS95, ONLY : DOT,GEMM,GEMV
  IMPLICIT NONE
  REAL(8), PARAMETER :: NR_PI = 3.1415926535897932384626433832d0
  REAL(8), PARAMETER :: NR_2PI =  6.2831853071795864769252867664d0
  REAL(8), PARAMETER :: NR_SQRT2PI =  2.506628274631000502415765284779d0
  REAL(8), PARAMETER :: NR_EPS = 2.220446049250313E-016
  REAL(8), PARAMETER :: NR_BIG = 1.0d300
  REAL(8), PARAMETER :: NR_SMALL = 1.0d-300

  INTEGER, PARAMETER :: NPAR_ARTH=16,NPAR2_ARTH=8
  INTEGER, PARAMETER :: NPAR_GEOP=4,NPAR2_GEOP=2
  INTEGER, PARAMETER :: NPAR_CUMSUM=16
  INTEGER, PARAMETER :: NPAR_CUMPROD=8
  INTEGER, PARAMETER :: NPAR_POLY=8
  INTEGER, PARAMETER :: NPAR_POLYTERM=8

  INTERFACE MMUL
     MODULE PROCEDURE MMUL1,MMUL2,MMUL3,MMUL4
  END INTERFACE MMUL

  INTERFACE array_copy
     MODULE PROCEDURE array_copy_r, array_copy_i
  END INTERFACE
  INTERFACE swap
     MODULE PROCEDURE swap_i,swap_r,swap_rv,masked_swap_rs,masked_swap_rv,masked_swap_rm
  END INTERFACE
  INTERFACE reallocate
     MODULE PROCEDURE reallocate_rv,reallocate_rm,reallocate_iv,reallocate_im,reallocate_hv
  END INTERFACE
  INTERFACE imaxloc
     MODULE PROCEDURE imaxloc_r,imaxloc_i
  END INTERFACE
  INTERFACE assert
     MODULE PROCEDURE assert1,assert2,assert3,assert4,assert_v
  END INTERFACE
  INTERFACE assert_eq
     MODULE PROCEDURE assert_eq2,assert_eq3,assert_eq4,assert_eqn
  END INTERFACE
  INTERFACE arth
     MODULE PROCEDURE arth_r,arth_i
  END INTERFACE
  INTERFACE geop
     MODULE PROCEDURE geop_r,geop_i,geop_rv
  END INTERFACE
  INTERFACE cumsum
     MODULE PROCEDURE cumsum_r,cumsum_i
  END INTERFACE
  INTERFACE poly
     MODULE PROCEDURE poly_rr,poly_rrv,poly_msk_rrv
  END INTERFACE
  INTERFACE outerdiff
     MODULE PROCEDURE outerdiff_r,outerdiff_i
  END INTERFACE

! My added types for use in some routines, this is important
  TYPE ARRAY1D
     REAL(8), POINTER :: a(:)
     INTEGER :: d1
  END TYPE ARRAY1D

  TYPE IARRAY1D
     INTEGER, POINTER :: a(:)
     INTEGER :: d1
  END TYPE IARRAY1D

  TYPE ARRAY2D
     REAL(8), POINTER :: a(:,:)
     INTEGER :: d1,d2
  END TYPE ARRAY2D
  
  TYPE ARRAY3D
     REAL(8), POINTER :: a(:,:,:)
     INTEGER :: d1,d2,d3
  END TYPE ARRAY3D

CONTAINS

  REAL(8) FUNCTION DOTP(a,b)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: a(:),b(:)
    INTEGER :: u
    u=ASSERT_EQ(SIZE(A),SIZE(B),"dotp")
    DOTP=dot(a,b)
  END FUNCTION DOTP

  FUNCTION MMUL1(A,B)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: A(:,:),B(:,:)
    REAL(8) :: MMUL1(SIZE(A,1),SIZE(B,2))
    INTEGER :: u
    u=ASSERT_EQ(SIZE(A,2),SIZE(B,1),"mmul1")
    CALL GEMM(A,B,MMUL1)    
  END FUNCTION MMUL1

  FUNCTION MMUL2(A,B)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: A(:,:),B(:)
    REAL(8) :: MMUL2(SIZE(A,1))
    INTEGER :: u
    u=ASSERT_EQ(SIZE(A,2),SIZE(B),"mmul2")
    CALL GEMV(A,B,MMUL2)    
  END FUNCTION MMUL2

  FUNCTION MMUL3(A,B)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: A(:),B(:,:)
    REAL(8) :: MMUL3(SIZE(B,2))
    INTEGER :: u
    u=ASSERT_EQ(SIZE(A),SIZE(B,1),"mmul3")
    CALL GEMV(TRANSPOSE(B),A,MMUL3)    
  END FUNCTION MMUL3

  FUNCTION MMUL4(A,B)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: A(:),B(:)
    REAL(8) :: MMUL4
    INTEGER :: u
    u=ASSERT_EQ(SIZE(A),SIZE(B),"mmul4")
    MMUL4=DOTP(A,B)
  END FUNCTION MMUL4

  SUBROUTINE array_copy_r(src,dest,n_copied,n_not_copied)
    REAL(8), DIMENSION(:), INTENT(IN) :: src
    REAL(8), DIMENSION(:), INTENT(OUT) :: dest
    INTEGER, INTENT(OUT) :: n_copied, n_not_copied
    n_copied=min(size(src),size(dest))
    n_not_copied=size(src)-n_copied
    dest(1:n_copied)=src(1:n_copied)
  END SUBROUTINE array_copy_r

  SUBROUTINE array_copy_i(src,dest,n_copied,n_not_copied)
    INTEGER, DIMENSION(:), INTENT(IN) :: src
    INTEGER, DIMENSION(:), INTENT(OUT) :: dest
    INTEGER, INTENT(OUT) :: n_copied, n_not_copied
    n_copied=min(size(src),size(dest))
    n_not_copied=size(src)-n_copied
    dest(1:n_copied)=src(1:n_copied)
  END SUBROUTINE array_copy_i

  SUBROUTINE swap_i(a,b)
    INTEGER, INTENT(INOUT) :: a,b
    INTEGER :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_i

  SUBROUTINE swap_r(a,b)
    REAL(8), INTENT(INOUT) :: a,b
    REAL(8) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_r

  SUBROUTINE swap_rv(a,b)
    REAL(8), DIMENSION(:), INTENT(INOUT) :: a,b
    REAL(8), DIMENSION(SIZE(a)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_rv

  SUBROUTINE masked_swap_rs(a,b,mask)
    REAL(8), INTENT(INOUT) :: a,b
    LOGICAL, INTENT(IN) :: mask
    REAL(8) :: swp
    if (mask) then
       swp=a
       a=b
       b=swp
    end if
  END SUBROUTINE masked_swap_rs

  SUBROUTINE masked_swap_rv(a,b,mask)
    REAL(8), DIMENSION(:), INTENT(INOUT) :: a,b
    LOGICAL, DIMENSION(:), INTENT(IN) :: mask
    REAL(8), DIMENSION(size(a)) :: swp
    where (mask)
       swp=a
       a=b
       b=swp
    end where
  END SUBROUTINE masked_swap_rv

  SUBROUTINE masked_swap_rm(a,b,mask)
    REAL(8), DIMENSION(:,:), INTENT(INOUT) :: a,b
    LOGICAL, DIMENSION(:,:), INTENT(IN) :: mask
    REAL(8), DIMENSION(size(a,1),size(a,2)) :: swp
    where (mask)
       swp=a
       a=b
       b=swp
    end where
  END SUBROUTINE masked_swap_rm

  FUNCTION reallocate_rv(p,n)
    REAL(8), DIMENSION(:), POINTER :: p, reallocate_rv
    INTEGER, INTENT(IN) :: n
    INTEGER :: nold,ierr
    allocate(reallocate_rv(n),stat=ierr)
    if (ierr /= 0) call &
         nrerror('reallocate_rv: problem in attempt to allocate memory')
    if (.not. associated(p)) RETURN
    nold=size(p)
    reallocate_rv(1:min(nold,n))=p(1:min(nold,n))
    deallocate(p)
  END FUNCTION reallocate_rv

  FUNCTION reallocate_iv(p,n)
    INTEGER, DIMENSION(:), POINTER :: p, reallocate_iv
    INTEGER, INTENT(IN) :: n
    INTEGER :: nold,ierr
    allocate(reallocate_iv(n),stat=ierr)
    if (ierr /= 0) call &
         nrerror('reallocate_iv: problem in attempt to allocate memory')
    if (.not. associated(p)) RETURN
    nold=size(p)
    reallocate_iv(1:min(nold,n))=p(1:min(nold,n))
    deallocate(p)
  END FUNCTION reallocate_iv

  FUNCTION reallocate_hv(p,n)
    CHARACTER(1), DIMENSION(:), POINTER :: p, reallocate_hv
    INTEGER, INTENT(IN) :: n
    INTEGER :: nold,ierr
    allocate(reallocate_hv(n),stat=ierr)
    if (ierr /= 0) call &
         nrerror('reallocate_hv: problem in attempt to allocate memory')
    if (.not. associated(p)) RETURN
    nold=size(p)
    reallocate_hv(1:min(nold,n))=p(1:min(nold,n))
    deallocate(p)
  END FUNCTION reallocate_hv

  FUNCTION reallocate_rm(p,n,m)
    REAL(8), DIMENSION(:,:), POINTER :: p, reallocate_rm
    INTEGER, INTENT(IN) :: n,m
    INTEGER :: nold,mold,ierr
    allocate(reallocate_rm(n,m),stat=ierr)
    if (ierr /= 0) call &
         nrerror('reallocate_rm: problem in attempt to allocate memory')
    if (.not. associated(p)) RETURN
    nold=size(p,1)
    mold=size(p,2)
    reallocate_rm(1:min(nold,n),1:min(mold,m))=&
         p(1:min(nold,n),1:min(mold,m))
    deallocate(p)
  END FUNCTION reallocate_rm

  FUNCTION reallocate_im(p,n,m)
    INTEGER, DIMENSION(:,:), POINTER :: p, reallocate_im
    INTEGER, INTENT(IN) :: n,m
    INTEGER :: nold,mold,ierr
    allocate(reallocate_im(n,m),stat=ierr)
    if (ierr /= 0) call &
         nrerror('reallocate_im: problem in attempt to allocate memory')
    if (.not. associated(p)) RETURN
    nold=size(p,1)
    mold=size(p,2)
    reallocate_im(1:min(nold,n),1:min(mold,m))=&
         p(1:min(nold,n),1:min(mold,m))
    deallocate(p)
  END FUNCTION reallocate_im

  FUNCTION ifirstloc(mask)
    LOGICAL, DIMENSION(:), INTENT(IN) :: mask
    INTEGER :: ifirstloc
    INTEGER, DIMENSION(1) :: loc
    loc=maxloc(merge(1,0,mask))
    ifirstloc=loc(1)
    if (.not. mask(ifirstloc)) ifirstloc=size(mask)+1
  END FUNCTION ifirstloc

  FUNCTION imaxloc_r(arr)
    REAL(8), DIMENSION(:), INTENT(IN) :: arr
    INTEGER :: imaxloc_r
    INTEGER, DIMENSION(1) :: imax
    imax=maxloc(arr(:))
    imaxloc_r=imax(1)
  END FUNCTION imaxloc_r

  FUNCTION imaxloc_i(iarr)
    INTEGER, DIMENSION(:), INTENT(IN) :: iarr
    INTEGER, DIMENSION(1) :: imax
    INTEGER :: imaxloc_i
    imax=maxloc(iarr(:))
    imaxloc_i=imax(1)
  END FUNCTION imaxloc_i

  FUNCTION iminloc(arr)
    REAL(8), DIMENSION(:), INTENT(IN) :: arr
    INTEGER, DIMENSION(1) :: imin
    INTEGER :: iminloc
    imin=minloc(arr(:))
    iminloc=imin(1)
  END FUNCTION iminloc

  SUBROUTINE assert1(n1,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    LOGICAL, INTENT(IN) :: n1
    if (.not. n1) then
       write (*,*) 'nrerror: an assertion failed with this tag:', &
            string
       STOP 'program terminated by assert1'
    end if
  END SUBROUTINE assert1

  SUBROUTINE assert2(n1,n2,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    LOGICAL, INTENT(IN) :: n1,n2
    if (.not. (n1 .and. n2)) then
       write (*,*) 'nrerror: an assertion failed with this tag:', &
            string
       STOP 'program terminated by assert2'
    end if
  END SUBROUTINE assert2

  SUBROUTINE assert3(n1,n2,n3,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    LOGICAL, INTENT(IN) :: n1,n2,n3
    if (.not. (n1 .and. n2 .and. n3)) then
       write (*,*) 'nrerror: an assertion failed with this tag:', &
            string
       STOP 'program terminated by assert3'
    end if
  END SUBROUTINE assert3

  SUBROUTINE assert4(n1,n2,n3,n4,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    LOGICAL, INTENT(IN) :: n1,n2,n3,n4
    if (.not. (n1 .and. n2 .and. n3 .and. n4)) then
       write (*,*) 'nrerror: an assertion failed with this tag:', &
            string
       STOP 'program terminated by assert4'
    end if
  END SUBROUTINE assert4

  SUBROUTINE assert_v(n,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    LOGICAL, DIMENSION(:), INTENT(IN) :: n
    if (.not. all(n)) then
       write (*,*) 'nrerror: an assertion failed with this tag:', &
            string
       STOP 'program terminated by assert_v'
    end if
  END SUBROUTINE assert_v

  FUNCTION assert_eq2(n1,n2,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN) :: n1,n2
    INTEGER :: assert_eq2
    if (n1 == n2) then
       assert_eq2=n1
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eq2'
    end if
  END FUNCTION assert_eq2

  FUNCTION assert_eq3(n1,n2,n3,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN) :: n1,n2,n3
    INTEGER :: assert_eq3
    if (n1 == n2 .and. n2 == n3) then
       assert_eq3=n1
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eq3'
    end if
  END FUNCTION assert_eq3

  FUNCTION assert_eq4(n1,n2,n3,n4,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN) :: n1,n2,n3,n4
    INTEGER :: assert_eq4
    if (n1 == n2 .and. n2 == n3 .and. n3 == n4) then
       assert_eq4=n1
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eq4'
    end if
  END FUNCTION assert_eq4

  FUNCTION assert_eqn(nn,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, DIMENSION(:), INTENT(IN) :: nn
    INTEGER :: assert_eqn
    if (all(nn(2:) == nn(1))) then
       assert_eqn=nn(1)
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eqn'
    end if
  END FUNCTION assert_eqn

  SUBROUTINE nrerror(string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    WRITE(*,*) string
    !PAUSE 'program terminated by nrerror'
    STOP
  END SUBROUTINE nrerror

  FUNCTION arth_r(first,increment,n)
    REAL(8), INTENT(IN) :: first,increment
    INTEGER, INTENT(IN) :: n
    REAL(8), DIMENSION(n) :: arth_r
    INTEGER :: k,k2
    REAL(8) :: temp
    if (n > 0) arth_r(1)=first
    if (n <= NPAR_ARTH) then
       do k=2,n
          arth_r(k)=arth_r(k-1)+increment
       end do
    else
       do k=2,NPAR2_ARTH
          arth_r(k)=arth_r(k-1)+increment
       end do
       temp=increment*NPAR2_ARTH
       k=NPAR2_ARTH
       do
          if (k >= n) exit
          k2=k+k
          arth_r(k+1:min(k2,n))=temp+arth_r(1:min(k,n-k))
          temp=temp+temp
          k=k2
       end do
    end if
  END FUNCTION arth_r

  FUNCTION arth_i(first,increment,n)
    INTEGER, INTENT(IN) :: first,increment,n
    INTEGER, DIMENSION(n) :: arth_i
    INTEGER :: k,k2,temp
    if (n > 0) arth_i(1)=first
    if (n <= NPAR_ARTH) then
       do k=2,n
          arth_i(k)=arth_i(k-1)+increment
       end do
    else
       do k=2,NPAR2_ARTH
          arth_i(k)=arth_i(k-1)+increment
       end do
       temp=increment*NPAR2_ARTH
       k=NPAR2_ARTH
       do
          if (k >= n) exit
          k2=k+k
          arth_i(k+1:min(k2,n))=temp+arth_i(1:min(k,n-k))
          temp=temp+temp
          k=k2
       end do
    end if
  END FUNCTION arth_i

  FUNCTION geop_r(first,factor,n)
    REAL(8), INTENT(IN) :: first,factor
    INTEGER, INTENT(IN) :: n
    REAL(8), DIMENSION(n) :: geop_r
    INTEGER :: k,k2
    REAL(8) :: temp
    if (n > 0) geop_r(1)=first
    if (n <= NPAR_GEOP) then
       do k=2,n
          geop_r(k)=geop_r(k-1)*factor
       end do
    else
       do k=2,NPAR2_GEOP
          geop_r(k)=geop_r(k-1)*factor
       end do
       temp=factor**NPAR2_GEOP
       k=NPAR2_GEOP
       do
          if (k >= n) exit
          k2=k+k
          geop_r(k+1:min(k2,n))=temp*geop_r(1:min(k,n-k))
          temp=temp*temp
          k=k2
       end do
    end if
  END FUNCTION geop_r

  FUNCTION geop_i(first,factor,n)
    INTEGER, INTENT(IN) :: first,factor,n
    INTEGER, DIMENSION(n) :: geop_i
    INTEGER :: k,k2,temp
    if (n > 0) geop_i(1)=first
    if (n <= NPAR_GEOP) then
       do k=2,n
          geop_i(k)=geop_i(k-1)*factor
       end do
    else
       do k=2,NPAR2_GEOP
          geop_i(k)=geop_i(k-1)*factor
       end do
       temp=factor**NPAR2_GEOP
       k=NPAR2_GEOP
       do
          if (k >= n) exit
          k2=k+k
          geop_i(k+1:min(k2,n))=temp*geop_i(1:min(k,n-k))
          temp=temp*temp
          k=k2
       end do
    end if
  END FUNCTION geop_i

  FUNCTION geop_rv(first,factor,n)
    REAL(8), DIMENSION(:), INTENT(IN) :: first,factor
    INTEGER, INTENT(IN) :: n
    REAL(8), DIMENSION(size(first),n) :: geop_rv
    INTEGER :: k,k2
    REAL(8), DIMENSION(size(first)) :: temp
    if (n > 0) geop_rv(:,1)=first(:)
    if (n <= NPAR_GEOP) then
       do k=2,n
          geop_rv(:,k)=geop_rv(:,k-1)*factor(:)
       end do
    else
       do k=2,NPAR2_GEOP
          geop_rv(:,k)=geop_rv(:,k-1)*factor(:)
       end do
       temp=factor**NPAR2_GEOP
       k=NPAR2_GEOP
       do
          if (k >= n) exit
          k2=k+k
          geop_rv(:,k+1:min(k2,n))=geop_rv(:,1:min(k,n-k))*&
               spread(temp,2,size(geop_rv(:,1:min(k,n-k)),2))
          temp=temp*temp
          k=k2
       end do
    end if
  END FUNCTION geop_rv

  RECURSIVE FUNCTION cumsum_r(arr,seed) RESULT(ans)
    REAL(8), DIMENSION(:), INTENT(IN) :: arr
    REAL(8), OPTIONAL, INTENT(IN) :: seed
    REAL(8), DIMENSION(size(arr)) :: ans
    INTEGER :: n,j
    REAL(8) :: sd
    n=size(arr)
    if (n == 0) RETURN
    sd=0.0d0
    if (present(seed)) sd=seed
    ans(1)=arr(1)+sd
    if (n < NPAR_CUMSUM) then
       do j=2,n
          ans(j)=ans(j-1)+arr(j)
       end do
    else
       ans(2:n:2)=cumsum_r(arr(2:n:2)+arr(1:n-1:2),sd)
       ans(3:n:2)=ans(2:n-1:2)+arr(3:n:2)
    end if
  END FUNCTION cumsum_r

  RECURSIVE FUNCTION cumsum_i(arr,seed) RESULT(ans)
    INTEGER, DIMENSION(:), INTENT(IN) :: arr
    INTEGER, OPTIONAL, INTENT(IN) :: seed
    INTEGER, DIMENSION(size(arr)) :: ans
    INTEGER :: n,j,sd
    n=size(arr)
    if (n == 0) RETURN
    sd=0
    if (present(seed)) sd=seed
    ans(1)=arr(1)+sd
    if (n < NPAR_CUMSUM) then
       do j=2,n
          ans(j)=ans(j-1)+arr(j)
       end do
    else
       ans(2:n:2)=cumsum_i(arr(2:n:2)+arr(1:n-1:2),sd)
       ans(3:n:2)=ans(2:n-1:2)+arr(3:n:2)
    end if
  END FUNCTION cumsum_i

  RECURSIVE FUNCTION cumprod(arr,seed) RESULT(ans)
    REAL(8), DIMENSION(:), INTENT(IN) :: arr
    REAL(8), OPTIONAL, INTENT(IN) :: seed
    REAL(8), DIMENSION(size(arr)) :: ans
    INTEGER :: n,j
    REAL(8) :: sd
    n=size(arr)
    if (n == 0) RETURN
    sd=1.0d0
    if (present(seed)) sd=seed
    ans(1)=arr(1)*sd
    if (n < NPAR_CUMPROD) then
       do j=2,n
          ans(j)=ans(j-1)*arr(j)
       end do
    else
       ans(2:n:2)=cumprod(arr(2:n:2)*arr(1:n-1:2),sd)
       ans(3:n:2)=ans(2:n-1:2)*arr(3:n:2)
    end if
  END FUNCTION cumprod

  FUNCTION poly_rr(x,coeffs)
    REAL(8), INTENT(IN) :: x
    REAL(8), DIMENSION(:), INTENT(IN) :: coeffs
    REAL(8) :: poly_rr
    REAL(8) :: pow
    REAL(8), DIMENSION(:), ALLOCATABLE :: vec
    INTEGER :: i,n,nn
    n=size(coeffs)
    if (n <= 0) then
       poly_rr=0.0d0
    else if (n < NPAR_POLY) then
       poly_rr=coeffs(n)
       do i=n-1,1,-1
          poly_rr=x*poly_rr+coeffs(i)
       end do
    else
       allocate(vec(n+1))
       pow=x
       vec(1:n)=coeffs
       do
          vec(n+1)=0.0d0
          nn=ishft(n+1,-1)
          vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
          if (nn == 1) exit
          pow=pow*pow
          n=nn
       end do
       poly_rr=vec(1)
       deallocate(vec)
    end if
  END FUNCTION poly_rr

  FUNCTION poly_rrv(x,coeffs)
    REAL(8), DIMENSION(:), INTENT(IN) :: coeffs,x
    REAL(8), DIMENSION(size(x)) :: poly_rrv
    INTEGER :: i,n,m
    m=size(coeffs)
    n=size(x)
    if (m <= 0) then
       poly_rrv=0.0d0
    else if (m < n .or. m < NPAR_POLY) then
       poly_rrv=coeffs(m)
       do i=m-1,1,-1
          poly_rrv=x*poly_rrv+coeffs(i)
       end do
    else
       do i=1,n
          poly_rrv(i)=poly_rr(x(i),coeffs)
       end do
    end if
  END FUNCTION poly_rrv

  FUNCTION poly_msk_rrv(x,coeffs,mask)
    REAL(8), DIMENSION(:), INTENT(IN) :: coeffs,x
    LOGICAL, DIMENSION(:), INTENT(IN) :: mask
    REAL(8), DIMENSION(size(x)) :: poly_msk_rrv
    poly_msk_rrv=unpack(poly_rrv(pack(x,mask),coeffs),mask,0.0d0)
  END FUNCTION poly_msk_rrv

  RECURSIVE FUNCTION poly_term(a,b) RESULT(u)
    REAL(8), DIMENSION(:), INTENT(IN) :: a
    REAL(8), INTENT(IN) :: b
    REAL(8), DIMENSION(size(a)) :: u
    INTEGER :: n,j
    n=size(a)
    if (n <= 0) RETURN
    u(1)=a(1)
    if (n < NPAR_POLYTERM) then
       do j=2,n
          u(j)=a(j)+b*u(j-1)
       end do
    else
       u(2:n:2)=poly_term(a(2:n:2)+a(1:n-1:2)*b,b*b)
       u(3:n:2)=a(3:n:2)+b*u(2:n-1:2)
    end if
  END FUNCTION poly_term

  FUNCTION outerdiv(a,b)
    REAL(8), DIMENSION(:), INTENT(IN) :: a,b
    REAL(8), DIMENSION(size(a),size(b)) :: outerdiv
    outerdiv = spread(a,dim=2,ncopies=size(b)) / &
         spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerdiv

  FUNCTION outersum(a,b)
    REAL(8), DIMENSION(:), INTENT(IN) :: a,b
    REAL(8), DIMENSION(size(a),size(b)) :: outersum
    outersum = spread(a,dim=2,ncopies=size(b)) + &
         spread(b,dim=1,ncopies=size(a))
  END FUNCTION outersum

  FUNCTION outerdiff_r(a,b)
    REAL(8), DIMENSION(:), INTENT(IN) :: a,b
    REAL(8), DIMENSION(size(a),size(b)) :: outerdiff_r
    outerdiff_r = spread(a,dim=2,ncopies=size(b)) - &
         spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerdiff_r

  FUNCTION outerdiff_i(a,b)
    INTEGER, DIMENSION(:), INTENT(IN) :: a,b
    INTEGER, DIMENSION(size(a),size(b)) :: outerdiff_i
    outerdiff_i = spread(a,dim=2,ncopies=size(b)) - &
         spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerdiff_i

  FUNCTION outerand(a,b)
    LOGICAL, DIMENSION(:), INTENT(IN) :: a,b
    LOGICAL, DIMENSION(size(a),size(b)) :: outerand
    outerand = spread(a,dim=2,ncopies=size(b)) .and. &
         spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerand

  SUBROUTINE scatter_add(dest,source,dest_index)
    REAL(8), DIMENSION(:), INTENT(OUT) :: dest
    REAL(8), DIMENSION(:), INTENT(IN) :: source
    INTEGER, DIMENSION(:), INTENT(IN) :: dest_index
    INTEGER :: m,n,j,i
    n=assert_eq2(size(source),size(dest_index),'scatter_add')
    m=size(dest)
    do j=1,n
       i=dest_index(j)
       if (i > 0 .and. i <= m) dest(i)=dest(i)+source(j)
    end do
  END SUBROUTINE scatter_add

  SUBROUTINE scatter_max(dest,source,dest_index)
    REAL(8), DIMENSION(:), INTENT(OUT) :: dest
    REAL(8), DIMENSION(:), INTENT(IN) :: source
    INTEGER, DIMENSION(:), INTENT(IN) :: dest_index
    INTEGER :: m,n,j,i
    n=assert_eq2(size(source),size(dest_index),'scatter_max')
    m=size(dest)
    do j=1,n
       i=dest_index(j)
       if (i > 0 .and. i <= m) dest(i)=max(dest(i),source(j))
    end do
  END SUBROUTINE scatter_max

  FUNCTION upper_triangle(j,k,extra)
    INTEGER, INTENT(IN) :: j,k
    INTEGER, OPTIONAL, INTENT(IN) :: extra
    LOGICAL, DIMENSION(j,k) :: upper_triangle
    INTEGER :: n
    n=0
    if (present(extra)) n=extra
    upper_triangle=(outerdiff_i(arth_i(1,1,j),arth_i(1,1,k)) < n)
  END FUNCTION upper_triangle

  FUNCTION lower_triangle(j,k,extra)
    INTEGER, INTENT(IN) :: j,k
    INTEGER, OPTIONAL, INTENT(IN) :: extra
    LOGICAL, DIMENSION(j,k) :: lower_triangle
    INTEGER :: n
    n=0
    if (present(extra)) n=extra
    lower_triangle=(outerdiff_i(arth_i(1,1,j),arth_i(1,1,k)) > -n)
  END FUNCTION lower_triangle

  FUNCTION vabs(v)
    REAL(8), DIMENSION(:), INTENT(IN) :: v
    REAL(8) :: vabs
    vabs=sqrt(dot_product(v,v))
  END FUNCTION vabs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! GAMMA SPECIAL FUNCTIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL(8) FUNCTION NR_gamminc(x,a)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: a,x
    IF ((x<0.0d0).OR.(a<0.0d0)) CALL NRERROR('arguments have to be postive: NR_gamminc')
    IF (x<a+1.0d0) THEN
       NR_gamminc=NR_gser(a,x)
    ElSE
       NR_gamminc=1.0d0-NR_gcf(a,x)
    END IF
  END FUNCTION NR_gamminc
  
  REAL(8) FUNCTION NR_gser(a,x,gln)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: a,x
    REAL(8), OPTIONAL, INTENT(OUT) :: gln
    INTEGER, PARAMETER :: ITMAX=1000
    REAL(8), PARAMETER :: EPS=2.220446049250313D-012
    INTEGER :: n
    REAL(8) :: ap,del,summ
    if (x == 0.0d0) then
       NR_gser=0.0d0
       RETURN
    end if
    ap=a
    summ=1.0d0/a
    del=summ
    do n=1,ITMAX
       ap=ap+1.0d0
       del=del*x/ap
       summ=summ+del
       if (abs(del) < abs(summ)*EPS) exit
    end do
    if (n > ITMAX) CALL NRERROR('a too large, ITMAX too small in NR_gser')
    if (present(gln)) then
       gln=LOG_GAMMA(a)
       NR_gser=summ*exp(-x+a*dlog(x)-gln)
    else
       NR_gser=summ*exp(-x+a*dlog(x)-LOG_GAMMA(a))
    end if
  END FUNCTION NR_gser
  
  REAL(8) FUNCTION NR_gcf(a,x,gln)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: a,x
    REAL(8), OPTIONAL, INTENT(OUT) :: gln
    INTEGER, PARAMETER :: ITMAX=1000
    REAL(8), PARAMETER :: EPS=2.220446049250313D-012,FPMIN=1.0D-250/EPS
    INTEGER :: i
    REAL(8) :: an,b,c,d,del,h
    if (x == 0.0d0) then
       NR_gcf=1.0d0
       RETURN
    end if
    b=x+1.0d0-a
    c=1.0d0/FPMIN
    d=1.0d0/b
    h=d
    do i=1,ITMAX
       an=-i*(i-a)
       b=b+2.0d0
       d=an*d+b
       if (abs(d) < FPMIN) d=FPMIN
       c=b+an/c
       if (abs(c) < FPMIN) c=FPMIN
       d=1.0d0/d
       del=d*c
       h=h*del
       if (abs(del-1.0d0) <= EPS) exit
    end do
    if (i > ITMAX) CALL NRERROR('a too large, ITMAX too small in NR_gcf')
    if (present(gln)) then
       gln=LOG_GAMMA(a)
       NR_gcf=exp(-x+a*log(x)-gln)*h
    else
       NR_gcf=exp(-x+a*log(x)-LOG_GAMMA(a))*h
    end if
  END FUNCTION NR_gcf
  
END MODULE nrutil

!MODULE FRPR_PRIVATE_F1DIM
!  IMPLICIT NONE
  
!  INTEGER :: ncom
!  REAL(8), DIMENSION(:), POINTER :: pcom,xicom
  
!CONTAINS  
!  FUNCTION f1dim(x,func)
!    IMPLICIT NONE
!    INTERFACE
!       FUNCTION func(p)
!         IMPLICIT NONE
!         REAL(8), DIMENSION(:), INTENT(IN) :: p
!         REAL(8) :: func
!       END FUNCTION func
!    END INTERFACE
!    REAL(8), INTENT(IN) :: x
!    REAL(8) :: f1dim
!    REAL(8), DIMENSION(:), ALLOCATABLE :: xt
!    ALLOCATE(xt(ncom))
!    xt(:)=pcom(:)+x*xicom(:)
!    f1dim=func(xt)
!    DEALLOCATE(xt)
!  END FUNCTION f1dim
  
!END MODULE FRPR_PRIVATE_F1DIM
