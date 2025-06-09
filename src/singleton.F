module singleton
  !-----------------------------------------------------------------------------
  ! Multivariate Fast Fourier Transform
  !
  ! Fortran 90 (ELF90) Implementation of Singleton's mixed-radix algorithm,
  ! RC Singleton, Stanford Research Institute, Sept. 1968.
  !
  ! Adapted from fftn.c, translated from Fortran 66 to C by Mark Olesen and
  ! John Beale.
  !
  ! Fourier transforms can be computed either in place, using assumed size
  ! arguments, or by generic function, using assumed shape arguments.
  !
  !
  ! Public:
  !
  !   fftkind                              kind parameter of complex arguments
  !                                        and function results.
  !
  !   fft(array, dim, inv)                 generic transform function
  !    COMPLEX(fftkind), DIMENSION(:,...,:), INTENT(IN)           :: array
  !    INTEGER,          DIMENSION(:),       INTENT(IN),  OPTIONAL:: dim
  !    LOGICAL,                              INTENT(IN),  OPTIONAL:: inv
  !
  !   fftn(array, shape, dim, inv, stat)   in place transform subroutine
  !    COMPLEX(fftkind), DIMENSION(*), INTENT(INOUT)        :: array
  !    INTEGER,          DIMENSION(:), INTENT(IN)           :: shape
  !    INTEGER,          DIMENSION(:), INTENT(IN),  OPTIONAL:: dim
  !    LOGICAL,                        INTENT(IN),  OPTIONAL:: inv
  !    INTEGER,                        INTENT(OUT), OPTIONAL:: stat
  !
  !
  ! Formal Parameters:
  !
  !   array    The complex array to be transformed. array can be of arbitrary
  !            rank (i.e. up to seven).
  !
  !   shape    With subroutine fftn, the shape of the array to be transformed
  !            has to be passed separately, since fftradix - the internal trans-
  !            formation routine - will treat array always as one dimensional.
  !            The product of elements in shape must be the number of
  !            elements in array.
  !            Although passing array with assumed shape would have been nicer,
  !            I prefered assumed size in order to prevent the compiler from
  !            using a copy-in-copy-out mechanism. That would generally be
  !            necessary with fftn passing array to fftradix and with fftn
  !            being prepared for accepting non consecutive array sections.
  !            Using assumed size, it's up to the user to pass an array argu-
  !            ment, that can be addressed as continous one dimensional array
  !            without copying. Otherwise, transformation will not really be
  !            performed in place.
  !            On the other hand, since the rank of array and the size of
  !            shape needn't match, fftn is appropriate for handling more than
  !            seven dimensions.
  !            As far as function fft is concerned all this doesn't matter,
  !            because the argument will be copied anyway. Thus no extra
  !            shape argument is needed for fft.
  !
  ! Optional Parameters:
  !
  !   dim      One dimensional integer array, containing the dimensions to be
  !            transformed. Default is (/1,...,N/) with N being the rank of
  !            array, i.e. complete transform. dim can restrict transformation
  !            to a subset of available dimensions. Its size must not exceed the
  !            rank of array or the size of shape respectivly.
  !
  !   inv      If .true., inverse transformation will be performed. Default is
  !            .false., i.e. forward transformation.
  !
  !   stat     If present, a system dependent nonzero status value will be
  !            returned in stat, if allocation of temporary storage failed.
  !            For functions, the integer variable status is used.
  !
  ! Scaling:
  !
  !   Transformation results will always be scaled by the square root of the
  !   product of sizes of each dimension in dim. (See examples below)
  !
  !
  ! Examples:
  !
  !   Let A be a L*M*N three dimensional complex array. Then
  !
  !     result = fft(A)
  !
  !   will produce a three dimensional transform, scaled by sqrt(L*M*N), while
  !
  !     call fftn(A, SHAPE(A))
  !
  !   will do the same in place.
  !
  !     result = fft(A, dim=(/1,3/))
  !
  !   will transform with respect to the first and the third dimension, scaled
  !   by sqrt(L*N).
  !
  !     result = fft(fft(A), inv=.true.)
  !
  !   should (approximately) reproduce A.
  !   With B having the same shape as A
  !
  !     result = fft(fft(A) * CONJG(fft(B)), inv=.true.)
  !
  !   will correlate A and B.
  !
  !
  ! Remarks:
  !
  !   Following changes have been introduced with respect to fftn.c:
  !   - complex arguments and results are of type complex, rather than
  !     real an imaginary part separately.
  !   - increment parameter (magnitude of isign) has been dropped,
  !     inc is always one, direction of transform is given by inv.     
  !   - maxf and maxp have been dropped. The amount of temporary storage
  !     needed is determined by the fftradix routine. Both fftn and fft
  !     can handle any size of array. (Maybe they take a lot of time and
  !     memory, but they will do it)
  !
  !   Redesigning fftradix in a way, that it handles assumed shape arrays
  !   would have been desirable. However, I found it rather hard to do this
  !   in an efficient way. Problems were:
  !   - to prevent stride multiplications when indexing arrays. At least our
  !     compiler was not clever enough to discover that in fact additions
  !     would do the job as well. On the other hand, I haven't been clever
  !     enough to find an implementation using array operations.
  !   - fftradix is rather large and different versions would be necessaray
  !     for each possible rank of array.
  !   Consequently, in place transformation still needs the argument stored
  !   in a consecutive bunch of memory and can't be performed on array
  !   sections like A(100:199:-3, 50:1020). Calling fftn with such sections
  !   will most probably imply copy-in-copy-out. However, the function fft
  !   works with everything it gets and should be convenient to use.
  !
  !   To enable this module to be used with ELF90 it appears to be necessary
  !   to allocate a 1-D work array into which the multi-dimensional array is
  !   copied, and then to copy the results back from the 1-D array to the
  !   multi-dimensional array ft.
  !
  !   Unfortunately, ELF90 will not allow a function to return more than one
  !   output variable.   The variable `stat' has been dropped from the function
  !   arguments.   Users should examine the value of the variable `status'
  !   instead.   This is a PUBLIC variable declared in this module.
  !
  ! Michael Steffens, 09.12.96, <Michael.Steffens@mbox.muk.uni-hannover.de>
  ! ELF90-compatible version by Alan Miller, 29 April 1997 & 6 June 1997
  ! amiller @ bigpond.net.au
  ! Restructured fftradix for better optimization by M. Steffens, 4 June 1997
  !-----------------------------------------------------------------------------
  implicit none

  private
  public:: fft, fftn, fftkind, status

  integer, parameter:: fftkind = KIND(0.0) !--- adjust here for other precisions
  integer, save     :: status = 0 !--- shifted to here as ELF90 does not allow
                                  !    arguments to be INTENT(OUT)

  real(fftkind), parameter:: sin60 = 0.86602540378443865_fftkind
  real(fftkind), parameter:: cos72 = 0.30901699437494742_fftkind
  real(fftkind), parameter:: sin72 = 0.95105651629515357_fftkind
  real(fftkind), parameter:: pi    = 3.14159265358979323_fftkind

  interface fft
     module procedure fft1d
     module procedure fft2d
     module procedure fft3d
     module procedure fft4d
     module procedure fft5d
     module procedure fft6d
     module procedure fft7d
  end interface


contains


  function fft1d(array, dim, inv) result(ft)
    !--- formal parameters
    complex(fftkind), dimension(:), intent(in)           :: array
    integer,          dimension(:), intent(in),  optional:: dim
    logical,                        intent(in),  optional:: inv
    !--- function result
    complex(fftkind), dimension(SIZE(array, 1)):: ft

    ft = array
    call fftn(ft, SHAPE(array), dim, inv = inv, stat = status)

    return
  end function fft1d


  function fft2d(array, dim, inv) result(ft)
    !--- formal parameters
    complex(fftkind), dimension(:,:), intent(in)           :: array
    integer,          dimension(:),   intent(in),  optional:: dim
    logical,                          intent(in),  optional:: inv
    !--- function result
    complex(fftkind), dimension(SIZE(array, 1), SIZE(array, 2)):: ft

    !--- Allocate 1-D array to be used by routine fftn
    complex(fftkind), dimension(:), allocatable :: work

    allocate( work(SIZE(array)) )
    work = RESHAPE(array, (/ SIZE(array) /))
    call fftn(work, SHAPE(array), dim, inv, stat = status)
    ft = RESHAPE(work, (/ SIZE(array, 1), SIZE(array, 2) /))
    deallocate( work )

    return
  end function fft2d


  function fft3d(array, dim, inv) result(ft)
    !--- formal parameters
    complex(fftkind), dimension(:,:,:), intent(in)           :: array
    integer,          dimension(:),     intent(in),  optional:: dim
    logical,                            intent(in),  optional:: inv
    !--- function result
    complex(fftkind), &
         dimension(SIZE(array, 1), SIZE(array, 2), SIZE(array, 3)):: ft

    !--- Allocate 1-D array to be used by routine fftn
    complex(fftkind), dimension(:), allocatable :: work

    allocate( work(SIZE(array)) )
    work = RESHAPE(array, (/ SIZE(array) /))
    call fftn(work, SHAPE(array), dim, inv, stat = status)
    ft = RESHAPE(work, (/ SIZE(array, 1), SIZE(array, 2), SIZE(array, 3) /))
    deallocate( work )

    return
  end function fft3d


  function fft4d(array, dim, inv) result(ft)
    !--- formal parameters
    complex(fftkind), dimension(:,:,:,:), intent(in)           :: array
    integer,          dimension(:),       intent(in),  optional:: dim
    logical,                              intent(in),  optional:: inv
    !--- function result
    complex(fftkind), dimension( &
         SIZE(array, 1), SIZE(array, 2), SIZE(array, 3), SIZE(array, 4)):: ft

    !--- Allocate 1-D array to be used by routine fftn
    complex(fftkind), dimension(:), allocatable :: work

    allocate( work(SIZE(array)) )
    work = RESHAPE(array, (/ SIZE(array) /))
    call fftn(work, SHAPE(array), dim, inv, stat = status)
    ft = RESHAPE(work, (/ SIZE(array, 1), SIZE(array, 2), SIZE(array, 3), &
                          SIZE(array, 4) /))
    deallocate( work )

    return
  end function fft4d


  function fft5d(array, dim, inv) result(ft)
    !--- formal parameters
    complex(fftkind), dimension(:,:,:,:,:), intent(in)           :: array
    integer,          dimension(:),         intent(in),  optional:: dim
    logical,                                intent(in),  optional:: inv
    !--- function result
    complex(fftkind), dimension( &
         SIZE(array, 1), SIZE(array, 2), SIZE(array, 3), SIZE(array, 4), &
         SIZE(array, 5)):: ft

    !--- Allocate 1-D array to be used by routine fftn
    complex(fftkind), dimension(:), allocatable :: work

    allocate( work(SIZE(array)) )
    work = RESHAPE(array, (/ SIZE(array) /))
    call fftn(work, SHAPE(array), dim, inv, stat = status)
    ft = RESHAPE(work, (/ SIZE(array, 1), SIZE(array, 2), SIZE(array, 3), &
                          SIZE(array, 4), SIZE(array, 5) /))
    deallocate( work )

    return
  end function fft5d


  function fft6d(array, dim, inv) result(ft)
    !--- formal parameters
    complex(fftkind), dimension(:,:,:,:,:,:), intent(in)           :: array
    integer,          dimension(:),           intent(in),  optional:: dim
    logical,                                  intent(in),  optional:: inv
    !--- function result
    complex(fftkind), dimension( &
         SIZE(array, 1), SIZE(array, 2), SIZE(array, 3), SIZE(array, 4), &
         SIZE(array, 5), SIZE(array, 6)):: ft

    !--- Allocate 1-D array to be used by routine fftn
    complex(fftkind), dimension(:), allocatable :: work

    allocate( work(SIZE(array)) )
    work = RESHAPE(array, (/ SIZE(array) /))
    call fftn(work, SHAPE(array), dim, inv, stat = status)
    ft = RESHAPE(work, (/ SIZE(array, 1), SIZE(array, 2), SIZE(array, 3), &
                          SIZE(array, 4), SIZE(array, 5), SIZE(array, 6) /))
    deallocate( work )

    return
  end function fft6d


  function fft7d(array, dim, inv) result(ft)
    !--- formal parameters
    complex(fftkind), dimension(:,:,:,:,:,:,:), intent(in)           :: array
    integer,          dimension(:),             intent(in),  optional:: dim
    logical,                                    intent(in),  optional:: inv
    !--- function result
    complex(fftkind), dimension( &
         SIZE(array, 1), SIZE(array, 2), SIZE(array, 3), SIZE(array, 4), &
         SIZE(array, 5), SIZE(array, 6), SIZE(array, 7)):: ft

    !--- Allocate 1-D array to be used by routine fftn
    complex(fftkind), dimension(:), allocatable :: work

    allocate( work(SIZE(array)) )
    work = RESHAPE(array, (/ SIZE(array) /))
    call fftn(work, SHAPE(array), dim, inv, stat = status)
    ft = RESHAPE(work, (/ SIZE(array, 1), SIZE(array, 2), SIZE(array, 3), &
                          SIZE(array, 4), SIZE(array, 5), SIZE(array, 6), &
                          SIZE(array, 7) /))
    deallocate( work )

    return
  end function fft7d


  subroutine fftn(array, shape, dim, inv, stat)
    !--- formal parameters
    complex(fftkind), dimension(:), intent(in out)       :: array
    integer,          dimension(:), intent(in)           :: shape
    integer,          dimension(:), intent(in),  optional:: dim
    logical,                        intent(in),  optional:: inv
    integer,                        intent(out), optional:: stat
    !--- local arrays
    integer, dimension(SIZE(shape)):: d
    !--- local scalars
    logical      :: inverse
    integer      :: i, ndim, ntotal
    real(fftkind):: scale

    !--- optional parameter settings
    if (PRESENT(inv)) then
       inverse = inv
    else
       inverse = .FALSE.
    end if
    if (PRESENT(dim)) then
       ndim = MIN(SIZE(dim), SIZE(d))
       d(1:ndim) = dim(1:ndim)
    else
       ndim = SIZE(d)
       d = (/(i, i = 1, SIZE(d))/)
    end if

    ntotal = PRODUCT(shape)
    scale = SQRT(1.0_fftkind / PRODUCT(shape(d(1:ndim))))
    do i = 1, ntotal
       array(i) = CMPLX(real(array(i)) * scale, AIMAG(array(i)) * scale, &
            KIND=fftkind)
    end do

    do i = 1, ndim
       call fftradix(array, ntotal, shape(d(i)), PRODUCT(shape(1:d(i))), &
            inverse)
    end do

    if (PRESENT(stat)) stat = status
    return
  end subroutine fftn


  subroutine fftradix(array, ntotal, npass, nspan, inv)
    !--- formal parameters
    complex(fftkind), dimension(:), intent(in out)       :: array
    integer,                        intent(in)           :: ntotal, npass, nspan
    logical,                        intent(in)           :: inv
    !--- local arrays
    integer,          dimension(BIT_SIZE(0))     :: factor
    !--- local scalars
    integer         :: maxfactor, nfactor, nsquare
    !--- intrinsics used

    if (npass <= 1) return

    call factorize(npass, factor, nfactor, nsquare)

    maxfactor = MAXVAL(factor(:nfactor))

    call transform(array, ntotal, npass, nspan, factor, nfactor, inv)
    call permute(array, ntotal, npass, nspan, &
            factor, nfactor, nsquare, maxfactor)

    return


  contains


    subroutine factorize(npass, factor, nfactor, nsquare)
      !--- formal parameters
      integer,               intent(in) :: npass
      integer, dimension(:), intent(out):: factor
      integer,               intent(out):: nfactor, nsquare
      !--- local scalars
      integer:: j, jj, k

      nfactor = 0
      k = npass
      do WHILE (MOD(k, 16) == 0) 
         nfactor = nfactor + 1
         factor(nfactor) = 4
         k = k / 16
      end do
      j = 3
      jj = 9
      do
         do WHILE (MOD(k, jj) == 0)
            nfactor = nfactor + 1
            factor(nfactor) = j
            k = k / jj
         end do
         j = j + 2
         jj = j * j
         if (jj > k) exit
      end do
      if (k <= 4) then
         nsquare = nfactor
         factor(nfactor + 1) = k
         if (k /= 1) nfactor = nfactor + 1
      else 
         if (k - ISHFT(k / 4, 2) == 0) then
            nfactor = nfactor + 1
            factor(nfactor) = 2
            k = k / 4
         end if
         nsquare = nfactor
         j = 2
         do
            if (MOD(k, j) == 0) then
               nfactor = nfactor + 1
               factor(nfactor) = j
               k = k / j
            end if
            j = ISHFT((j + 1) / 2, 1) + 1
            if (j > k) exit
         end do
      end if
      if (nsquare > 0) then
         j = nsquare
         do
            nfactor = nfactor + 1
            factor(nfactor) = factor(j)
            j = j - 1
            if (j==0) exit
         end do
      end if

      return
    end subroutine factorize


    subroutine transform(array, ntotal, npass, nspan, &
         factor, nfactor, inv) !-- compute fourier transform
      !--- formal parameters
      complex(fftkind), dimension(:), intent(in out):: array
      integer,                        intent(in)    :: ntotal, npass, nspan
      integer,          dimension(:), intent(in)    :: factor
      integer,                        intent(in)    :: nfactor
      logical,                        intent(in)    :: inv
      !--- local scalars
      integer         :: ii, ispan
      integer         :: j, jc, jf, jj
      integer         :: k, kk, kspan, k1, k2, k3, k4
      integer         :: nn, nt
      real(fftkind)   :: s60, c72, s72, pi2, radf
      real(fftkind)   :: c1, s1, c2, s2, c3, s3, cd, sd, ak
      complex(fftkind):: cc, cj, ck, cjp, cjm, ckp, ckm
      !--- local arrays
      complex(fftkind), dimension(:), allocatable  :: ctmp
      real(fftkind),    dimension(:), allocatable  :: sine, cosine
      integer         :: maxfactor

      maxfactor = MAXVAL(factor(:nfactor))
      allocate(ctmp(maxfactor), sine(maxfactor), cosine(maxfactor), STAT=status)
      if (status /= 0) return

      c72 = cos72
      if (inv) then
         s72 = sin72
         s60 = sin60
         pi2 = pi
      else
         s72 = -sin72
         s60 = -sin60
         pi2 = -pi
      end if

      nt = ntotal
      nn = nt - 1
      kspan = nspan
      jc = nspan / npass
      radf = pi2 * jc
      pi2 = pi2 * 2.0_fftkind !-- use 2 PI from here on

      ii = 0
      jf = 0
      do
         sd = radf / kspan
         cd = SIN(sd)
         cd = 2.0_fftkind * cd * cd
         sd = SIN(sd + sd)
         kk = 1
         ii = ii + 1

         select case (factor(ii))
         case (2)
            !-- transform for factor of 2 (including rotation factor)
            kspan = kspan / 2
            k1 = kspan + 2
            do
               do
                  k2 = kk + kspan
                  ck = array(k2)
                  array(k2) = array(kk)-ck
                  array(kk) = array(kk) + ck
                  kk = k2 + kspan
                  if (kk > nn) exit
               end do
               kk = kk - nn
               if (kk > jc) exit
            end do
            if (kk > kspan) then
               deallocate(ctmp, sine, cosine, STAT=status)
               return
            end if
            do
               c1 = 1.0_fftkind - cd
               s1 = sd
               do
                  do
                     do
                        k2 = kk + kspan
                        ck = array(kk) - array(k2)
                        array(kk) = array(kk) + array(k2)
                        array(k2) = ck * CMPLX(c1, s1, KIND=fftkind)
                        kk = k2 + kspan
                        if (kk >= nt) exit
                     end do
                     k2 = kk - nt
                     c1 = -c1
                     kk = k1 - k2
                     if (kk <= k2) exit
                  end do
                  ak = c1 - (cd * c1 + sd * s1)
                  s1 = sd * c1 - cd * s1 + s1
                  c1 = 2.0_fftkind - (ak * ak + s1 * s1)
                  s1 = s1 * c1
                  c1 = c1 * ak
                  kk = kk + jc
                  if (kk >= k2) exit
               end do
               k1 = k1 + 1 + 1
               kk = (k1 - kspan) / 2 + jc
               if (kk > jc + jc) exit
            end do

         case (4) !-- transform for factor of 4
            ispan = kspan
            kspan = kspan / 4

            do
               c1 = 1.0_fftkind
               s1 = 0.0_fftkind
               do
                  do
                     k1 = kk + kspan
                     k2 = k1 + kspan
                     k3 = k2 + kspan
                     ckp = array(kk) + array(k2)
                     ckm = array(kk) - array(k2)
                     cjp = array(k1) + array(k3)
                     cjm = array(k1) - array(k3)
                     array(kk) = ckp + cjp
                     cjp = ckp - cjp
                     if (inv) then
                        ckp = ckm + CMPLX(-AIMAG(cjm), real(cjm), KIND=fftkind)
                        ckm = ckm + CMPLX(AIMAG(cjm), -real(cjm), KIND=fftkind)
                     else
                        ckp = ckm + CMPLX(AIMAG(cjm), -real(cjm), KIND=fftkind)
                        ckm = ckm + CMPLX(-AIMAG(cjm), real(cjm), KIND=fftkind)
                     end if
                     !-- avoid useless multiplies
                     if (s1 == 0.0_fftkind) then
                        array(k1) = ckp
                        array(k2) = cjp
                        array(k3) = ckm
                     else
                        array(k1) = ckp * CMPLX(c1, s1, KIND=fftkind)
                        array(k2) = cjp * CMPLX(c2, s2, KIND=fftkind)
                        array(k3) = ckm * CMPLX(c3, s3, KIND=fftkind)
                     end if
                     kk = k3 + kspan
                     if (kk > nt) exit
                  end do

                  c2 = c1 - (cd * c1 + sd * s1)
                  s1 = sd * c1 - cd * s1 + s1
                  c1 = 2.0_fftkind - (c2 * c2 + s1 * s1)
                  s1 = s1 * c1
                  c1 = c1 * c2
                  !-- values of c2, c3, s2, s3 that will get used next time
                  c2 = c1 * c1 - s1 * s1
                  s2 = 2.0_fftkind * c1 * s1
                  c3 = c2 * c1 - s2 * s1
                  s3 = c2 * s1 + s2 * c1
                  kk = kk - nt + jc
                  if (kk > kspan) exit
               end do
               kk = kk - kspan + 1
               if (kk > jc) exit
            end do
            if (kspan == jc) then
               deallocate(ctmp, sine, cosine, STAT=status)
               return
            end if

         case default
            !-- transform for odd factors
            k = factor(ii)
            ispan = kspan
            kspan = kspan / k

            select case (k)
            case (3) !-- transform for factor of 3 (optional code)
               do
                  do
                     k1 = kk + kspan
                     k2 = k1 + kspan
                     ck = array(kk)
                     cj = array(k1) + array(k2)
                     array(kk) = ck + cj
                     ck = ck - CMPLX( &
                          0.5_fftkind * real (cj), &
                          0.5_fftkind * AIMAG(cj), &
                          KIND=fftkind)
                     cj = CMPLX( &
                          (real (array(k1)) - real (array(k2))) * s60, &
                          (AIMAG(array(k1)) - AIMAG(array(k2))) * s60, &
                          KIND=fftkind)
                     array(k1) = ck + CMPLX(-AIMAG(cj), real(cj), KIND=fftkind)
                     array(k2) = ck + CMPLX(AIMAG(cj), -real(cj), KIND=fftkind)
                     kk = k2 + kspan
                     if (kk >= nn) exit
                  end do
                  kk = kk - nn
                  if (kk > kspan) exit
               end do

            case (5) !-- transform for factor of 5 (optional code)
               c2 = c72 * c72 - s72 * s72
               s2 = 2.0_fftkind * c72 * s72
               do
                  do
                     k1 = kk + kspan
                     k2 = k1 + kspan
                     k3 = k2 + kspan
                     k4 = k3 + kspan
                     ckp = array(k1) + array(k4)
                     ckm = array(k1) - array(k4)
                     cjp = array(k2) + array(k3)
                     cjm = array(k2) - array(k3)
                     cc = array(kk)
                     array(kk) = cc + ckp + cjp
                     ck = CMPLX(real(ckp) * c72, AIMAG(ckp) * c72, &
                          KIND=fftkind) + &
                          CMPLX(real(cjp) * c2,  AIMAG(cjp) * c2,  &
                          KIND=fftkind) + cc
                     cj = CMPLX(real(ckm) * s72, AIMAG(ckm) * s72, &
                          KIND=fftkind) + &
                          CMPLX(real(cjm) * s2,  AIMAG(cjm) * s2,  &
                          KIND=fftkind)
                     array(k1) = ck + CMPLX(-AIMAG(cj), real(cj), KIND=fftkind)
                     array(k4) = ck + CMPLX(AIMAG(cj), -real(cj), KIND=fftkind)
                     ck = CMPLX(real(ckp) * c2,  AIMAG(ckp) * c2,  &
                          KIND=fftkind) + &
                          CMPLX(real(cjp) * c72, AIMAG(cjp) * c72, &
                          KIND=fftkind) + cc
                     cj = CMPLX(real(ckm) * s2,  AIMAG(ckm) * s2,  &
                          KIND=fftkind) - &
                          CMPLX(real(cjm) * s72, AIMAG(cjm) * s72, &
                          KIND=fftkind)
                     array(k2) = ck + CMPLX(-AIMAG(cj), real(cj), KIND=fftkind)
                     array(k3) = ck + CMPLX(AIMAG(cj), -real(cj), KIND=fftkind)
                     kk = k4 + kspan
                     if (kk >= nn) exit
                  end do
                  kk = kk - nn
                  if (kk > kspan) exit
               end do

            case default
               if (k /= jf) then
                  jf = k
                  s1 = pi2 / k
                  c1 = COS(s1)
                  s1 = SIN(s1)
                  cosine (jf) = 1.0_fftkind
                  sine (jf) = 0.0_fftkind
                  j = 1
                  do
                     cosine (j) = cosine (k) * c1 + sine (k) * s1
                     sine (j) = cosine (k) * s1 - sine (k) * c1
                     k = k-1
                     cosine (k) = cosine (j)
                     sine (k) = -sine (j)
                     j = j + 1
                     if (j >= k) exit
                  end do
               end if
               do
                  do
                     k1 = kk
                     k2 = kk + ispan
                     cc = array(kk)
                     ck = cc
                     j = 1
                     k1 = k1 + kspan
                     do
                        k2 = k2 - kspan
                        j = j + 1
                        ctmp(j) = array(k1) + array(k2)
                        ck = ck + ctmp(j)
                        j = j + 1
                        ctmp(j) = array(k1) - array(k2)
                        k1 = k1 + kspan
                        if (k1 >= k2) exit
                     end do
                     array(kk) = ck
                     k1 = kk
                     k2 = kk + ispan
                     j = 1
                     do
                        k1 = k1 + kspan
                        k2 = k2 - kspan
                        jj = j
                        ck = cc
                        cj = (0.0_fftkind, 0.0_fftkind)
                        k = 1
                        do
                           k = k + 1
                           ck = ck + CMPLX( &
                                real (ctmp(k)) * cosine(jj), &
                                AIMAG(ctmp(k)) * cosine(jj), KIND=fftkind)
                           k = k + 1
                           cj = cj + CMPLX( &
                                real (ctmp(k)) * sine(jj), &
                                AIMAG(ctmp(k)) * sine(jj), KIND=fftkind)
                           jj = jj + j
                           if (jj > jf) jj = jj - jf
                           if (k >= jf) exit
                        end do
                        k = jf - j
                        array(k1) = ck + CMPLX(-AIMAG(cj), real(cj), &
                             KIND=fftkind)
                        array(k2) = ck + CMPLX(AIMAG(cj), -real(cj), &
                             KIND=fftkind)
                        j = j + 1
                        if (j >= k) exit
                     end do
                     kk = kk + ispan
                     if (kk > nn) exit
                  end do
                  kk = kk - nn
                  if (kk > kspan) exit
               end do

            end select
            !--  multiply by rotation factor (except for factors of 2 and 4)
            if (ii == nfactor) then
               deallocate(ctmp, sine, cosine, STAT=status)
               return
            end if

            kk = jc + 1
            do
               c2 = 1.0_fftkind - cd
               s1 = sd
               do
                  c1 = c2
                  s2 = s1
                  kk = kk + kspan
                  do
                     do
                        array(kk) = CMPLX(c2, s2, KIND=fftkind) * array(kk)
                        kk = kk + ispan
                        if (kk > nt) exit
                     end do
                     ak = s1 * s2
                     s2 = s1 * c2 + c1 * s2
                     c2 = c1 * c2 - ak
                     kk = kk - nt + kspan
                     if (kk > ispan) exit
                  end do
                  c2 = c1 - (cd * c1 + sd * s1)
                  s1 = s1 + sd * c1 - cd * s1
                  c1 = 2.0_fftkind - (c2 * c2 + s1 * s1)
                  s1 = s1 * c1
                  c2 = c2 * c1
                  kk = kk - ispan + jc
                  if (kk > kspan) exit
               end do
               kk = kk - kspan + jc + 1
               if (kk > jc + jc) exit
            end do

         end select
      end do

      deallocate(ctmp, sine, cosine, STAT=status)
      return
    end subroutine transform


    subroutine permute(array, ntotal, npass, nspan, &
         factor, nfactor, nsquare, maxfactor)
      !--- formal parameters
      complex(fftkind), dimension(:), intent(in out):: array
      integer,                        intent(in)    :: ntotal, npass, nspan
      integer,          dimension(:), intent(in out):: factor
      integer,                        intent(in)    :: nfactor, nsquare
      integer,                        intent(in)    :: maxfactor
      !--- local scalars
      integer         :: ii, ispan
      integer         :: j, jc, jj
      integer         :: k, kk, kspan, kt, k1, k2, k3
      integer         :: nn, nperm, nt
      complex(fftkind):: ck
      !--- local arrays
      complex(fftkind), dimension(:), allocatable  :: ctmp
      integer         , dimension(:), allocatable  :: perm

      allocate(ctmp(maxfactor), STAT=status)
      if (status /= 0) return

      if (nfactor - ISHFT(nsquare, 1) > 0) then
         nperm = MAX(nfactor + 1, PRODUCT(factor(nsquare+1: nfactor-nsquare)) - 1)
      else
         nperm = nfactor + 1
      end if
      allocate(perm(nperm), STAT=status)
      if (status /= 0) return

      !--  permute the results to normal order---done in two stages
      !--  permutation for square factors of n

      nt = ntotal
      nn = nt - 1
      kt = nsquare
      kspan = nspan
      jc = nspan / npass

      perm (1) = nspan
      if (kt > 0) then
         k = kt + kt + 1
         if (nfactor < k) k = k - 1
         j = 1
         perm (k + 1) = jc
         do
            perm (j + 1) = perm (j) / factor(j)
            perm (k) = perm (k + 1) * factor(j)
            j = j + 1
            k = k - 1
            if (j >= k) exit
         end do
         k3 = perm (k + 1)
         kspan = perm (2)
         kk = jc + 1
         k2 = kspan + 1
         j = 1

         if (npass /= ntotal) then
            permute_multi: do
               do
                  do
                     k = kk + jc
                     do
                        !-- swap array(kk) <> array(k2)
                        ck = array(kk)
                        array(kk) = array(k2)
                        array(k2) = ck
                        kk = kk + 1
                        k2 = k2 + 1
                        if (kk >= k) exit
                     end do
                     kk = kk + nspan - jc
                     k2 = k2 + nspan - jc
                     if (kk >= nt) exit
                  end do
                  kk = kk - nt + jc
                  k2 = k2 - nt + kspan
                  if (k2 >= nspan) exit
               end do
               do
                  do
                     k2 = k2 - perm (j)
                     j = j + 1
                     k2 = perm (j + 1) + k2
                     if (k2 <= perm (j)) exit
                  end do
                  j = 1
                  do
                     if (kk < k2) cycle permute_multi
                     kk = kk + jc
                     k2 = k2 + kspan
                     if (k2 >= nspan) exit
                  end do
                  if (kk >= nspan) exit
               end do
               exit
            end do permute_multi
         else
            permute_single: do
               do
                  !-- swap array(kk) <> array(k2)
                  ck = array(kk)
                  array(kk) = array(k2)
                  array(k2) = ck
                  kk = kk + 1
                  k2 = k2 + kspan
                  if (k2 >= nspan) exit
               end do
               do
                  do
                     k2 = k2 - perm (j)
                     j = j + 1
                     k2 = perm (j + 1) + k2
                     if (k2 <= perm (j)) exit
                  end do
                  j = 1
                  do
                     if (kk < k2) cycle permute_single
                     kk = kk + 1
                     k2 = k2 + kspan
                     if (k2 >= nspan) exit
                  end do
                  if (kk >= nspan) exit
               end do
               exit
            end do permute_single
         end if
         jc = k3
      end if

      if (ISHFT(kt, 1) + 1 >= nfactor) then
         deallocate(perm, ctmp)
         return
      end if

      ispan = perm (kt + 1)
      !-- permutation for square-free factors of n
      j = nfactor - kt
      factor(j + 1) = 1
      do
         factor(j) = factor(j) * factor(j+1)
         j = j - 1
         if (j == kt) exit
      end do
      kt = kt + 1
      nn = factor(kt) - 1
      j = 0
      jj = 0
      do
         k = kt + 1
         k2 = factor(kt)
         kk = factor(k)
         j = j + 1
         if (j > nn) exit !-- exit infinite loop
         jj = jj + kk
         do WHILE (jj >= k2)
            jj = jj - k2
            k2 = kk
            k = k + 1
            kk = factor(k)
            jj = jj + kk
         end do
         perm (j) = jj
      end do
      !--  determine the permutation cycles of length greater than 1
      j = 0
      do
         do
            j = j + 1
            kk = perm(j)
            if (kk >= 0) exit
         end do
         if (kk /= j) then
            do
               k = kk
               kk = perm (k)
               perm (k) = -kk
               if (kk == j) exit
            end do
            k3 = kk
         else
            perm (j) = -j
            if (j == nn) exit !-- exit infinite loop
         end if
      end do
      !--  reorder a and b, following the permutation cycles
      do
         j = k3 + 1
         nt = nt - ispan
         ii = nt - 1 + 1
         if (nt < 0) exit !-- exit infinite loop
         do
            do
               j = j-1
               if (perm(j) >= 0) exit
            end do
            jj = jc
            do
               kspan = jj
               if (jj > maxfactor) kspan = maxfactor
               jj = jj - kspan
               k = perm(j)
               kk = jc * k + ii + jj
               k1 = kk + kspan
               k2 = 0
               do
                  k2 = k2 + 1
                  ctmp(k2) = array(k1)
                  k1 = k1 - 1
                  if (k1 == kk) exit
               end do
               do
                  k1 = kk + kspan
                  k2 = k1 - jc * (k + perm(k))
                  k = -perm(k)
                  do
                     array(k1) = array(k2)
                     k1 = k1 - 1
                     k2 = k2 - 1
                     if (k1 == kk) exit
                  end do
                  kk = k2
                  if (k == j) exit
               end do
               k1 = kk + kspan
               k2 = 0
               do
                  k2 = k2 + 1
                  array(k1) = ctmp(k2)
                  k1 = k1 - 1
                  if (k1 == kk) exit
               end do
               if (jj == 0) exit
            end do
            if (j == 1) exit
         end do
      end do

      deallocate(perm, ctmp)
      return
    end subroutine permute

  end subroutine fftradix

end module singleton
