
! ---

subroutine lapack_diag_non_sym_new(n, A, WR, WI, VL, VR)

  BEGIN_DOC
  !
  ! You enter with a general non hermitian matrix A(n,n) 
  !
  ! You get out with the real WR and imaginary part WI of the eigenvalues 
  !
  ! Eigvalue(n) = WR(n) + i * WI(n)
  !
  ! And the left VL and right VR eigenvectors 
  !
  ! VL(i,j) = <i|Psi_left(j)>  :: projection on the basis element |i> on the jth left  eigenvector 
  !
  ! VR(i,j) = <i|Psi_right(j)> :: projection on the basis element |i> on the jth right eigenvector 
  !
  END_DOC

  implicit none
  integer,          intent(in)  :: n
  double precision, intent(in)  :: A(n,n)
  double precision, intent(out) :: WR(n), WI(n), VL(n,n), VR(n,n)

  character*1                   :: JOBVL,JOBVR,BALANC,SENSE
  integer                       :: ILO, IHI
  integer                       :: lda, ldvl, ldvr, LWORK, INFO
  double precision              :: ABNRM
  integer,          allocatable :: IWORK(:)
  double precision, allocatable :: WORK(:), SCALE_array(:), RCONDE(:), RCONDV(:)

  JOBVL  = "V" ! computes the left  eigenvectors 
  JOBVR  = "V" ! computes the right eigenvectors 
  BALANC = "B" ! Diagonal scaling and Permutation for optimization
  SENSE  = "B"
  lda  = n
  ldvl = n
  ldvr = n
  allocate(WORK(1),SCALE_array(n),RCONDE(n),RCONDV(n),IWORK(2*n-2))
  LWORK = -1 ! to ask for the optimal size of WORK
  call dgeevx(BALANC,JOBVL,JOBVR,SENSE,&  ! CHARACTERS 
              n,A,lda,                 &  ! MATRIX TO DIAGONALIZE
              WR,WI,                   &  ! REAL AND IMAGINARY PART OF EIGENVALUES 
              VL,ldvl,VR,ldvr,         &  ! LEFT AND RIGHT EIGENVECTORS 
              ILO,IHI,SCALE_array,ABNRM,RCONDE,RCONDV, & ! OUTPUTS OF OPTIMIZATION
              WORK,LWORK,IWORK,INFO)

  !if(INFO.gt.0)then
  ! print*,'dgeev failed !!',INFO
  if( INFO.ne.0 ) then
    print *, 'dgeevx failed !!', INFO
    stop
  endif

  LWORK = max(int(work(1)), 1) ! this is the optimal size of WORK 
  deallocate(WORK)
  allocate(WORK(LWORK))
  ! Actual dnon_hrmt_real_diag_newiagonalization 
  call dgeevx(BALANC,JOBVL,JOBVR,SENSE,&  ! CHARACTERS 
              n,A,lda,                 &  ! MATRIX TO DIAGONALIZE
              WR,WI,                   &  ! REAL AND IMAGINARY PART OF EIGENVALUES 
              VL,ldvl,VR,ldvr,         &  ! LEFT AND RIGHT EIGENVECTORS 
              ILO,IHI,SCALE_array,ABNRM,RCONDE,RCONDV, & ! OUTPUTS OF OPTIMIZATION
              WORK,LWORK,IWORK,INFO)

  !if(INFO.ne.0)then
  ! print*,'dgeev failed !!',INFO
  if( INFO.ne.0 ) then
    print *, 'dgeevx failed !!', INFO
    stop
  endif

  deallocate( WORK, SCALE_array, RCONDE, RCONDV, IWORK )

end subroutine lapack_diag_non_sym_new

! ---

subroutine lapack_diag_non_sym(n, A, WR, WI, VL, VR)


  BEGIN_DOC
  ! You enter with a general non hermitian matrix A(n,n) 
  !
  ! You get out with the real WR and imaginary part WI of the eigenvalues 
  !
  ! Eigvalue(n) = WR(n) + i * WI(n)
  !
  ! And the left VL and right VR eigenvectors 
  !
  ! VL(i,j) = <i|Psi_left(j)>  :: projection on the basis element |i> on the jth left  eigenvector 
  !
  ! VR(i,j) = <i|Psi_right(j)> :: projection on the basis element |i> on the jth right eigenvector 
  END_DOC

  implicit none
  integer, intent(in) :: n
  double precision, intent(in) :: A(n,n)
  double precision, intent(out):: WR(n),WI(n),VL(n,n),VR(n,n)
  character*1 :: JOBVL,JOBVR
  integer     :: lda,ldvl,ldvr,LWORK,INFO
  double precision, allocatable :: WORK(:)
  integer :: n_good

  JOBVL = "V" ! computes the left  eigenvectors 
  JOBVR = "V" ! computes the right eigenvectors 

  lda  = n
  ldvl = n
  ldvr = n
  allocate(WORK(1))
  LWORK = -1 ! to ask for the optimal size of WORK
  call dgeev('V','V',n,A,lda,WR,WI,VL,ldvl,VR,ldvr,WORK,LWORK,INFO)
  if(INFO.gt.0)then
   print*,'dgeev failed !!',INFO
   stop
  endif
  LWORK = max(int(work(1)), 1) ! this is the optimal size of WORK 
  deallocate(WORK)
  allocate(WORK(LWORK))
  ! Actual diagonalization 
  call dgeev('V','V',n,A,lda,WR,WI,VL,ldvl,VR,ldvr,WORK,LWORK,INFO)
  if(INFO.ne.0)then
   print*,'dgeev failed !!',INFO
   stop
  endif

end subroutine lapack_diag_non_sym

! ---

subroutine impose_biorthog_qr(m, n, Vl, Vr, S)

  implicit none 
  integer, intent(in)             :: m, n
  double precision, intent(inout) :: Vl(m,n), Vr(m,n), S(n,n)

  integer                         :: i, j
  integer                         :: LWORK, INFO
  double precision, allocatable   :: TAU(:), WORK(:)
  double precision, allocatable   :: R(:,:), tmp(:,:)

  print *, ' apply QR decomposition ...'

  ! -------------------------------------------------------------------------------------
  !                           QR factorization of S: S = Q x R

  allocate( TAU(n), WORK(1) )

  LWORK = -1
  call dgeqrf(n, n, S, n, TAU, WORK, LWORK, INFO)
  if(INFO .ne. 0) then
    print*,'dgeqrf failed !!', INFO
    stop
  endif

  LWORK = max(n, int(WORK(1)))
  deallocate(WORK)

  allocate( WORK(LWORK) )
  call dgeqrf(n, n, S, n, TAU, WORK, LWORK, INFO)
  if(INFO .ne. 0) then
    print*,'dgeqrf failed !!', INFO
    stop
  endif

  ! save the upper triangular R
  allocate( R(n,n) )
  R(:,:) = S(:,:)

  ! get Q
  LWORK = -1
  call dorgqr(n, n, n, S, n, TAU, WORK, LWORK, INFO)
  if(INFO .ne. 0) then
    print*,'dorgqr failed !!', INFO
    stop
  endif

  LWORK = max(n, int(WORK(1)))
  deallocate(WORK)

  allocate( WORK(LWORK) )
  call dorgqr(n, n, n, S, n, TAU, WORK, LWORK, INFO)
  if(INFO .ne. 0) then
    print*,'dorgqr failed !!', INFO
    stop
  endif

  deallocate( WORK, TAU )

  !
  ! -------------------------------------------------------------------------------------

  ! ---

  ! -------------------------------------------------------------------------------------
  !                               get bi-orhtog left & right vectors:
  !                                           Vr' = Vr x inv(R) 
  !                                           Vl' = inv(Q) x Vl =  Q.T   x Vl 

  ! Q.T x Vl, where Q = S

  allocate( tmp(n,m) )
  call dgemm( 'T', 'T', n, m, n, 1.d0        &
            , S, size(S, 1), Vl, size(Vl, 1) &
            , 0.d0, tmp, size(tmp, 1) )

  do i = 1, n
    do j = 1, m
      Vl(j,i) = tmp(i,j)
    enddo
  enddo
  deallocate(tmp)

  ! ---

  ! inv(R) 
  !print *, ' inversing upper triangular matrix ...'
  call dtrtri("U", "N", n, R, n, INFO)
  if(INFO .ne. 0) then
    print*,'dtrtri failed !!', INFO
    stop
  endif
  !print *, ' inversing upper triangular matrix OK' 

  do i = 1, n-1
    do j = i+1, n
      R(j,i) = 0.d0
    enddo
  enddo

  ! Vr x inv(R) 
  allocate( tmp(m,n) )
  call dgemm( 'N', 'N', m, n, n, 1.d0        &
            , Vr, size(Vr, 1), R, size(R, 1) &
            , 0.d0, tmp, size(tmp, 1) )
  deallocate( R )

  do i = 1, n
    do j = 1, m
      Vr(j,i) = tmp(j,i)
    enddo
  enddo
  deallocate(tmp)


  return
end subroutine impose_biorthog_qr

! ---

subroutine impose_biorthog_lu(m, n, Vl, Vr, S)

  implicit none 
  integer, intent(in)             :: m, n
  double precision, intent(inout) :: Vl(m,n), Vr(m,n), S(n,n)

  integer                         :: i, j
  integer                         :: INFO
  double precision                :: nrm
  integer,          allocatable   :: IPIV(:)
  double precision, allocatable   :: L(:,:), tmp(:,:), vectmp(:)
  !double precision, allocatable   :: T(:,:), ll(:,:), rr(:,:), tt(:,:)

  !allocate( T(n,n) )
  !T(:,:) = S(:,:)

  print *, ' apply LU decomposition ...'

  ! -------------------------------------------------------------------------------------
  !                           LU factorization of S: S = P x L x U

  allocate( IPIV(n) )

  call dgetrf(n, n, S, n, IPIV, INFO)
  if(INFO .ne. 0) then
    print*, 'dgetrf failed !!', INFO
    stop
  endif

  ! check | S - P x L x U |
  !allocate( ll(n,n), rr(n,n), tmp(n,n) )
  !ll = S
  !rr = S
  !do i = 1, n-1
  !  ll(i,i) = 1.d0
  !  do j = i+1, n
  !    ll(i,j) = 0.d0
  !    rr(j,i) = 0.d0
  !  enddo
  !enddo
  !ll(n,n) = 1.d0
  !call dgemm( 'N', 'N', n, n, n, 1.d0          &
  !          , ll, size(ll, 1), rr, size(rr, 1) &
  !          , 0.d0, tmp, size(tmp, 1) )
  ! deallocate(ll, rr)
  !allocate( vectmp(n) )
  !do j = n-1, 1, -1
  !  i = IPIV(j)
  !  if(i.ne.j) then
  !    print *, j, i
  !    vectmp(:) = tmp(i,:)
  !    tmp(i,:)  = tmp(j,:)
  !    tmp(j,:)  = vectmp(:)
  !  endif
  !enddo
  !deallocate( vectmp )
  !nrm = 0.d0
  !do i = 1, n
  !  do j = 1, n
  !    nrm += dabs(tmp(j,i) - T(j,i))
  !  enddo
  !enddo
  !deallocate( tmp )
  !print*, '|L.T x R - S| =', nrm
  !stop

  ! ------
  ! inv(L) 
  ! ------

  allocate( L(n,n) )
  L(:,:) = S(:,:)

  call dtrtri("L", "U", n, L, n, INFO)
  if(INFO .ne. 0) then
    print*,  'dtrtri failed !!', INFO
    stop
  endif
  do i = 1, n-1
    L(i,i) = 1.d0
    do j = i+1, n
      L(i,j) = 0.d0
    enddo
  enddo
  L(n,n) = 1.d0

  ! ------
  ! inv(U) 
  ! ------
  
  call dtrtri("U", "N", n, S, n, INFO)
  if(INFO .ne. 0) then
    print*,  'dtrtri failed !!', INFO
    stop
  endif

  do i = 1, n-1
    do j = i+1, n
      S(j,i) = 0.d0
    enddo
  enddo

  !
  ! -------------------------------------------------------------------------------------

  ! ---

  ! -------------------------------------------------------------------------------------
  !                               get bi-orhtog left & right vectors:
  !                                           Vr' = Vr x inv(U) 
  !                                           Vl' = inv(L) x inv(P) x Vl

  ! inv(P) x Vl
  allocate( vectmp(n) )
  do j = n-1, 1, -1
    i = IPIV(j)
    if(i.ne.j) then
      vectmp(:) = L(:,j)
      L(:,j)    = L(:,i)
      L(:,i)    = vectmp(:)
    endif
  enddo
  deallocate( vectmp )

  ! Vl'
  allocate( tmp(m,n) )
  call dgemm( 'N', 'T', m, n, n, 1.d0        &
            , Vl, size(Vl, 1), L, size(L, 1) &
            , 0.d0, tmp, size(tmp, 1) )
  deallocate(L)

  Vl = tmp
  deallocate(tmp)

  ! ---

  ! Vr x inv(U) 
  allocate( tmp(m,n) )
  call dgemm( 'N', 'N', m, n, n, 1.d0        &
            , Vr, size(Vr, 1), S, size(S, 1) &
            , 0.d0, tmp, size(tmp, 1) )
  Vr = tmp
  deallocate(tmp)

  !allocate( tmp(n,n) )
  !call dgemm( 'T', 'N', n, n, m, 1.d0          &
  !          , Vl, size(Vl, 1), Vr, size(Vr, 1) &
  !          , 0.d0, tmp, size(tmp, 1) )
  !nrm = 0.d0
  !do i = 1, n
  !  do j = 1, n
  !    nrm += dabs(tmp(j,i))
  !  enddo
  !enddo
  !deallocate( tmp )
  !print*, '|L.T x R| =', nrm
  !stop

  return
end subroutine impose_biorthog_lu

! ---

subroutine impose_biorthog_inv(m, n, Vl, Vr)

  implicit none 
  integer,          intent(in)  :: m, n
  double precision, intent(in)  :: Vr(m,n)
  double precision, intent(out) :: Vl(m,n)

  integer                       :: i, j
  integer                       :: INFO, LWORK
  integer,          allocatable :: IPIV(:)
  double precision, allocatable :: WORK(:), Mtmp(:,:)


  print *, ' compute inv of right eigenvectors...'

  if(m .ne. n) then
    print*, 'm not equal n !'
    stop
  endif

  allocate( Mtmp(m,m) )
  Mtmp = Vr

  allocate( IPIV(m) )
  do i = 1, m
    IPIV(i) = i
    !Mtmp(i,i) = Mtmp(i,i) + 1d-10
  enddo

  allocate(WORK(1))
  LWORK = -1
  call dgetri(m, Mtmp, size(Mtmp, 1), IPIV, WORK, LWORK, INFO)
  if(INFO .ne. 0) then
    print *, ' problem in dgetri 1 !', INFO
    stop
  endif
  LWORK = WORK(1)
  deallocate( WORK )

  allocate(WORK(LWORK))
  call dgetri(m, Mtmp, size(Mtmp, 1), IPIV, WORK, LWORK, INFO)
  if(INFO .ne. 0) then
    print *, ' problem in dgetri 2 !', INFO
    stop
  endif
  deallocate(WORK)

  do i = 1, m
    do j = 1, m
      Vl(j,i) = Mtmp(i,j)
    enddo
  enddo
  deallocate( Mtmp )

  return
end subroutine impose_biorthog_inv

! ---

subroutine check_EIGVEC(n, m, A, eigval, leigvec, reigvec)

  implicit none
  integer,          intent(in)  :: n, m
  double precision, intent(in)  :: A(n,n), eigval(m), leigvec(n,m), reigvec(n,m)
 
  integer                       :: i, j
  double precision              :: tmp, tmp_all, tmp_nrm
  double precision              :: V_nrm
  double precision, allocatable :: Mtmp(:,:)
 
  allocate( Mtmp(n,m) )
  
  ! ---

  print *, ' check right eigvec : '

  call dgemm( 'N', 'N', n, m, n, 1.d0                  &
            , A, size(A, 1), reigvec, size(reigvec, 1) &
            , 0.d0, Mtmp, size(Mtmp, 1) )

  V_nrm   = 0.d0
  tmp_nrm = 0.d0
  tmp_all = 0.d0
  do j = 1, m
    tmp = 0.d0
    do i = 1, n
      tmp     = tmp     + dabs(Mtmp(i,j) - eigval(j) * reigvec(i,j))
      tmp_nrm = tmp_nrm + dabs(Mtmp(i,j))
      V_nrm   = V_nrm   + dabs(reigvec(i,j))
    enddo
    tmp_all = tmp_all + tmp
    print *, j, tmp
  enddo
  print *, ' err estim = ', tmp_all/tmp_nrm
  print *, ' CR norm   = ', V_nrm 

  Mtmp = 0.d0

  ! ---

  print *, ' check left eigvec : '

  call dgemm( 'T', 'N', n, m, n, 1.d0                  &
            , A, size(A, 1), leigvec, size(leigvec, 1) &
            , 0.d0, Mtmp, size(Mtmp, 1) )

  V_nrm   = 0.d0
  tmp_nrm = 0.d0
  tmp_all = 0.d0
  do j = 1, m
    tmp = 0.d0
    do i = 1, n
      tmp     = tmp     + dabs(Mtmp(i,j) - eigval(j) * leigvec(i,j))
      tmp_nrm = tmp_nrm + dabs(Mtmp(i,j))
      V_nrm   = V_nrm   + dabs(leigvec(i,j))
    enddo
    tmp_all = tmp_all + tmp
    print *, j, tmp
  enddo
  print *, ' err estim = ', tmp_all/tmp_nrm
  print *, ' CL norm   = ', V_nrm 

  ! ---

  deallocate( Mtmp )

end subroutine check_EIGVEC

! ---




