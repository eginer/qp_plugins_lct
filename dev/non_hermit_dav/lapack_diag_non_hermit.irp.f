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

  integer,          intent(in)  :: n
  double precision, intent(in)  :: A(n,n)
  double precision, intent(out) :: WR(n), WI(n), VL(n,n), VR(n,n)

  integer                       :: lda, ldvl, ldvr, LWORK, INFO
  double precision, allocatable :: Atmp(:,:), WORK(:)

  lda  = n
  ldvl = n
  ldvr = n

  allocate( Atmp(n,n) )
  Atmp(1:n,1:n) = A(1:n,1:n)

  allocate(WORK(1))
  LWORK = -1 ! to ask for the optimal size of WORK
  call dgeev('V', 'V', n, Atmp, lda, WR, WI, VL, ldvl, VR, ldvr, WORK, LWORK, INFO)
  if(INFO.gt.0)then
    print*,'dgeev failed !!',INFO
    stop
  endif
  LWORK = max(int(WORK(1)), 1) ! this is the optimal size of WORK 
  deallocate(WORK)

  allocate(WORK(LWORK))

  ! Actual diagonalization 
  call dgeev('V', 'V', n, Atmp, lda, WR, WI, VL, ldvl, VR, ldvr, WORK, LWORK, INFO)
  if(INFO.ne.0) then
    print*,'dgeev failed !!', INFO
    stop
  endif

  deallocate(Atmp, WORK)

end subroutine lapack_diag_non_sym


subroutine non_sym_diag_inv_right(n,A,leigvec,reigvec,n_real_eigv,eigval)
 implicit none
 BEGIN_DOC
! routine which returns the sorted EIGENVALUES and the corresponding LEFT/RIGHT eigenvetors 
!
! of a non hermitian matrix A(n,n)
!
! n_real_eigv is the number of real eigenvalues, which might be smaller than the dimension "n" 
!
! THE LEFT EIGENVECTORS ARE TAKEN AS THE TRANSPOSED OF THE INVERSE MATRIX OF THE RIGHT EIGENVECTORS
!
! THIS GUARANTEES BI-ORTHONORMALITY 
 END_DOC
 integer, intent(in) :: n
 double precision, intent(in) :: A(n,n)
 double precision, intent(out) :: reigvec(n,n),leigvec(n,n),eigval(n)
 double precision, allocatable :: Aw(:,:)
 integer, intent(out) :: n_real_eigv
 print*,'Computing the left/right eigenvectors ...'
 character*1 :: JOBVL,JOBVR
 JOBVL = "V" ! computes the left  eigenvectors 
 JOBVR = "V" ! computes the right eigenvectors 
 double precision, allocatable :: WR(:),WI(:),Vl(:,:),VR(:,:),S(:,:),inv_reigvec(:,:)
 integer :: i,j
 integer :: n_good
 integer, allocatable :: list_good(:), iorder(:)
 double precision :: thr
 thr = 1.d-10
 ! Eigvalue(n) = WR(n) + i * WI(n)
 allocate(WR(n),WI(n),VL(n,n),VR(n,n),Aw(n,n))
 Aw = A
 do i = 1, n
  do j = i+1, n
   if(dabs(Aw(j,j)-Aw(i,i)).lt.thr)then
     Aw(j,j)+= thr
     Aw(i,i)-= thr
     Aw(j,i) = 0.d0
     Aw(i,j) = Aw(j,i)
   endif
  enddo
 enddo
 call lapack_diag_non_sym(n,Aw,WR,WI,VL,VR)
 ! You track the real eigenvalues 
 n_good = 0
 do i = 1, n
  if(dabs(WI(i)).lt.thr)then
   n_good += 1
  else
   print*,'Found an imaginary component to eigenvalue'
   print*,'Re(i) + Im(i)',WR(i),WI(i)
  endif
 enddo
 n_real_eigv = n_good 
 allocate(iorder(n))
 do i = 1, n
   eigval(i) = WR(i)
   iorder(i) = i
 enddo
 ! You sort the real eigenvalues 
 call dsort(eigval,iorder,n_good)
 do i = 1, n
  do j = 1, n
   reigvec(j,i) = VR(j,iorder(i))
   leigvec(j,i) = VL(j,iorder(i))
  enddo
 enddo
 allocate(inv_reigvec(n,n))
 call get_pseudo_inverse(reigvec,n,n,n,inv_reigvec,n,thr)
 do i = 1, n_real_eigv
  do j = 1, n
   leigvec(j,i) = inv_reigvec(i,j)
  enddo
 enddo
 allocate( S(n,n) )

  ! S = VL x VR
  call dgemm( 'T', 'N', n, n, n, 1.d0                              &
            , leigvec, size(leigvec, 1), reigvec, size(reigvec, 1) &
            , 0.d0, S, size(S, 1) )
! do i = 1,n
!  write(*,'(100(F10.5,X))')S(:,i)
! enddo
!call lapack_diag_non_sym(n,S,WR,WI,VL,VR)
!print*,'Eigenvalues of S'
!do i = 1, n
! print*,WR(i),dabs(WI(i))
!enddo
! call get_inv_half_svd(S, n_real_eigv, inv_reigvec)

  double precision :: accu_d,accu_nd
  accu_nd = 0.d0
  accu_d = 0.d0
  do i = 1, n
    do j = 1, n
      if(i==j) then
       accu_d += S(j,i) * S(j,i)
      else
       accu_nd = accu_nd + S(j,i) * S(j,i)
      endif
    enddo
  enddo
  accu_nd = dsqrt(accu_nd)

  print*,'accu_nd = ',accu_nd
  if( accu_nd .lt. 1d-10 ) then
    ! L x R is already bi-orthogonal
    !print *, ' L & T bi-orthogonality: ok'
    return
  else
   print*,'WARNING PB with bi-orthonormality!!'
   print*,'it might come from the fact that you obtained complex eigenvalues'
   print*,'n = ',n
   print*,'n_real_eigv = ',n_real_eigv
   print*,'If the imaginary parts of the eigenvalues are really small (1^-10)'
   print*,'You can use the routine non_hrmt_diag_split_degen which might give you only real eigenvalues'
   stop
  endif
end

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

subroutine lapack_diag_non_sym_right(n, A, WR, WI, VR)

  implicit none

  integer,          intent(in)  :: n
  double precision, intent(in)  :: A(n,n)
  double precision, intent(out) :: WR(n), WI(n), VR(n,n)

  integer                       :: i, lda, ldvl, ldvr, LWORK, INFO
  double precision, allocatable :: Atmp(:,:), WORK(:), VL(:,:)

  lda  = n
  ldvl = 1
  ldvr = n

  allocate( Atmp(n,n), VL(1,1) )
  Atmp(1:n,1:n) = A(1:n,1:n)

  allocate(WORK(1))
  LWORK = -1
  call dgeev('N', 'V', n, Atmp, lda, WR, WI, VL, ldvl, VR, ldvr, WORK, LWORK, INFO)
  if(INFO.gt.0)then
    print*,'dgeev failed !!',INFO
    stop
  endif

  LWORK = max(int(WORK(1)), 1) ! this is the optimal size of WORK 
  deallocate(WORK)

  allocate(WORK(LWORK))

  ! Actual diagonalization 
  call dgeev('N', 'V', n, Atmp, lda, WR, WI, VL, ldvl, VR, ldvr, WORK, LWORK, INFO)
  if(INFO.ne.0) then
    print*,'dgeev failed !!', INFO
    stop
  endif

  deallocate(Atmp, WORK, VL)

  print *, ' JOBL = F'
  print *, ' eigenvalues'
  do i = 1, n
    write(*, '(1000(F16.10,X))') WR(i), WI(i)
  enddo
  print *, ' right eigenvect' 
  do i = 1, n
    write(*, '(1000(F16.10,X))') VR(:,i)
  enddo

end subroutine lapack_diag_non_sym_right



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

subroutine check_degen(n, m, eigval, leigvec, reigvec)

  implicit none
  integer,          intent(in)    :: n, m
  double precision, intent(in)    :: eigval(m)
  double precision, intent(inout) :: leigvec(n,m), reigvec(n,m)
 
  integer                         :: i, j
  double precision                :: ei, ej, de, de_thr, accu_nd
  double precision, allocatable   :: S(:,:)

  de_thr = 1d-7

  do i = 1, m-1
    ei = eigval(i)

    do j = i+1, m
      ej = eigval(j)
      de = dabs(ei - ej)

      if(de .lt. de_thr) then

        leigvec(:,i) = 0.d0
        leigvec(:,j) = 0.d0
        leigvec(i,i) = 1.d0
        leigvec(j,j) = 1.d0

        reigvec(:,i) = 0.d0
        reigvec(:,j) = 0.d0
        reigvec(i,i) = 1.d0
        reigvec(j,j) = 1.d0

      endif

    enddo
  enddo

  ! ---

  allocate( S(m,m) )

  ! S = VL x VR
  call dgemm( 'T', 'N', m, m, n, 1.d0                              &
            , leigvec, size(leigvec, 1), reigvec, size(reigvec, 1) &
            , 0.d0, S, size(S, 1) )

  accu_nd = 0.d0
  do i = 1, m
    do j = 1, m
      if(i==j) cycle
      accu_nd = accu_nd + S(j,i) * S(j,i)
    enddo
  enddo
  accu_nd = dsqrt(accu_nd)

  deallocate( S )

  print *, ' check_degen: L & T bi-orthogonality: ok'
  print *, ' accu_nd = ', accu_nd

  if( accu_nd .lt. 1d-8 ) then
    return
  else
    stop
  endif

end subroutine check_degen

! ---

subroutine rotate_degen_eigvec(n, C)

  implicit none

  integer,          intent(in)    :: n
  double precision, intent(inout) :: C(n,n)

  integer                         :: i
  double precision, allocatable   :: S(:,:), tmp(:,:)
  
  allocate(S(n,n))
  ! S = C.T x C
  call dgemm( 'T', 'N', n, n, n, 1.d0      &
            , C, size(C, 1), C, size(C, 1) &
            , 0.d0, S, size(S, 1) )

  print *, ' eigenvec overlap: '
  do i = 1, n
    write(*, '(1000(F16.10,X))') S(i,:)
  enddo
 
  call get_halfinv_svd(n, S)

  allocate(tmp(n,n))
  ! C <-- C x S^-1/2
  call dgemm( 'N', 'N', n, n, n, 1.d0      &
            , C, size(C, 1), S, size(S, 1) &
            , 0.d0, tmp, size(tmp, 1) )

  C(1:n,1:n) = tmp(1:n,1:n) 

  deallocate(S, tmp)

end subroutine rotate_degen_eigvec

! ---

subroutine get_halfinv_svd(n, S)

  implicit none

  integer,          intent(in)    :: n
  double precision, intent(inout) :: S(n,n)

  integer                         :: num_linear_dependencies
  integer                         :: i, j, k
  double precision                :: accu_d, accu_nd, thresh
  double precision, parameter     :: threshold = 1.d-6
  double precision, allocatable   :: U(:,:), Vt(:,:), D(:)
  double precision, allocatable   :: S0(:,:), Stmp(:,:), Stmp2(:,:)

  allocate( S0(n,n) )
  S0(1:n,1:n) = S(1:n,1:n)

  allocate(U(n,n), Vt(n,n), D(n))
  call svd(S, n, U, n, D, Vt, n, n, n)

  num_linear_dependencies = 0
  do i = 1, n
    if(abs(D(i)) <= threshold) then
      D(i) = 0.d0
      num_linear_dependencies += 1
    else
      ASSERT (D(i) > 0.d0)
      D(i) = 1.d0 / dsqrt(D(i))
    endif
  enddo
  write(*,*) ' linear dependencies', num_linear_dependencies

  S(:,:) = 0.d0
  do k = 1, n
    if(D(k) /= 0.d0) then
      do j = 1, n
        do i = 1, n
          S(i,j) = S(i,j) + U(i,k) * D(k) * Vt(k,j)
        enddo
      enddo
    endif
  enddo
  deallocate(U, D, Vt)

  allocate( Stmp(n,n), Stmp2(n,n) )
  Stmp  = 0.d0
  Stmp2 = 0.d0
  ! S^-1/2 x S
  call dgemm( 'N', 'N', n, n, n, 1.d0        &
            , S, size(S, 1), S0, size(S0, 1) &
            , 0.d0, Stmp, size(Stmp, 1) )
  ! ( S^-1/2 x S ) x S^-1/2
  call dgemm( 'N', 'N', n, n, n, 1.d0            &
            , Stmp, size(Stmp, 1), S, size(S, 1) &
            , 0.d0, Stmp2, size(Stmp2, 1) )

  accu_nd = 0.d0
  accu_d  = 0.d0
  thresh  = 1.d-10
  do i = 1, n
    do j = 1, n
      if(i==j) then
       accu_d += Stmp2(j,i)
      else 
       accu_nd = accu_nd + Stmp2(j,i) * Stmp2(j,i)
      endif
    enddo
  enddo
  accu_nd = dsqrt(accu_nd)
  if( accu_nd.gt.thresh .or. dabs(accu_d-dble(n)).gt.thresh) then
    print*, ' after S^-1/2: sum of off-diag S elements = ', accu_nd
    print*, ' after S^-1/2: sum of     diag S elements = ', accu_d
    do i = 1, n
      write(*,'(1000(F16.10,X))') Stmp2(i,:)
    enddo
    stop
  endif

  deallocate(S0, Stmp, Stmp2)

end subroutine get_halfinv_svd

! ---

subroutine non_sym_diag_inv_right_split_degen(n,A,leigvec,reigvec,n_real_eigv,eigval)
 implicit none
  BEGIN_DOC
  !
  ! routine which returns the sorted REAL EIGENVALUES ONLY and corresponding LEFT/RIGHT eigenvetors 
  !
  ! of a non hermitian matrix A(n,n)
  !
  ! n_real_eigv is the number of real eigenvalues, which might be smaller than the dimension "n" 
  !
  END_DOC

  integer,          intent(in)  :: n
  double precision, intent(in)  :: A(n,n)
  integer,          intent(out) :: n_real_eigv
  double precision, intent(out) :: reigvec(n,n), leigvec(n,n), eigval(n)
  double precision              :: reigvec_tmp(n,n), leigvec_tmp(n,n)

  integer                       :: i, j, k , iteration,l
  integer                       :: n_good
  double precision              :: shift,shift_current
  double precision              :: r,thr
  integer,          allocatable :: list_good(:), iorder_origin(:),iorder(:)
  double precision, allocatable :: WR(:), WI(:), Vl(:,:), VR(:,:),S(:,:)
  double precision, allocatable :: Aw(:,:),diag_elem(:),A_save(:,:)
  double precision, allocatable :: im_part(:),re_part(:),inv_reigvec(:,:)


  print*,'Computing the left/right eigenvectors ...'
  print*,'Using the degeneracy splitting algorithm'


  ! pre-processing the matrix :: sorting by diagonal elements
  allocate(diag_elem(n),iorder_origin(n),A_save(n,n))
  do i = 1, n
   iorder_origin(i) = i
   diag_elem(i) = A(i,i)
  enddo
  call dsort(diag_elem, iorder_origin, n)
  do i = 1, n
   do j = 1, n
    A_save(j,i) = A(iorder_origin(j),iorder_origin(i))
   enddo
  enddo

  shift = 1.d-15
  shift_current = shift
  iteration = 1 
  logical :: good_ortho
  good_ortho = .False.
  do while(n_real_eigv.ne.n.or. .not.good_ortho)
    if(shift.gt.1.d-3)then
     print*,'shift > 1.d-3 !!'
     print*,'Your matrix intrinsically contains complex eigenvalues'
     stop
    endif
    print*,'***** iteration = ',iteration
    print*,'shift = ',shift
    allocate(WR(n), WI(n), VL(n,n), VR(n,n), Aw(n,n))
    Aw = A_save
    do i = 1, n
     do j = 1, n
      if(dabs(Aw(j,i)).lt.shift)then
       Aw(j,i) = 0.d0
      endif
     enddo
    enddo
    call lapack_diag_non_sym(n,Aw,WR,WI,VL,VR)
    allocate(im_part(n),iorder(n))
    do i = 1, n
     im_part(i) = -dabs(WI(i))
     iorder(i) = i
    enddo
    call dsort(im_part, iorder, n)
    shift_current = max(10.d0 * dabs(im_part(1)),shift)
    print*,'Largest imaginary part found in eigenvalues = ',im_part(1)
    print*,'Splitting the degeneracies by ',shift_current
    Aw = A_save
    call split_matrix_degen(Aw,n,shift_current)
    deallocate( im_part, iorder )
    call lapack_diag_non_sym(n,Aw,WR,WI,VL,VR)
    ! You track the real eigenvalues 
    n_good = 0
    do i = 1, n
      if(dabs(WI(i)).lt.1.d-20)then
        n_good += 1
      else
        print*,'Found an imaginary component to eigenvalue'
        print*,'Re(i) + Im(i)',WR(i),WI(i)
      endif
    enddo
    allocate( list_good(n_good), iorder(n_good) )
    n_good = 0
    do i = 1, n
      if(dabs(WI(i)).lt.1.d-20)then
        n_good += 1
        list_good(n_good) = i
        eigval(n_good) = WR(i)
      endif
    enddo
    deallocate( WR, WI )
    
    n_real_eigv = n_good 
    do i = 1, n_good
      iorder(i) = i
    enddo
    
    ! You sort the real eigenvalues 
    call dsort(eigval, iorder, n_good)
    
    reigvec(:,:) = 0.d0 
    leigvec(:,:) = 0.d0 
    do i = 1, n_real_eigv
      do j = 1, n
        reigvec_tmp(j,i) = VR(j,list_good(iorder(i)))
        leigvec_tmp(j,i) = Vl(j,list_good(iorder(i)))
      enddo
    enddo
    
    if(n_real_eigv == n)then
      allocate(S(n,n))
      call check_bi_ortho(reigvec_tmp,leigvec_tmp,n,S,accu_nd)
      print*,'accu_nd = ',accu_nd
      double precision :: accu_nd
      good_ortho = accu_nd .lt. 1.d-10
      deallocate(S)
    endif
 
   deallocate( list_good, iorder )
   deallocate( VL, VR, Aw)
   shift *= 10.d0
   iteration += 1
  enddo

  double precision :: shift_degen
  integer,allocatable ::list_degen(:,:)
  integer :: n_list,n_degen,i_degen
  allocate(list_degen(2,n))
  shift_degen = max(1.d-9,shift_current)
 
  call give_degeneracies(eigval,n,shift_degen,n_list,list_degen)
  print*,'passed the give_degeneracies'
  do i = 1, n
   print*,'e(i) = ',i,eigval(i)
  enddo
  do i = 1, n_list
   i_degen = list_degen(1,i)
   n_degen = list_degen(2,i)
   print*,'i_degen,n_degen',i_degen,n_degen
   double precision, allocatable :: mat_tmp(:,:)
   allocate(inv_reigvec(n_degen,n),mat_tmp(n,n_degen))
   k = 1
   do j = i_degen,i_degen+n_degen-1
    mat_tmp(:,k) =  reigvec_tmp(:,j)
    k+= 1
   enddo 
   thr = 1.d-10
   print*,'in the pseudo_inverse'
   call get_pseudo_inverse(mat_tmp,n,n,n_degen,inv_reigvec,n_degen,thr)
   print*,'passed the pseudo_inverse'
   k = 1
   do j = i_degen,i_degen+n_degen-1
    do l = 1, n
     leigvec_tmp(l,j) = inv_reigvec(k,l)
    enddo
    k+= 1
   enddo
   double precision, allocatable :: test(:,:)
   allocate(test(n,n))
   test = 0.d0
   do l = i_degen,i_degen+n_degen-1
    do j = i_degen,i_degen+n_degen-1
     do k = 1, n
      test(j,l) += leigvec_tmp(k,j) * reigvec_tmp(k,l)
     enddo
    enddo
   enddo
   print*,'id mat'
   do l = 1, n
    write(*,'(100(F16.10,X))')test(:,l)
   enddo
   deallocate(test)
   deallocate(inv_reigvec,mat_tmp)
  enddo
  deallocate(list_degen)
  print*,'passed all the inversions'

  allocate( S(n,n) )
  call check_bi_ortho(reigvec_tmp,leigvec_tmp,n,S,accu_nd)
  print*,'accu_nd = ',accu_nd
  good_ortho = accu_nd .lt. 1.d-10
  deallocate(S)
  print*,'accu_nd = ',accu_nd
  if( accu_nd .gt. 1d-10 ) then
   print*,'WARNING PB with bi-orthonormality!!'
   print*,'it might come from the fact that you obtained complex eigenvalues'
   print*,'n = ',n
   print*,'n_real_eigv = ',n_real_eigv
   print*,'If the imaginary parts of the eigenvalues are really small (1^-10)'
   print*,'You can use the routine non_hrmt_diag_split_degen which might give you only real eigenvalues'
   stop
  endif

  do i = 1, n
   do j = 1, n
    reigvec(iorder_origin(j),i) = reigvec_tmp(j,i)
    leigvec(iorder_origin(j),i) = leigvec_tmp(j,i)
   enddo
  enddo

end
