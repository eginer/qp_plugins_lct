
! ---

subroutine non_hrmt_real_diag_new(n,A,leigvec,reigvec,n_real_eigv,eigval)

  BEGIN_DOC
  !
  ! routine which returns the sorted REAL EIGENVALUES ONLY and corresponding LEFT/RIGHT eigenvetors 
  !
  ! of a non hermitian matrix A(n,n)
  !
  ! n_real_eigv is the number of real eigenvalues, which might be smaller than the dimension "n" 
  !
  END_DOC

  implicit none

  integer,          intent(in)  :: n
  double precision, intent(in)  :: A(n,n)
  integer,          intent(out) :: n_real_eigv
  double precision, intent(out) :: reigvec(n,n), leigvec(n,n), eigval(n)

  integer                       :: i, j
  integer                       :: n_good
  double precision              :: thr
  double precision              :: r
  integer,          allocatable :: list_good(:), iorder(:)
  double precision, allocatable :: WR(:), WI(:), Vl(:,:), VR(:,:)
  double precision, allocatable :: Aw(:,:)


  thr = 1.d-5

  print*,'Computing the left/right eigenvectors ...'

  ! Eigvalue(n) = WR(n) + i * WI(n)
  allocate(WR(n), WI(n), VL(n,n), VR(n,n), Aw(n,n))
  Aw = A
  do i = 1, n
    call RANDOM_NUMBER(r)
    Aw(i,i) += thr * r
  enddo
  call lapack_diag_non_sym_new(n,Aw,WR,WI,VL,VR)
  deallocate( Aw )

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

  allocate( list_good(n_good), iorder(n_good) )
  n_good = 0
  do i = 1, n
    if(dabs(WI(i)).lt.thr)then
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
      reigvec(j,i) = VR(j,list_good(iorder(i)))
      leigvec(j,i) = Vl(j,list_good(iorder(i)))
    enddo
  enddo

  deallocate( list_good, iorder )
  deallocate( VL, VR )

end subroutine non_hrmt_real_diag_new

! ---

subroutine lapack_diag_non_sym_new(n,A,WR,WI,VL,VR)

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

subroutine lapack_diag_non_sym(n,A,WR,WI,VL,VR)
 implicit none
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
 integer, intent(in) :: n
 double precision, intent(in) :: A(n,n)
 double precision, intent(out):: WR(n),WI(n),VL(n,n),VR(n,n)
  character*1 :: JOBVL,JOBVR
  JOBVL = "V" ! computes the left  eigenvectors 
  JOBVR = "V" ! computes the right eigenvectors 
  integer     :: lda,ldvl,ldvr,LWORK,INFO
  double precision, allocatable :: WORK(:)
  integer :: n_good
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
end 


subroutine non_hrmt_real_diag(n,A,leigvec,reigvec,n_real_eigv,eigval)
 implicit none
 BEGIN_DOC
! routine which returns the sorted REAL EIGENVALUES ONLY and corresponding LEFT/RIGHT eigenvetors 
!
! of a non hermitian matrix A(n,n)
!
! n_real_eigv is the number of real eigenvalues, which might be smaller than the dimension "n" 
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
 double precision, allocatable :: WR(:),WI(:),Vl(:,:),VR(:,:)
 integer :: i,j
 integer :: n_good
 integer, allocatable :: list_good(:), iorder(:)
 double precision :: thr
 thr = 1.d-5
 ! Eigvalue(n) = WR(n) + i * WI(n)
 allocate(WR(n),WI(n),VL(n,n),VR(n,n),Aw(n,n))
 Aw = A
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
 allocate(list_good(n_good),iorder(n_good))
 n_good = 0
 do i = 1, n
  if(dabs(WI(i)).lt.thr)then
   n_good += 1
   list_good(n_good) = i
   eigval(n_good) = WR(i)
  endif
 enddo
 n_real_eigv = n_good 
 do i = 1, n_good
  iorder(i) = i
 enddo
 ! You sort the real eigenvalues 
 call dsort(eigval,iorder,n_good)
! print*,'n_real_eigv = ',n_real_eigv
! print*,'n           = ',n
 do i = 1, n_real_eigv
!  print*,i,'eigval(i) = ',eigval(i) 
  do j = 1, n
   reigvec(j,i) = VR(j,list_good(iorder(i)))
   leigvec(j,i) = Vl(j,list_good(iorder(i)))
  enddo
!  write(*,'(X,100(F16.10,X))')reigvec(:,i)
!  write(*,'(X,100(F16.10,X))')leigvec(:,i)
 enddo
end

subroutine lapack_diag_general_non_sym(n,A,B,WR,beta,WI,VL,VR)
 implicit none
  BEGIN_DOC
! You enter with a general non hermitian matrix A(n,n) and another B(n,n)
!
! You get out with the real WR and imaginary part WI of the eigenvalues 
!
! Eigvalue(n) = (WR(n) + i * WI(n))/(beta(n))
!
! And the left VL and right VR eigenvectors 
!
! VL(i,j) = <i|Psi_left(j)>  :: projection on the basis element |i> on the jth left  eigenvector 
!
! VR(i,j) = <i|Psi_right(j)> :: projection on the basis element |i> on the jth right eigenvector 
  END_DOC
 integer, intent(in) :: n
 double precision, intent(in) :: A(n,n),B(n,n)
 double precision, intent(out):: WR(n),WI(n),beta(n),VL(n,n),VR(n,n)
  character*1 :: JOBVL,JOBVR
  JOBVL = "V" ! computes the left  eigenvectors 
  JOBVR = "V" ! computes the right eigenvectors 
  integer     :: lda,ldvl,ldvr,LWORK,INFO
  double precision, allocatable :: WORK(:)
  integer :: n_good
  lda  = n
  ldvl = n
  ldvr = n
  allocate(WORK(1))
  LWORK = -1 ! to ask for the optimal size of WORK
  call dggev('V','V',n,A,lda,B,lda,WR,WI,beta,VL,ldvl,VR,ldvr,WORK,LWORK,INFO)
  if(INFO.gt.0)then
   print*,'dgeev failed !!',INFO
   stop
  endif
  LWORK = max(int(work(1)), 1) ! this is the optimal size of WORK 
  deallocate(WORK)
  allocate(WORK(LWORK))
  ! Actual diagonalization 
  call dggev('V','V',n,A,lda,B,lda,WR,WI,beta,VL,ldvl,VR,ldvr,WORK,LWORK,INFO)
  if(INFO.ne.0)then
   print*,'dgeev failed !!',INFO
   stop
  endif
end 

subroutine non_hrmt_general_real_diag(n,A,B,reigvec,leigvec,n_real_eigv,eigval)
 implicit none
 BEGIN_DOC
! routine which returns the sorted REAL EIGENVALUES ONLY and corresponding LEFT/RIGHT eigenvetors 
!
! of a non hermitian matrix A(n,n) and B(n,n) 
!
! A reigvec = eigval * B * reigvec
!
! (A)^\dagger leigvec = eigval * B * leigvec
!
! n_real_eigv is the number of real eigenvalues, which might be smaller than the dimension "n" 
 END_DOC
 integer, intent(in) :: n
 double precision, intent(in) :: A(n,n),B(n,n)
 double precision, intent(out) :: reigvec(n,n),leigvec(n,n),eigval(n)
 double precision, allocatable :: Aw(:,:),Bw(:,:)
 integer, intent(out) :: n_real_eigv
 print*,'Computing the left/right eigenvectors ...'
 character*1 :: JOBVL,JOBVR
 JOBVL = "V" ! computes the left  eigenvectors 
 JOBVR = "V" ! computes the right eigenvectors 
 double precision, allocatable :: WR(:),WI(:),Vl(:,:),VR(:,:),beta(:)
 integer :: i,j
 integer :: n_good
 integer, allocatable :: list_good(:), iorder(:)
 ! Eigvalue(n) = WR(n) + i * WI(n)
 allocate(WR(n),WI(n),VL(n,n),VR(n,n),Aw(n,n),beta(n),Bw(n,n))
 Aw = A
 Bw = B
 call lapack_diag_general_non_sym(n,A,B,WR,beta,WI,VL,VR)
 ! You track the real eigenvalues 
 n_good = 0
 do i = 1, n
  print*,'beta(i) = ',beta(i)
  if(dabs(WI(i)).lt.1.d-12)then
   n_good += 1
  else
!   print*,'Found an imaginary component to eigenvalue'
!   print*,'Re(i) + Im(i)',WR(i),WI(i)
  endif
 enddo
 allocate(list_good(n_good),iorder(n_good))
 n_good = 0
 do i = 1, n
  if(dabs(WI(i)).lt.1.d-12)then
   n_good += 1
   list_good(n_good) = i
   eigval(n_good) = WR(i)/beta(i)
  endif
 enddo
 n_real_eigv = n_good 
 do i = 1, n_good
  iorder(i) = i
 enddo
 ! You sort the real eigenvalues 
 call dsort(eigval,iorder,n_good)
 print*,'n_real_eigv = ',n_real_eigv
 print*,'n           = ',n
 do i = 1, n_real_eigv
  print*,i,'eigval(i) = ',eigval(i) 
  do j = 1, n
   reigvec(j,i) = VR(j,list_good(iorder(i)))
   leigvec(j,i) = Vl(j,list_good(iorder(i)))
  enddo
 enddo
end

! ---

subroutine non_hrmt_bieig(n, A, leigvec, reigvec, n_real_eigv, eigval)

  BEGIN_DOC
  ! 
  ! routine which returns the sorted REAL EIGENVALUES ONLY and corresponding LEFT/RIGHT eigenvetors 
  ! of a non hermitian matrix A(n,n)
  !
  ! n_real_eigv is the number of real eigenvalues, which might be smaller than the dimension "n" 
  !
  END_DOC

  implicit none
  integer,          intent(in)  :: n
  double precision, intent(in)  :: A(n,n)
  integer,          intent(out) :: n_real_eigv
  double precision, intent(out) :: reigvec(n,n), leigvec(n,n), eigval(n)

  integer                       :: i, j
  integer                       :: n_good
  double precision              :: thr
  double precision              :: accu_nd

  integer,          allocatable :: list_good(:), iorder(:)
  double precision, allocatable :: Aw(:,:)
  double precision, allocatable :: WR(:), WI(:), VL(:,:), VR(:,:)
  double precision, allocatable :: S(:,:)


  ! -------------------------------------------------------------------------------------
  !

  print *, 'Computing the left/right eigenvectors ...'

  allocate( WR(n), WI(n), VL(n,n), VR(n,n), Aw(n,n) )
  Aw(:,:) = A(:,:)

  call lapack_diag_non_sym_new(n, Aw, WR, WI, VL, VR)

  deallocate( Aw )

  !
  ! -------------------------------------------------------------------------------------

  ! ---

  ! -------------------------------------------------------------------------------------
  !                  track & sort the real eigenvalues 

  n_good = 0
  thr    = 1.d-5
  do i = 1, n
    if( dabs(WI(i)).lt.thr ) then
      n_good += 1
    else
      print*, 'Found an imaginary component to eigenvalue on i = ', i
      print*, 'Re(i) + Im(i)', WR(i), WI(i)
    endif
  enddo

  allocate( list_good(n_good), iorder(n_good) )

  n_good = 0
  do i = 1, n
    if( dabs(WI(i)).lt.thr ) then
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
  call dsort(eigval, iorder, n_good)
      
  reigvec(:,:) = 0.d0 
  leigvec(:,:) = 0.d0 
  do i = 1, n_real_eigv
    do j = 1, n
      reigvec(j,i) = VR(j,list_good(iorder(i)))
      leigvec(j,i) = VL(j,list_good(iorder(i)))
    enddo
  enddo

  deallocate( list_good, iorder )
  deallocate( VL, VR )

  !
  ! -------------------------------------------------------------------------------------

  ! ---

  ! -------------------------------------------------------------------------------------
  !                               check bi-orthogonality

  allocate( S(n_real_eigv,n_real_eigv) )

  ! S = VL x VR
  call dgemm( 'T', 'N', n_real_eigv, n_real_eigv, n, 1.d0          &
            , leigvec, size(leigvec, 1), reigvec, size(reigvec, 1) &
            , 0.d0, S, size(S, 1) )

  accu_nd = 0.d0
  do i = 1, n_real_eigv
    do j = 1, n_real_eigv
      if(i==j) cycle
      accu_nd = accu_nd + S(j,i) * S(j,i)
    enddo
  enddo
  accu_nd = dsqrt(accu_nd)

  if( accu_nd .lt. 1d-8 ) then
    ! L x R is already bi-orthogonal

    print *, ' L & T bi-orthogonality: ok'
    deallocate( S )
    return

  else
    ! impose bi-orthogonality 

    print *, ' L & T bi-orthogonality: not imposed yet'
    print *, ' accu_nd = ', accu_nd
    call impose_biorthog_qr( n, n_real_eigv, leigvec, reigvec, S )
    deallocate( S )
  
  endif

  !
  ! -------------------------------------------------------------------------------------

  return

end subroutine non_hrmt_bieig

! ---
subroutine non_hrmt_bieig_random_diag(n, A, leigvec, reigvec, n_real_eigv, eigval)

  BEGIN_DOC
  ! 
  ! routine which returns the sorted REAL EIGENVALUES ONLY and corresponding LEFT/RIGHT eigenvetors 
  ! of a non hermitian matrix A(n,n)
  !
  ! n_real_eigv is the number of real eigenvalues, which might be smaller than the dimension "n" 
  !
  END_DOC
  implicit none
  integer,          intent(in)  :: n
  double precision, intent(in)  :: A(n,n)
  integer,          intent(out) :: n_real_eigv
  double precision, intent(out) :: reigvec(n,n), leigvec(n,n), eigval(n)

  integer                       :: i, j
  integer                       :: n_good
  double precision              :: thr
  double precision              :: accu_nd

  integer,          allocatable :: list_good(:), iorder(:)
  double precision, allocatable :: Aw(:,:)
  double precision, allocatable :: WR(:), WI(:), VL(:,:), VR(:,:)
  double precision, allocatable :: S(:,:)
  double precision :: r


  ! -------------------------------------------------------------------------------------
  !

  print *, 'Computing the left/right eigenvectors ...'
  allocate( WR(n), WI(n), VL(n,n), VR(n,n), Aw(n,n) )

  Aw(:,:) = A(:,:)
  call lapack_diag_non_sym_new(n, Aw, WR, WI, VL, VR)

  thr    = 1.d-12
  double precision, allocatable :: im_part(:)
  n_good = 0
  do i = 1, n
    if( dabs(WI(i)).lt.thr ) then
      n_good += 1
    else
      print*, 'Found an imaginary component to eigenvalue on i = ', i
      print*, 'Re(i) + Im(i)', WR(i), WI(i)
    endif
  enddo
  print*,'n_good = ',n_good
  if(n_good .lt. n)then
   print*,'Removing degeneracies to remove imaginary parts'
   allocate(im_part(n),iorder(n))
   r = 0.d0
   do i = 1, n
     im_part(i) = -dabs(WI(i))
     iorder(i) = i
   enddo
   call dsort(im_part,iorder,n) 
   thr = 10.d0 * dabs(im_part(1))
   print*,'adding random numbers on the diagonal of magnitude ',thr
   Aw(:,:) = A(:,:)
   do i = 1, n
     call RANDOM_NUMBER(r)
     print*,'r = ',r*thr
     Aw(i,i) += thr * r
   enddo
   print*,'Rediagonalizing the matrix with random numbers'
   call lapack_diag_non_sym_new(n, Aw, WR, WI, VL, VR)
   deallocate(im_part,iorder)
  endif
  deallocate( Aw )

  !
  ! -------------------------------------------------------------------------------------

  ! ---

  ! -------------------------------------------------------------------------------------
  !                  track & sort the real eigenvalues 

  n_good = 0
  thr    = 1.d-5
  do i = 1, n
    if( dabs(WI(i)).lt.thr ) then
      n_good += 1
    else
      print*, 'Found an imaginary component to eigenvalue on i = ', i
      print*, 'Re(i) + Im(i)', WR(i), WI(i)
    endif
  enddo
  print*,'n_good = ',n_good
  allocate( list_good(n_good), iorder(n_good) )

  n_good = 0
  do i = 1, n
    if( dabs(WI(i)).lt.thr ) then
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
  call dsort(eigval, iorder, n_good)
      
  reigvec(:,:) = 0.d0 
  leigvec(:,:) = 0.d0 
  do i = 1, n_real_eigv
    do j = 1, n
      reigvec(j,i) = VR(j,list_good(iorder(i)))
      leigvec(j,i) = VL(j,list_good(iorder(i)))
    enddo
  enddo

  deallocate( list_good, iorder )
  deallocate( VL, VR )

  !
  ! -------------------------------------------------------------------------------------

  ! ---

  ! -------------------------------------------------------------------------------------
  !                               check bi-orthogonality

  allocate( S(n_real_eigv,n_real_eigv) )

  ! S = VL x VR
  call dgemm( 'T', 'N', n_real_eigv, n_real_eigv, n, 1.d0          &
            , leigvec, size(leigvec, 1), reigvec, size(reigvec, 1) &
            , 0.d0, S, size(S, 1) )

  accu_nd = 0.d0
  do i = 1, n_real_eigv
    do j = 1, n_real_eigv
      if(i==j) cycle
      accu_nd = accu_nd + S(j,i) * S(j,i)
    enddo
  enddo
  accu_nd = dsqrt(accu_nd)

  if( accu_nd .lt. 1d-8 ) then
    ! L x R is already bi-orthogonal

    print *, ' L & T bi-orthogonality: ok'
    deallocate( S )
    return

  else
    ! impose bi-orthogonality 

    print *, ' L & T bi-orthogonality: not imposed yet'
    print *, ' accu_nd = ', accu_nd
    call impose_biorthog_qr( n, n_real_eigv, leigvec, reigvec, S )
    deallocate( S )
  
  endif

  !
  ! -------------------------------------------------------------------------------------

  return

end 

subroutine non_hrmt_bieig_real_im(n, A, leigvec, reigvec, n_real_eigv, eigval)

  BEGIN_DOC
  ! 
  ! routine which returns the EIGENVALUES sorted the REAL part and corresponding LEFT/RIGHT eigenvetors 
  ! of a non hermitian matrix A(n,n)
  !
  ! n_real_eigv is the number of real eigenvalues, which might be smaller than the dimension "n" 
  !
  END_DOC
  implicit none
  integer,          intent(in)  :: n
  double precision, intent(in)  :: A(n,n)
  integer,          intent(out) :: n_real_eigv
  double precision, intent(out) :: reigvec(n,n), leigvec(n,n), eigval(n)

  integer                       :: i, j
  integer                       :: n_bad
  double precision              :: thr
  double precision              :: accu_nd

  integer,          allocatable :: iorder(:)
  double precision, allocatable :: Aw(:,:)
  double precision, allocatable :: WR(:), WI(:), VL(:,:), VR(:,:)
  double precision, allocatable :: S(:,:)
  double precision :: r


  ! -------------------------------------------------------------------------------------
  !

  print *, 'Computing the left/right eigenvectors ...'
  allocate( WR(n), WI(n), VL(n,n), VR(n,n), Aw(n,n), iorder(n))

  Aw(:,:) = A(:,:)
  call lapack_diag_non_sym_new(n, Aw, WR, WI, VL, VR)

  ! -------------------------------------------------------------------------------------
  !                  track & sort the real eigenvalues 

  i = 1
  thr    = 1.d-5
  n_real_eigv = 0
  do while (i.le.n) 
    if( dabs(WI(i)).gt.thr ) then
      print*, 'Found an imaginary component to eigenvalue on i = ', i
      print*, 'Re(i) , Im(i)  ', WR(i), WI(i)
      iorder(i) = i
      eigval(i) = WR(i)
      i+=1
      print*, 'Re(i+1),Im(i+1)',WR(i), WI(i)
      iorder(i) = i
      eigval(i) = WR(i)
      i+=1
    else  
      n_real_eigv += 1
      iorder(i) = i
      eigval(i) = WR(i)
      i+=1
    endif
  enddo
  call dsort(eigval, iorder, n)
  reigvec(:,:) = 0.d0 
  leigvec(:,:) = 0.d0 
  do i = 1, n
    do j = 1, n
      reigvec(j,i) = VR(j,iorder(i))
      leigvec(j,i) = VL(j,iorder(i))
    enddo
  enddo

  deallocate( iorder )
  deallocate( VL, VR )

  !
  ! -------------------------------------------------------------------------------------

  ! ---

  ! -------------------------------------------------------------------------------------
  !                               check bi-orthogonality

  allocate( S(n,n) )

  ! S = VL x VR
  call dgemm( 'T', 'N', n, n, n, 1.d0          &
            , leigvec, size(leigvec, 1), reigvec, size(reigvec, 1) &
            , 0.d0, S, size(S, 1) )

  accu_nd = 0.d0
  do i = 1, n
    do j = 1, n
      if(i==j) cycle
      accu_nd = accu_nd + S(j,i) * S(j,i)
    enddo
  enddo
  accu_nd = dsqrt(accu_nd)

  if( accu_nd .lt. 1d-8 ) then
    ! L x R is already bi-orthogonal

    print *, ' L & T bi-orthogonality: ok'
    deallocate( S )
    return

  else
    ! impose bi-orthogonality 

    print *, ' L & T bi-orthogonality: not imposed yet'
    print *, ' accu_nd = ', accu_nd
!    call impose_biorthog_qr( n, n, leigvec, reigvec, S )
    deallocate( S )
  
  endif

  !
  ! -------------------------------------------------------------------------------------

  return

end 

subroutine non_hrmt_real_im(n, A, leigvec, reigvec, n_real_eigv, eigval)

  BEGIN_DOC
  ! 
  ! routine which returns the EIGENVALUES sorted the REAL part and corresponding LEFT/RIGHT eigenvetors 
  ! of a non hermitian matrix A(n,n)
  !
  ! n_real_eigv is the number of real eigenvalues, which might be smaller than the dimension "n" 
  !
  END_DOC
  implicit none
  integer,          intent(in)  :: n
  double precision, intent(in)  :: A(n,n)
  integer,          intent(out) :: n_real_eigv
  double precision, intent(out) :: reigvec(n,n), leigvec(n,n), eigval(n)

  integer                       :: i, j
  integer                       :: n_bad
  double precision              :: thr
  double precision              :: accu_nd

  integer,          allocatable :: iorder(:)
  double precision, allocatable :: Aw(:,:)
  double precision, allocatable :: WR(:), WI(:), VL(:,:), VR(:,:)
  double precision, allocatable :: S(:,:)
  double precision :: r,eps


  ! -------------------------------------------------------------------------------------
  !

  print *, 'Computing the left/right eigenvectors ...'
  allocate( WR(n), WI(n), VL(n,n), VR(n,n), Aw(n,n), iorder(n))

  Aw(:,:) = A(:,:)
  call lapack_diag_non_sym_new(n, Aw, WR, WI, VL, VR)
  double precision, allocatable :: im_part(:)
  allocate(im_part(n))
!  do i = 1, n
!   iorder(i) = i
!   im_part(i) = -dabs(WI(i))
!  enddo
!  call dsort(im_part, iorder, n)
!  eps = dabs(im_part(1)) * 10.d0
!  Aw = A
!  do i = 1, n
!    call RANDOM_NUMBER(r)
!    Aw(i,i) += eps * r
!  enddo
!  call lapack_diag_non_sym_new(n, Aw, WR, WI, VL, VR)

  ! -------------------------------------------------------------------------------------
  !                  track & sort the real eigenvalues 

  i = 1
  thr    = 1.d-5
  n_real_eigv = 0
  do while (i.le.n) 
    if( dabs(WI(i)).gt.thr ) then
      print*, 'Found an imaginary component to eigenvalue on i = ', i
      print*, 'Re(i) , Im(i)  ', WR(i), WI(i)
      iorder(i) = i
      eigval(i) = WR(i)
      i+=1
      print*, 'Re(i+1),Im(i+1)',WR(i), WI(i)
      iorder(i) = i
      eigval(i) = WR(i)
      i+=1
    else  
      n_real_eigv += 1
      iorder(i) = i
      eigval(i) = WR(i)
      i+=1
    endif
  enddo
  call dsort(eigval, iorder, n)
  reigvec(:,:) = 0.d0 
  leigvec(:,:) = 0.d0 
  do i = 1, n
    do j = 1, n
      reigvec(j,i) = VR(j,iorder(i))
      leigvec(j,i) = VL(j,iorder(i))
    enddo
  enddo

  deallocate( iorder )
  deallocate( VL, VR )

  !
  ! -------------------------------------------------------------------------------------

  ! ---

  ! -------------------------------------------------------------------------------------
  !                               check bi-orthogonality

  allocate( S(n,n) )

  ! S = VL x VR
  call dgemm( 'T', 'N', n, n, n, 1.d0          &
            , leigvec, size(leigvec, 1), reigvec, size(reigvec, 1) &
            , 0.d0, S, size(S, 1) )

  accu_nd = 0.d0
  do i = 1, n
    do j = 1, n
      if(i==j) cycle
      accu_nd = accu_nd + S(j,i) * S(j,i)
    enddo
  enddo
  accu_nd = dsqrt(accu_nd)

  if( accu_nd .lt. 1d-8 ) then
    ! L x R is already bi-orthogonal

    print *, ' L & T bi-orthogonality: ok'
    deallocate( S )
    return

  else
    ! impose bi-orthogonality 

    print *, ' L & T bi-orthogonality: not ok !'
    print *, ' accu_nd = ', accu_nd
    deallocate( S )
  
  endif

  !
  ! -------------------------------------------------------------------------------------

  return

end 

subroutine impose_biorthog_qr( m, n, Vl, Vr, S )

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
  call dgeqrf( n, n, S, n, TAU, WORK, LWORK, INFO )
  if( INFO.ne.0 ) then
    print*,'dgeqrf failed !!', INFO
    stop
  endif

  LWORK = max(n, int(WORK(1)))
  deallocate(WORK)

  allocate( WORK(LWORK) )
  call dgeqrf( n, n, S, n, TAU, WORK, LWORK, INFO )
  if( INFO.ne.0 ) then
    print*,'dgeqrf failed !!', INFO
    stop
  endif

  ! save the upper triangular R
  allocate( R(n,n) )
  R(:,:) = S(:,:)

  ! get Q
  LWORK = -1
  call dorgqr( n, n, n, S, n, TAU, WORK, LWORK, INFO)
  if( INFO.ne.0 ) then
    print*,'dorgqr failed !!', INFO
    stop
  endif

  LWORK = max(n, int(WORK(1)))
  deallocate(WORK)

  allocate( WORK(LWORK) )
  call dorgqr( n, n, n, S, n, TAU, WORK, LWORK, INFO )
  if( INFO.ne.0 ) then
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
  call dtrtri( "U", "N", n, R, n, INFO )
  if( INFO.ne.0 ) then
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

