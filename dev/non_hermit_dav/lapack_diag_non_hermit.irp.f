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
!    if(Aw(j,i) * A(i,j) .lt.0d0  )then
!     if(dabs(Aw(j,i) * A(i,j)).lt.thr**(1.5d0))then
!      print*,Aw(j,j),Aw(i,i)
!      print*,Aw(j,i) , A(i,j)
      Aw(j,i) = 0.d0
      Aw(i,j) = Aw(j,i)
!     endif
!    endif
   endif
  enddo
 enddo
 call lapack_diag_non_sym(n,Aw,WR,WI,VL,VR)
 ! You track the real eigenvalues 
 n_good = 0
! do i = 1, n
!  write(*,'(100(F16.12,X))')A(:,i)
! enddo
 do i = 1, n
  print*,'Im part of lambda = ',dabs(WI(i))
  if(dabs(WI(i)).lt.thr)then
   n_good += 1
  else
   print*,'Found an imaginary component to eigenvalue'
   print*,'Re(i) + Im(i)',WR(i),WI(i)
   write(*,'(100(F10.5,X))')VR(:,i)
   write(*,'(100(F10.5,X))')VR(:,i+1)
   write(*,'(100(F10.5,X))')VL(:,i)
   write(*,'(100(F10.5,X))')VL(:,i+1)
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
 do i = 1, n_real_eigv
  do j = 1, n
   reigvec(j,i) = VR(j,list_good(iorder(i)))
   leigvec(j,i) = VL(j,list_good(iorder(i)))
  enddo
 enddo
 allocate(inv_reigvec(n_real_eigv,n_real_eigv))
! call get_pseudo_inverse(reigvec,n_real_eigv,n_real_eigv,n_real_eigv,inv_reigvec,n_real_eigv,thr)
! do i = 1, n_real_eigv
!  do j = 1, n
!   leigvec(j,i) = inv_reigvec(i,j)
!  enddo
! enddo
 allocate( S(n_real_eigv,n_real_eigv) )

  ! S = VL x VR
  call dgemm( 'T', 'N', n_real_eigv, n_real_eigv, n_real_eigv, 1.d0                              &
            , leigvec, size(leigvec, 1), reigvec, size(reigvec, 1) &
            , 0.d0, S, size(S, 1) )
   do i = 1,n_real_eigv
    write(*,'(100(F10.5,X))')S(:,i)
   enddo
! call lapack_diag_non_sym(n,S,WR,WI,VL,VR)
! print*,'Eigenvalues of S'
! do i = 1, n
!  print*,WR(i),dabs(WI(i))
! enddo
  call dgemm( 'T', 'N', n_real_eigv, n_real_eigv, n_real_eigv, 1.d0                              &
            , leigvec, size(leigvec, 1), reigvec, size(reigvec, 1) &
            , 0.d0, S, size(S, 1) )
! call get_inv_half_svd(S, n_real_eigv, inv_reigvec)

  double precision :: accu_d,accu_nd
  accu_nd = 0.d0
  accu_d = 0.d0
  do i = 1, n_real_eigv
    do j = 1, n_real_eigv
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
   print*,'PB with bi-orthonormality!!'
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
  double precision              :: shift,shift_current
  double precision              :: r,thr
  integer,          allocatable :: list_good(:), iorder(:)
  double precision, allocatable :: WR(:), WI(:), Vl(:,:), VR(:,:)
  double precision, allocatable :: Aw(:,:)
  double precision, allocatable :: im_part(:)



  print*,'Computing the left/right eigenvectors ...'

  ! Eigvalue(n) = WR(n) + i * WI(n)
  shift = 1.d-10
  do while(n_real_eigv.ne.n.or.shift.gt.1.d-3)
   allocate(WR(n), WI(n), VL(n,n), VR(n,n), Aw(n,n))
   Aw = A
   call lapack_diag_non_sym(n,Aw,WR,WI,VL,VR)
   allocate(im_part(n), iorder(n))
   do i = 1, n
    im_part(i) = -dabs(WI(i))
    iorder(i) = i
   enddo
   shift_current = max(10.d0 * dabs(im_part(1)),shift)
   print*,'adding random number of magnitude ',shift_current
   Aw = A
   do i = 1, n
     call RANDOM_NUMBER(r)
     Aw(i,i) += shift_current * r
   enddo
   deallocate( im_part, iorder )
   call lapack_diag_non_sym(n,Aw,WR,WI,VL,VR)
 
   ! You track the real eigenvalues 
   thr = 1.d-10
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
   deallocate( VL, VR, Aw)
   shift *= 10.d0
  enddo
  if(shift.gt.1.d-3)then
   print*,'shift > 1.d-3 !!'
   print*,'Your matrix intrinsically contains complex eigenvalues'
  endif

end subroutine non_hrmt_real_diag_new

! ---


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
 do i = 1, n_real_eigv
  do j = 1, n
   reigvec(j,i) = VR(j,list_good(iorder(i)))
   leigvec(j,i) = Vl(j,list_good(iorder(i)))
  enddo
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

  double precision, allocatable :: Atmp(:,:)

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

  call lapack_diag_non_sym(n, Aw, WR, WI, VL, VR)
  !call lapack_diag_non_sym_new(n, Aw, WR, WI, VL, VR)

  deallocate( Aw )

  !
  ! -------------------------------------------------------------------------------------

  ! ---

  ! -------------------------------------------------------------------------------------
  !                  track & sort the real eigenvalues 

  n_good = 0
  thr    = 1.d-10
  do i = 1, n
    if(dabs(WI(i)) .lt. thr) then
      n_good += 1
    else
      print*, 'Found an imaginary component to eigenvalue on i = ', i
      print*, 'Re(i) + Im(i)', WR(i), WI(i)
    endif
  enddo

  allocate(list_good(n_good), iorder(n_good))

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

  !print *, ' check_EIGVEC before QR:'
  !call check_EIGVEC(n, n_real_eigv, A, eigval, leigvec, reigvec)

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
    print *, ' accu_nd = ', accu_nd
    deallocate( S )
    return

  else
    ! impose bi-orthogonality 

    print *, ' L & T bi-orthogonality: not imposed yet'
    print *, ' accu_nd = ', accu_nd
    !call impose_biorthog_qr(n, n_real_eigv, leigvec, reigvec, S)
    call impose_biorthog_lu(n, n_real_eigv, leigvec, reigvec, S)
    deallocate( S )

    !print *, ' check_EIGVEC after QR:'
    !call check_EIGVEC(n, n_real_eigv, A, eigval, leigvec, reigvec)
  
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

! ---

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
  double precision :: r

  ! -------------------------------------------------------------------------------------
  !

  print *, 'Computing the left/right eigenvectors ...'
  allocate( WR(n), WI(n), VL(n,n), VR(n,n), Aw(n,n), iorder(n))

  Aw(:,:) = A(:,:)
   do i = 1, n
     call RANDOM_NUMBER(r)
     Aw(i,i) += 10.d-10* r
   enddo
  call lapack_diag_non_sym(n, Aw, WR, WI, VL, VR)

  ! -------------------------------------------------------------------------------------
  !                  track & sort the real eigenvalues 

  i = 1
  thr    = 1.d-15
  n_real_eigv = 0
  do while (i.le.n) 
!    print*,i,dabs(WI(i))
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

  !print *, ' check_EIGVEC :'
  !call check_EIGVEC(n, n_real_eigv, A, eigval, leigvec, reigvec)

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

  deallocate( S )

end subroutine non_hrmt_real_im

! ---

subroutine non_hrmt_generalized_real_im(n, A, B, leigvec, reigvec, n_real_eigv, eigval)

  BEGIN_DOC
  ! 
  ! routine which returns the EIGENVALUES sorted the REAL part and corresponding LEFT/RIGHT eigenvetors 
  ! for A R = lambda B R and A^\dagger L = lambda B^\dagger L
  !
  ! n_real_eigv is the number of real eigenvalues, which might be smaller than the dimension "n" 
  !
  END_DOC

  implicit none
  integer,          intent(in)  :: n
  double precision, intent(in)  :: A(n,n),B(n,n)
  integer,          intent(out) :: n_real_eigv
  double precision, intent(out) :: reigvec(n,n), leigvec(n,n), eigval(n)

  integer                       :: i, j
  integer                       :: n_bad
  double precision              :: thr
  double precision              :: accu_nd

  integer,          allocatable :: iorder(:)
  double precision, allocatable :: Aw(:,:),Bw(:,:)
  double precision, allocatable :: WR(:), WI(:), VL(:,:), VR(:,:), beta(:)
  double precision, allocatable :: S(:,:)
  double precision :: r

  ! -------------------------------------------------------------------------------------
  !

  print *, 'Computing the left/right eigenvectors ...'
  allocate( WR(n), WI(n), VL(n,n), VR(n,n), Aw(n,n), Bw(n,n),iorder(n),beta(n))

  Aw(:,:) = A(:,:)
  Bw(:,:) = B(:,:)
  call lapack_diag_general_non_sym(n,Aw,Bw,WR,beta,WI,VL,VR)

  ! -------------------------------------------------------------------------------------
  !                  track & sort the real eigenvalues 

  i = 1
  thr    = 1.d-10
  n_real_eigv = 0
  do while (i.le.n) 
    if( dabs(WI(i)).gt.thr ) then
      print*, 'Found an imaginary component to eigenvalue on i = ', i
      print*, 'Re(i) , Im(i)  ', WR(i), WI(i)
      iorder(i) = i
      eigval(i) = WR(i)/(beta(i) + 1.d-10)
      i+=1
      print*, 'Re(i+1),Im(i+1)',WR(i), WI(i)
      iorder(i) = i
      eigval(i) = WR(i)/(beta(i) + 1.d-10)
      i+=1
    else  
      n_real_eigv += 1
      iorder(i) = i
      eigval(i) = WR(i)/(beta(i) + 1.d-10)
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

  !print *, ' check_EIGVEC :'
  !call check_EIGVEC(n, n_real_eigv, A, eigval, leigvec, reigvec)

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

  deallocate( S )

end 

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

subroutine non_hrmt_bieig_fullvect(n, A, leigvec, reigvec, n_real_eigv, eigval)

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

  integer,          allocatable :: iorder(:)
  double precision, allocatable :: Aw(:,:)
  double precision, allocatable :: WR(:), WI(:), VL(:,:), VR(:,:)
  double precision, allocatable :: S(:,:)
  double precision, allocatable :: eigval_sorted(:)


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

  allocate( eigval_sorted(n), iorder(n) )

  n_good = 0
  thr    = 1.d-10

  do i = 1, n

    iorder(i) = i
    eigval_sorted(i) = WR(i)

    if(dabs(WI(i)) .gt. thr) then
      print*, ' Found an imaginary component to eigenvalue on i = ', i
      print*, ' Re(i) + Im(i)', WR(i), WI(i)
    else
      n_good += 1
    endif

  enddo

  n_real_eigv = n_good 

  call dsort(eigval_sorted, iorder, n)
      
  reigvec(:,:) = 0.d0 
  leigvec(:,:) = 0.d0 
  do i = 1, n
    eigval(i) = WR(i)
    do j = 1, n
      reigvec(j,i) = VR(j,iorder(i))
      leigvec(j,i) = VL(j,iorder(i))
    enddo
  enddo

  deallocate( eigval_sorted, iorder )
  deallocate( WR, WI )
  deallocate( VL, VR )

  !
  ! -------------------------------------------------------------------------------------

  ! ---

  ! -------------------------------------------------------------------------------------
  !                               check bi-orthogonality

  allocate( S(n,n) )

  ! S = VL x VR
  call dgemm( 'T', 'N', n, n, n, 1.d0                              &
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

    !print *, ' L & T bi-orthogonality: ok'
    deallocate( S )
    return

  else
    ! impose bi-orthogonality 

    !print *, ' L & T bi-orthogonality: not imposed yet'
    !print *, ' accu_nd = ', accu_nd
    call impose_biorthog_qr(n, n, leigvec, reigvec, S)
    deallocate( S )
  
  endif

  !
  ! -------------------------------------------------------------------------------------

  return

end subroutine non_hrmt_bieig_fullvect

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

      if(de .gt. de_thr) then

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

