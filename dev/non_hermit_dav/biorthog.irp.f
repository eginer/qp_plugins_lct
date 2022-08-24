
subroutine non_hrmt_diag_split_degen(n, A, leigvec, reigvec, n_real_eigv, eigval)

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
  double precision              :: reigvec_tmp(n,n), leigvec_tmp(n,n)

  integer                       :: i, j, n_degen,k , iteration
  integer                       :: n_good
  double precision              :: shift,shift_current
  double precision              :: r,thr
  integer,          allocatable :: list_good(:), iorder_origin(:),iorder(:)
  double precision, allocatable :: WR(:), WI(:), Vl(:,:), VR(:,:),S(:,:)
  double precision, allocatable :: Aw(:,:),diag_elem(:),A_save(:,:)
  double precision, allocatable :: im_part(:),re_part(:)


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
  do i = 1, n
   do j = 1, n
    reigvec(iorder_origin(j),i) = reigvec_tmp(j,i)
    leigvec(iorder_origin(j),i) = leigvec_tmp(j,i)
   enddo
  enddo

end subroutine non_hrmt_diag_split_degen


! ---

subroutine non_hrmt_real_diag_new(n, A, leigvec, reigvec, n_real_eigv, eigval)

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

subroutine non_hrmt_real_diag(n, A, leigvec, reigvec, n_real_eigv, eigval)

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

  integer                       :: i, j, n_good
  double precision              :: thr, threshold, accu_d, accu_nd
  integer,          allocatable :: list_good(:), iorder(:)
  double precision, allocatable :: Aw(:,:)
  double precision, allocatable :: WR(:), WI(:), Vl(:,:), VR(:,:), S(:,:), S_inv_half_tmp(:,:)

  print*, ' Computing the left/right eigenvectors ...'

  ! Eigvalue(n) = WR(n) + i * WI(n)
  allocate(WR(n), WI(n), VL(n,n), VR(n,n), Aw(n,n))
  Aw = A
  call lapack_diag_non_sym(n, Aw, WR, WI, VL, VR)

  ! ---
  ! You track the real eigenvalues 

  thr = 1d-15

  n_good = 0
  do i = 1, n
    if(dabs(WI(i)).lt.thr) then
      n_good += 1
    else
      print*, ' Found an imaginary component to eigenvalue'
      print*, ' Re(i) + Im(i)', WR(i), WI(i)
    endif
  enddo

  allocate(list_good(n_good), iorder(n_good))
  n_good = 0
  do i = 1, n
    if(dabs(WI(i)).lt.thr) then
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
  call dsort(eigval, iorder, n_good)
  do i = 1, n_real_eigv
    do j = 1, n
      reigvec(j,i) = VR(j,list_good(iorder(i)))
      leigvec(j,i) = Vl(j,list_good(iorder(i)))
    enddo
  enddo

  ! ---

  allocate( S(n_real_eigv,n_real_eigv), S_inv_half_tmp(n_real_eigv,n_real_eigv) )

  ! S = VL x VR
  call dgemm( 'T', 'N', n_real_eigv, n_real_eigv, n_real_eigv, 1.d0 &
            , leigvec, size(leigvec, 1), reigvec, size(reigvec, 1)  &
            , 0.d0, S, size(S, 1) )

  accu_nd = 0.d0
  accu_d  = 0.d0
  do i = 1, n_real_eigv
    do j = 1, n_real_eigv
      if(i==j) then
        accu_d += S(j,i)
      else 
        accu_nd = accu_nd + S(j,i) * S(j,i)
      endif
    enddo
  enddo
  accu_nd = dsqrt(accu_nd)

  threshold = 1.d-15 
  if( (accu_nd .gt. threshold) .or. (dabs(accu_d-dble(n_real_eigv)) .gt. threshold) ) then

    print*, ' sum of off-diag S elements = ', accu_nd
    print*, ' Should be zero '
    print*, ' sum of     diag S elements = ', accu_d
    print*, ' Should be ',n
    print*, ' Not bi-orthonormal !!'
    print*, ' Notice that if you are interested in ground state it is not a problem :)'
  endif

end subroutine non_hrmt_real_diag

! ---

subroutine lapack_diag_general_non_sym(n,A,B,WR,beta,WI,VL,VR)

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

 implicit none
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

! ---

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
  double precision, allocatable :: WR(:), WI(:), VL(:,:), VR(:,:)
  double precision, allocatable :: S(:,:)


  ! -------------------------------------------------------------------------------------
  !

  print *, ' Computing the left/right eigenvectors ...'

  allocate( WR(n), WI(n), VL(n,n), VR(n,n) ) 
  
  print *, ' fock matrix'
  do i = 1, n
    write(*, '(1000(F16.10,X))') A(i,:)
  enddo

  !call lapack_diag_non_sym_right(n, A, WR, WI, VR)
  call lapack_diag_non_sym(n, A, WR, WI, VL, VR)
  !call lapack_diag_non_sym_new(n, A, WR, WI, VL, VR)

  print *, ' eigenvalues'
  do i = 1, n
    write(*, '(1000(F16.10,X))') WR(i), WI(i)
  enddo
  print *, ' right eigenvect bef' 
  do i = 1, n
    write(*, '(1000(F16.10,X))') VR(:,i)
  enddo
  print *, ' left eigenvect bef'
  do i = 1, n
    write(*, '(1000(F16.10,X))') VL(:,i)
  enddo

  !call rotate_degen_eigvec(n, VR)
  !call rotate_degen_eigvec(n, VL)
  !print *, ' right eigenvect aft' 
  !do i = 1, n
  !  write(*, '(1000(F16.10,X))') VR(:,i)
  !enddo
  !print *, ' left eigenvect aft'
  !do i = 1, n
  !  write(*, '(1000(F16.10,X))') VL(:,i)
  !enddo

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

  print *, ' eigenvalues'
  do i = 1, n
    write(*, '(1000(F16.10,X))') eigval(i)
  enddo
  print *, ' right eigenvect aft ord' 
  do i = 1, n
    write(*, '(1000(F16.10,X))') reigvec(:,i)
  enddo
  print *, ' left eigenvect aft ord'
  do i = 1, n
    write(*, '(1000(F16.10,X))') leigvec(:,i)
  enddo

  !print *, ' check_EIGVEC before QR:'
  !call check_EIGVEC(n, n_real_eigv, A, eigval, leigvec, reigvec)

  !
  ! -------------------------------------------------------------------------------------

  ! ---

  call check_degen(n, n_real_eigv, eigval, leigvec, reigvec)


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
    call impose_biorthog_qr(n, n_real_eigv, leigvec, reigvec, S)
    !call impose_biorthog_lu(n, n_real_eigv, leigvec, reigvec, S)
    deallocate( S )

    print *, ' right eigenvect aft QR' 
    do i = 1, n
      write(*, '(1000(F16.10,X))') reigvec(:,i)
    enddo
    print *, ' left eigenvect aft QR'
    do i = 1, n
      write(*, '(1000(F16.10,X))') leigvec(:,i)
    enddo

    !call check_EIGVEC(n, n_real_eigv, A, eigval, leigvec, reigvec)

    stop
  
  endif

  !
  ! -------------------------------------------------------------------------------------

  return

end subroutine non_hrmt_bieig

! ---

subroutine non_hrmt_bieiginv(n, A, leigvec, reigvec, n_real_eigv, eigval)

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

  integer                       :: i, j, n_good
  double precision              :: thr, accu_nd, r

  integer,          allocatable :: list_good(:), iorder(:)
  double precision, allocatable :: Aw(:,:)
  double precision, allocatable :: WR(:), WI(:), VL(:,:), VR(:,:)
  double precision, allocatable :: S(:,:)


  ! -------------------------------------------------------------------------------------
  !

  print *, 'Computing the left/right eigenvectors ...'

  allocate( WR(n), WI(n), VL(n,n), VR(n,n), Aw(n,n) )
  Aw(:,:) = A(:,:)
  do i = 1, n
    call RANDOM_NUMBER(r)
    Aw(i,i) += 1d-10 * r
  enddo

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

  ! S = VL x VR
  allocate( S(n_real_eigv,n_real_eigv) )
  call dgemm( 'T', 'N', n_real_eigv, n_real_eigv, n, 1.d0          &
            , leigvec, size(leigvec, 1), reigvec, size(reigvec, 1) &
            , 0.d0, S, size(S, 1) )
  accu_nd = 0.d0
  do i = 1, n_real_eigv
    do j = 1, n_real_eigv
      if(i == j) cycle
      accu_nd = accu_nd + S(j,i) * S(j,i)
    enddo
  enddo
  deallocate( S )

  accu_nd = dsqrt(accu_nd)

  if(accu_nd .lt. 1d-8) then

    print *, ' L & T bi-orthogonality: ok'
    print *, ' accu_nd = ', accu_nd
    return

  else

    print *, ' L & T bi-orthogonality: not imposed yet'
    print *, ' accu_nd = ', accu_nd
    call impose_biorthog_inv(n, n_real_eigv, leigvec, reigvec)

    ! S = VL x VR
    allocate( S(n_real_eigv,n_real_eigv) )
    call dgemm( 'T', 'N', n_real_eigv, n_real_eigv, n, 1.d0          &
              , leigvec, size(leigvec, 1), reigvec, size(reigvec, 1) &
              , 0.d0, S, size(S, 1) )
    accu_nd = 0.d0
    do i = 1, n_real_eigv
      do j = 1, n_real_eigv
        if(i == j) cycle
        accu_nd = accu_nd + S(j,i) * S(j,i)
      enddo
    enddo
    deallocate( S )
    accu_nd = dsqrt(accu_nd)
    if(accu_nd .lt. 1d-8) then
      print *, ' L & T bi-orthogonality: ok'
      print *, ' accu_nd = ', accu_nd
      stop
    endif
    print *, ' accu_nd = ', accu_nd

    !print *, ' check_EIGVEC after QR:'
    !call check_EIGVEC(n, n_real_eigv, A, eigval, leigvec, reigvec)
  
  endif

  !
  ! -------------------------------------------------------------------------------------

  return

end subroutine non_hrmt_bieiginv

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

  if(accu_nd .lt. 1d-8) then
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

end subroutine non_hrmt_generalized_real_im

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


subroutine split_matrix_degen(aw,n,shift)
 implicit none
 BEGIN_DOC
 ! subroutines that splits the degeneracies of a matrix by adding a splitting of magnitude thr * n_degen/2
 !
 ! WARNING !! THE MATRIX IS ASSUMED TO BE PASSED WITH INCREASING DIAGONAL ELEMENTS
 END_DOC
 double precision,intent(inout) :: Aw(n,n)
 double precision,intent(in)    :: shift
 integer, intent(in) :: n
 integer :: i,j,n_degen
 i=1
 do while(i.lt.n)
  if(dabs(Aw(i,i)-Aw(i+1,i+1)).lt.shift)then
   j=1
   do while(dabs(Aw(i,i)-Aw(i+j,i+j)).lt.shift.and.i+j.le.n)
    j+=1
   enddo
   n_degen = j
   j=0
   do while(dabs(Aw(i+j,i+j)-Aw(i+j+1,i+j+1)).lt.shift.and.i+j+1.lt.n)
    Aw(i+j,i+j) += (j-n_degen/2) * shift
    j+=1
   enddo
   Aw(i+n_degen-1,i+n_degen-1) += (n_degen-1-n_degen/2) * shift
   i+=n_degen
  else 
   i+=1
  endif
 enddo

end

subroutine check_bi_ortho(reigvec,leigvec,n,S,accu_nd)
 implicit none
 integer, intent(in) :: n
 double precision,intent(in) :: reigvec(n,n),leigvec(n,n)
 double precision, intent(out) :: S(n,n),accu_nd
 BEGIN_DOC
! retunrs the overlap matrix S = Leigvec^T Reigvec 
!
! and the square root of the sum of the squared off-diagonal elements of S
 END_DOC
 integer :: i,j
  ! S = VL x VR
  call dgemm( 'T', 'N', n, n, n, 1.d0 &
            , leigvec, size(leigvec, 1), reigvec, size(reigvec, 1)  &
            , 0.d0, S, size(S, 1) )
  accu_nd = 0.d0
  do i = 1, n
    do j = 1, n
      if(i.ne.j) then
        accu_nd = accu_nd + S(j,i) * S(j,i)
      endif
    enddo
  enddo
  accu_nd = dsqrt(accu_nd)

end
