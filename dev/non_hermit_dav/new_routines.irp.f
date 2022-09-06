subroutine non_hrmt_diag_split_degen_bi_orthog(n, A, leigvec, reigvec, n_real_eigv, eigval)

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
  double precision, allocatable :: reigvec_tmp(:,:), leigvec_tmp(:,:)

  integer                       :: i, j, n_degen,k , iteration
  double precision              :: shift_current
  double precision              :: r,thr,accu_d, accu_nd
  integer,          allocatable :: iorder_origin(:),iorder(:)
  double precision, allocatable :: WR(:), WI(:), Vl(:,:), VR(:,:),S(:,:)
  double precision, allocatable :: Aw(:,:),diag_elem(:),A_save(:,:)
  double precision, allocatable :: im_part(:),re_part(:)
  double precision :: accu,thr_cut


  thr_cut = 1.d-15
  print*,'Computing the left/right eigenvectors ...'
  print*,'Using the degeneracy splitting algorithm'
 ! initialization
  shift_current = 1.d-15
  iteration = 0 
  print*,'***** iteration = ',iteration


  ! pre-processing the matrix :: sorting by diagonal elements
  allocate(reigvec_tmp(n,n), leigvec_tmp(n,n))
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

  allocate(WR(n), WI(n), VL(n,n), VR(n,n), Aw(n,n))
  allocate(im_part(n),iorder(n))
  allocate( S(n,n) )


  Aw = A_save
  call cancel_small_elmts(aw,n,thr_cut)
  call lapack_diag_non_sym(n,Aw,WR,WI,VL,VR)
  do i = 1, n
   im_part(i) = -dabs(WI(i))
   iorder(i) = i
  enddo
  call dsort(im_part, iorder, n)
  n_real_eigv = 0
  do i = 1, n
    if(dabs(WI(i)).lt.1.d-20)then
      n_real_eigv += 1
    else
!      print*,'Found an imaginary component to eigenvalue'
!      print*,'Re(i) + Im(i)',WR(i),WI(i)
    endif
  enddo
  if(n_real_eigv.ne.n)then
   shift_current = max(10.d0 * dabs(im_part(1)),shift_current*10.d0)
   print*,'Largest imaginary part found in eigenvalues = ',im_part(1)
   print*,'Splitting the degeneracies by ',shift_current
  else
   print*,'All eigenvalues are real !'
  endif


  do while(n_real_eigv.ne.n)
   iteration += 1
   print*,'***** iteration = ',iteration
   if(shift_current.gt.1.d-3)then
    print*,'shift_current > 1.d-3 !!'
    print*,'Your matrix intrinsically contains complex eigenvalues'
    stop
   endif
   Aw = A_save
   call cancel_small_elmts(Aw,n,thr_cut)
   call split_matrix_degen(Aw,n,shift_current)
   call lapack_diag_non_sym(n,Aw,WR,WI,VL,VR)
   n_real_eigv = 0
   do i = 1, n
     if(dabs(WI(i)).lt.1.d-20)then
       n_real_eigv+= 1
     else
!       print*,'Found an imaginary component to eigenvalue'
!       print*,'Re(i) + Im(i)',WR(i),WI(i)
     endif
   enddo
   if(n_real_eigv.ne.n)then
    do i = 1, n
     im_part(i) = -dabs(WI(i))
     iorder(i) = i
    enddo
    call dsort(im_part, iorder, n)
    shift_current = max(10.d0 * dabs(im_part(1)),shift_current*10.d0)
    print*,'Largest imaginary part found in eigenvalues = ',im_part(1)
    print*,'Splitting the degeneracies by ',shift_current
   else
    print*,'All eigenvalues are real !'
   endif

  enddo

!!! ONCE ALL EIGENVALUES ARE REAL ::: CHECK BI-ORTHONORMALITY
  n_real_eigv = n
  reigvec_tmp(:,:) = 0.d0 
  leigvec_tmp(:,:) = 0.d0 
  do i = 1, n
    eigval(i) = WR(i)
    do j = 1, n
      reigvec_tmp(j,i) = VR(j,i)
      leigvec_tmp(j,i) = Vl(j,i)
    enddo
!    print*,'WR(i) = ',WR(i)
  enddo
!  print*,'Checking the eigenvectors ...'
!  print*,'Thr for eigenvectors = ',shift_current
!  call check_EIGVEC(n, n, Aw, eigval, leigvec_tmp, reigvec_tmp,shift_current)
  ! -------------------------------------------------------------------------------------
  !                               check bi-orthogonality
  call check_biorthog(n, n_real_eigv, leigvec_tmp, reigvec_tmp, accu_d, accu_nd, S)
  print *, ' accu_nd bi-orthog = ', accu_nd
  
  if( accu_nd .lt. 1d-10 ) then
  
    print *, ' bi-orthogonality: ok'
  
  else
  
    print *, ' '
    print *, ' bi-orthogonality: not imposed yet'
    print *, ' '
  
    ! ---
  
    print *, ' '
    print *, ' orthog between degen eigenvect' 
    print *, ' '
  
    call impose_orthog_degen_eigvec(n, eigval, reigvec_tmp)
!!   print *, ' right eigenvect aft orthog' 
!!   do i = 1, n
!!     write(*, '(1000(F16.10,X))') reigvec_tmp(:,i)
!!   enddo
  
    call impose_orthog_degen_eigvec(n, eigval, leigvec_tmp)
!!   print *, ' left eigenvect aft orthog' 
!!   do i = 1, n
!!     write(*, '(1000(F16.10,X))') leigvec_tmp(:,i)
!!   enddo
  
  
    call check_biorthog(n, n, leigvec_tmp, reigvec_tmp, accu_d, accu_nd, S)
    if( accu_nd .lt. 1d-10 ) then
      print *, ' bi-orthogonality: ok'
    endif
  endif

!  print*,'Rechecking the eigenvectors ...'
!  print*,'Thr for eigenvectors = ',shift_current
!  call check_EIGVEC(n, n, A_save, eigval, leigvec_tmp, reigvec_tmp,shift_current)
 
  do i = 1, n
   do j = 1, n
    VR(iorder_origin(j),i) = reigvec_tmp(j,i)
    VL(iorder_origin(j),i) = leigvec_tmp(j,i)
   enddo
  enddo

  eigval = 0.d0
  do i = 1, n
   iorder(i) = i
   accu = 0.d0
   do j = 1, n
    accu += VL(j,i) * VR(j,i) 
    do k = 1, n
     eigval(i) +=  VL(j,i) * A(j,k) * VR(k,i) 
    enddo
   enddo
   eigval(i) *= 1.d0/accu
!   print*,'eigval(i) = ',eigval(i)
  enddo
  ! sorting the 
  call dsort(eigval, iorder, n)
!  print*,'Rechecking the eigenvectors again ...'
!  print*,'Thr for eigenvectors = ',shift_current
!  call check_EIGVEC(n, n, A, eigval, VL, VR,shift_current)
!  call check_biorthog(n, n, VL, VR, accu_d, accu_nd, S)
!  print*,'accu_nd = ',accu_nd
  do i = 1, n
   do j = 1, n
    reigvec(j,i) = VR(j,iorder(i))
    leigvec(j,i) = VL(j,iorder(i))
   enddo
  enddo
  print*,'Checking for final reigvec/leigvec'
  print*,'Thr for eigenvectors = ',shift_current
  call check_EIGVEC(n, n, A, eigval, leigvec, reigvec,shift_current)
  call check_biorthog(n, n, leigvec, reigvec, accu_d, accu_nd, S)
  print *, ' accu_nd bi-orthog = ', accu_nd
  
  if( accu_nd .lt. 1d-10 ) then
    print *, ' bi-orthogonality: ok'
  else 
   print*,'Something went wrong in non_hrmt_diag_split_degen_bi_orthog'
   print*,'Eigenvectors are not bi orthonormal ..'
   print*,'accu_nd = ',accu_nd
   stop
  endif

end subroutine non_hrmt_diag_split_degen


