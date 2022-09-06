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
  integer                       :: n_good
  double precision              :: shift,shift_current
  double precision              :: r,thr,accu_d, accu_nd
  integer,          allocatable :: list_good(:), iorder_origin(:),iorder(:)
  double precision, allocatable :: WR(:), WI(:), Vl(:,:), VR(:,:),S(:,:)
  double precision, allocatable :: Aw(:,:),diag_elem(:),A_save(:,:)
  double precision, allocatable :: im_part(:),re_part(:)
  double precision :: accu


  print*,'Computing the left/right eigenvectors ...'
  print*,'Using the degeneracy splitting algorithm'


  ! pre-processing the matrix :: sorting by diagonal elements
  allocate(reigvec_tmp(n,n), leigvec_tmp(n,n))
  allocate(diag_elem(n),iorder_origin(n),A_save(n,n))
  do i = 1, n
   iorder_origin(i) = i
   diag_elem(i) = A(i,i)
   print*,'diag_elem(i) = ',i,diag_elem(i)
  enddo
  call dsort(diag_elem, iorder_origin, n)
  do i = 1, n
  print*,i,iorder_origin(i),diag_elem(i)
!   iorder_origin(i) = i
   do j = 1, n
    A_save(j,i) = A(iorder_origin(j),iorder_origin(i))
   enddo
  enddo

  allocate(WR(n), WI(n), VL(n,n), VR(n,n), Aw(n,n))
  allocate(im_part(n),iorder(n))
  allocate( S(n,n) )
  allocate( list_good(n))

  shift = 1.d-15
  shift_current = shift
  iteration = 1 
  do while(n_real_eigv.ne.n)
   if(shift.gt.1.d-3)then
    print*,'shift > 1.d-3 !!'
    print*,'Your matrix intrinsically contains complex eigenvalues'
    stop
   endif
   print*,'***** iteration = ',iteration
   print*,'shift = ',shift
   Aw = A_save
   do i = 1, n
    do j = 1, n
     if(dabs(Aw(j,i)).lt.shift)then
      Aw(j,i) = 0.d0
     endif
    enddo
   enddo
   call lapack_diag_non_sym(n,Aw,WR,WI,VL,VR)
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
   call lapack_diag_non_sym(n,Aw,WR,WI,VL,VR)
   ! You track the real eigenvalues 
   n_good = 0
   do i = 1, n
     print*,'WR(i) = ',WR(i)
     if(dabs(WI(i)).lt.1.d-20)then
       n_good += 1
     else
       print*,'Found an imaginary component to eigenvalue'
       print*,'Re(i) + Im(i)',WR(i),WI(i)
     endif
   enddo
   call check_EIGVEC(n, n, Aw, WR, VL, VR)

   if(n_good == n)then
!!!!! ONCE ALL EIGENVALUES ARE REAL ::: CHECK BI-ORTHONORMALITY
    print*,'All eigenvalues are real !'
    n_real_eigv = n
    reigvec_tmp(:,:) = 0.d0 
    leigvec_tmp(:,:) = 0.d0 
    do i = 1, n
      eigval(i) = WR(i)
      do j = 1, n
        reigvec_tmp(j,i) = VR(j,i)
        leigvec_tmp(j,i) = Vl(j,i)
      enddo
      print*,'WR(i) = ',WR(i)
    enddo
    call check_EIGVEC(n, n, Aw, eigval, leigvec_tmp, reigvec_tmp)
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
      print *, ' right eigenvect aft orthog' 
      do i = 1, n
        write(*, '(1000(F16.10,X))') reigvec_tmp(:,i)
      enddo
  
      call impose_orthog_degen_eigvec(n, eigval, leigvec_tmp)
      print *, ' left eigenvect aft orthog' 
      do i = 1, n
        write(*, '(1000(F16.10,X))') leigvec_tmp(:,i)
      enddo
  
  
      call check_biorthog(n, n, leigvec_tmp, reigvec_tmp, accu_d, accu_nd, S)
      if( accu_nd .lt. 1d-10 ) then
        print *, ' bi-orthogonality: ok'
      endif
    endif
   endif
   shift *= 10.d0
   iteration += 1
  enddo
  deallocate( im_part, list_good, iorder )
  deallocate(WR, WI, VL, VR, Aw)
  deallocate( S )

  allocate(S(n,n),WR(n),iorder(n),VR(n,n),VL(n,n))
  print*,'Final check '
  print*,'First with A_save'
  call check_EIGVEC(n, n, A_save, eigval, leigvec_tmp, reigvec_tmp)
  WR = 0.d0
  double precision :: tmp1,tmp2
 
  do i = 1, n
   do j = 1, n
    VR(iorder_origin(j),i) = reigvec_tmp(j,i)
    VL(iorder_origin(j),i) = leigvec_tmp(j,i)
   enddo
  enddo

  do i = 1, n
   iorder(i) = i
   tmp1 = 0.d0
   tmp2 = 0.d0
   accu = 0.d0
   do j = 1, n
    accu += VL(j,i) * VR(j,i)
    tmp1 += reigvec_tmp(j,i) * leigvec_tmp(j,i)
    do k = 1, n
     WR(i) +=  VL(j,i) * A(j,k) * VR(k,i) 
     tmp2 += leigvec_tmp(j,i) * A_save(j,k) * reigvec_tmp(k,i)
    enddo
   enddo
   tmp2 *= 1.d0/tmp1
   WR(i) *= 1.d0/accu
   print*,'WR(i)/eigval/error ',WR(i),eigval(i),dabs(WR(i)-eigval(i))
   print*,'tmp  /eigval/error ',tmp2 ,eigval(i),dabs(tmp2 -eigval(i))
  enddo
  ! sorting the 
  call check_EIGVEC(n, n, A, WR, VL, VR)
  call check_biorthog(n, n, VL, VR, accu_d, accu_nd, S)
  print*,'accu_nd = ',accu_nd
  call dsort(WR, iorder, n)
  do i = 1, n
   print*,'WR(i) = ',WR(i)
   do j = 1, n
    reigvec(j,i) = VR(j,iorder(i))
    leigvec(j,i) = VL(j,iorder(i))
   enddo
  enddo
  print*,'Checking for final reigvec/leigvec'
  call check_EIGVEC(n, n, A, WR, leigvec, reigvec)
  call check_biorthog(n, n, leigvec, reigvec, accu_d, accu_nd, S)
  print*,'accu_nd = ',accu_nd
  
     if( accu_nd .lt. 1d-10 ) then
       print *, ' bi-orthogonality: ok'
     endif
 stop

end subroutine non_hrmt_diag_split_degen


