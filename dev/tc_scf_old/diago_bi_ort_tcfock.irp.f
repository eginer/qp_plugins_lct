 BEGIN_PROVIDER [ double precision, fock_tc_reigvec_mo, (mo_num, mo_num)]
&BEGIN_PROVIDER [ double precision, fock_tc_leigvec_mo, (mo_num, mo_num)]
&BEGIN_PROVIDER [ double precision, eigval_fock_tc_mo, (mo_num)]
&BEGIN_PROVIDER [ double precision, overlap_fock_tc_eigvec_mo, (mo_num, mo_num)]

  BEGIN_DOC
  ! EIGENVECTORS OF FOCK MATRIX ON THE MO BASIS and their OVERLAP
  END_DOC

  implicit none
  integer                       :: n_real_tc 
  integer                       :: i, k, l
  double precision              :: accu_d, accu_nd, accu_tmp
  double precision              :: norm
  double precision, allocatable :: eigval_right_tmp(:)

  allocate( eigval_right_tmp(mo_num) )

  PROVIDE Fock_matrix_tc_mo_tot

!  print*,'Fock matrix'
!  do i = 1, mo_num
!   write(*,'(100(F16.10,X))')Fock_matrix_tc_mo_tot(:,i)
!  enddo
!   call non_hrmt_diag_split_degen( mo_num, Fock_matrix_tc_mo_tot &

  call non_hrmt_bieig( mo_num, Fock_matrix_tc_mo_tot          &
                     , fock_tc_leigvec_mo, fock_tc_reigvec_mo & 
                     , n_real_tc, eigval_right_tmp )
  !if(max_ov_tc_scf)then
  ! call non_hrmt_fock_mat( mo_num, Fock_matrix_tc_mo_tot &
  !                    , fock_tc_leigvec_mo, fock_tc_reigvec_mo & 
  !                    , n_real_tc, eigval_right_tmp )
  !else 
  ! call non_hrmt_diag_split_degen_bi_orthog( mo_num, Fock_matrix_tc_mo_tot &
  !                    , fock_tc_leigvec_mo, fock_tc_reigvec_mo & 
  !                    , n_real_tc, eigval_right_tmp )
  !endif

!  if(n_real_tc .ne. mo_num)then
!   print*,'n_real_tc ne mo_num ! ',n_real_tc
!   stop
!  endif

  eigval_fock_tc_mo = eigval_right_tmp
!  print*,'Eigenvalues of Fock_matrix_tc_mo_tot'
!  do i = 1, mo_num
!    print*, i, eigval_fock_tc_mo(i)
!  enddo
!  deallocate( eigval_right_tmp )

  ! L.T x R 
  call dgemm( "T", "N", mo_num, mo_num, mo_num, 1.d0          &
            , fock_tc_leigvec_mo, size(fock_tc_leigvec_mo, 1) &
            , fock_tc_reigvec_mo, size(fock_tc_reigvec_mo, 1) &
            , 0.d0, overlap_fock_tc_eigvec_mo, size(overlap_fock_tc_eigvec_mo, 1) )

  accu_d  = 0.d0
  accu_nd = 0.d0
  do i = 1, mo_num
    do k = 1, mo_num
      if(i==k) then
        accu_tmp = overlap_fock_tc_eigvec_mo(k,i)
        accu_d += accu_tmp 
      else
        accu_tmp = overlap_fock_tc_eigvec_mo(k,i)
        accu_nd += accu_tmp * accu_tmp
        if(dabs(overlap_fock_tc_eigvec_mo(k,i)).gt.1.d-10)then
         print*,'k,i',k,i,overlap_fock_tc_eigvec_mo(k,i)
        endif
      endif
    enddo 
  enddo
  accu_nd = dsqrt(accu_nd)

  if( accu_nd .gt. 1d-8 ) then
    print *, ' bi-orthog failed'
    print*,'accu_nd MO = ', accu_nd
    print*,'overlap_fock_tc_eigvec_mo = '
    do i = 1, mo_num
      write(*,'(100(F16.10,X))') overlap_fock_tc_eigvec_mo(i,:)
    enddo
   stop
  endif

  if( dabs(accu_d - dble(mo_num)) .gt. 1e-7 ) then
    print *, 'mo_num     = ', mo_num 
    print *, 'accu_d  MO = ', accu_d
    print *, 'normalizing vectors ...'
    do i = 1, mo_num
      norm = dsqrt(dabs(overlap_fock_tc_eigvec_mo(i,i)))
      if( norm.gt.1e-7 ) then
        do k = 1, mo_num
          fock_tc_reigvec_mo(k,i) *= 1.d0/norm
          fock_tc_leigvec_mo(k,i) *= 1.d0/norm
        enddo
      endif
    enddo
    call dgemm( "T", "N", mo_num, mo_num, mo_num, 1.d0          &
              , fock_tc_leigvec_mo, size(fock_tc_leigvec_mo, 1) &
              , fock_tc_reigvec_mo, size(fock_tc_reigvec_mo, 1) &
              , 0.d0, overlap_fock_tc_eigvec_mo, size(overlap_fock_tc_eigvec_mo, 1) )
  endif
 
END_PROVIDER 

 BEGIN_PROVIDER [ double precision, fock_tc_reigvec_ao, (ao_num, mo_num)]
&BEGIN_PROVIDER [ double precision, fock_tc_leigvec_ao, (ao_num, mo_num)]
&BEGIN_PROVIDER [ double precision, overlap_fock_tc_eigvec_ao, (mo_num, mo_num) ]

  BEGIN_DOC
  ! EIGENVECTORS OF FOCK MATRIX ON THE AO BASIS and their OVERLAP
  !
  ! THE OVERLAP SHOULD BE THE SAME AS overlap_fock_tc_eigvec_mo
  END_DOC

  implicit none
  integer                       :: i, j, k, q, p
  double precision              :: accu, accu_d
  double precision, allocatable :: tmp(:,:)


!  ! MO_R x R
   call dgemm( 'N', 'N', ao_num, mo_num, mo_num, 1.d0          &
             , mo_r_coef, size(mo_r_coef, 1)                   &
             , fock_tc_reigvec_mo, size(fock_tc_reigvec_mo, 1) &
             , 0.d0, fock_tc_reigvec_ao, size(fock_tc_reigvec_ao, 1) )

   ! MO_L x L
   call dgemm( 'N', 'N', ao_num, mo_num, mo_num, 1.d0          &
             , mo_l_coef, size(mo_l_coef, 1)                   &
             , fock_tc_leigvec_mo, size(fock_tc_leigvec_mo, 1) &
             , 0.d0, fock_tc_leigvec_ao, size(fock_tc_leigvec_ao, 1) )
 
  allocate( tmp(mo_num,ao_num) )

  ! tmp <-- L.T x S_ao
  call dgemm( "T", "N", mo_num, ao_num, ao_num, 1.d0                                           &
            , fock_tc_leigvec_ao, size(fock_tc_leigvec_ao, 1), ao_overlap, size(ao_overlap, 1) &
            , 0.d0, tmp, size(tmp, 1) )

  ! S <-- tmp x R
  call dgemm( "N", "N", mo_num, mo_num, ao_num, 1.d0                             &
            , tmp, size(tmp, 1), fock_tc_reigvec_ao, size(fock_tc_reigvec_ao, 1) &
            , 0.d0, overlap_fock_tc_eigvec_ao, size(overlap_fock_tc_eigvec_ao, 1) )

  deallocate( tmp )

  ! ---
  double precision :: norm
  do i = 1, mo_num
   norm = 1.d0/dsqrt(dabs(overlap_fock_tc_eigvec_ao(i,i)))
   do j = 1, mo_num
    fock_tc_reigvec_ao(j,i) *= norm
    fock_tc_leigvec_ao(j,i) *= norm
   enddo
  enddo

  allocate( tmp(mo_num,ao_num) )

  ! tmp <-- L.T x S_ao
  call dgemm( "T", "N", mo_num, ao_num, ao_num, 1.d0                                           &
            , fock_tc_leigvec_ao, size(fock_tc_leigvec_ao, 1), ao_overlap, size(ao_overlap, 1) &
            , 0.d0, tmp, size(tmp, 1) )

  ! S <-- tmp x R
  call dgemm( "N", "N", mo_num, mo_num, ao_num, 1.d0                             &
            , tmp, size(tmp, 1), fock_tc_reigvec_ao, size(fock_tc_reigvec_ao, 1) &
            , 0.d0, overlap_fock_tc_eigvec_ao, size(overlap_fock_tc_eigvec_ao, 1) )

  deallocate( tmp )

END_PROVIDER

