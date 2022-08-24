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

! call non_hrmt_generalized_real_im(mo_num, Fock_matrix_tc_mo_tot, overlap_bi_ortho &
   !call non_hrmt_real_im( mo_num, Fock_matrix_tc_mo_tot &
   !call non_hrmt_real_diag( mo_num, Fock_matrix_tc_mo_tot &
!   call non_hrmt_real_diag_new( mo_num, Fock_matrix_tc_mo_tot &
!   call non_hrmt_bieig( mo_num, Fock_matrix_tc_mo_tot &
!   call non_hrmt_bieig_real_im( mo_num, Fock_matrix_tc_mo_tot &
!   call non_hrmt_bieiginv( mo_num, Fock_matrix_tc_mo_tot &
   call non_hrmt_diag_split_degen( mo_num, Fock_matrix_tc_mo_tot &
                     , fock_tc_leigvec_mo, fock_tc_reigvec_mo & 
                     , n_real_tc, eigval_right_tmp )

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
!
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

  accu_d = 0.d0
  accu = 0.d0
  do i = 1, mo_num
    accu_d += overlap_fock_tc_eigvec_ao(i,i)
    do j = 1, mo_num
      accu += dabs(overlap_fock_tc_eigvec_ao(j,i) - overlap_fock_tc_eigvec_mo(j,i))
    enddo
  enddo
  print*,'Trace of the overlap_fock_tc_eigvec_ao = ',accu_d
  print*,'mo_num                                 = ',mo_num
  accu = accu / dble(mo_num**2)

  if(dabs(accu).gt.1.d-10) then
    print*,'Warning !! '
    print*,'overlap_fock_tc_eigvec_ao and overlap_fock_tc_eigvec_mo are different at '
    print*,'an average of ', accu

!     print*,'overlap_fock_tc_eigvec_ao = '
!     do i = 1, mo_num
!       write(*,'(100(F16.10,X))')overlap_fock_tc_eigvec_ao(i,:)
!     enddo
!     print*,'overlap_fock_tc_eigvec_mo = '
!     do i = 1, mo_num
!       write(*,'(100(F16.10,X))')overlap_fock_tc_eigvec_mo(i,:)
!     enddo
!     stop

  endif

END_PROVIDER

! ---

 BEGIN_PROVIDER [ double precision, bi_ortho_fock_tc_reigvec_mo, (mo_num, mo_num)]
&BEGIN_PROVIDER [ double precision, bi_ortho_fock_tc_leigvec_mo, (mo_num, mo_num)]
&BEGIN_PROVIDER [ double precision, eigval_bi_ortho_fock_tc_mo, (mo_num)]
 implicit none
 integer :: n_real_tc
 call non_hrmt_real_diag_new( mo_num, overlap_fock_tc_eigvec_ao&
! call non_hrmt_bieig_real_im( mo_num, overlap_fock_tc_eigvec_ao&
                    , bi_ortho_fock_tc_leigvec_mo, bi_ortho_fock_tc_reigvec_mo & 
                    , n_real_tc, eigval_bi_ortho_fock_tc_mo)
 do i = 1, mo_num
  print*,eigval_bi_ortho_fock_tc_mo(i)
 enddo
 integer :: i,j,k,l
 double precision, allocatable :: tmp(:,:)
 allocate(tmp(mo_num,mo_num))
 tmp = 0.d0
 print*,'bi_ortho_fock_tc_leigvec_mo'
 do i = 1, mo_num
  write(*,'(100(F16.10,X))')bi_ortho_fock_tc_leigvec_mo(:,i)
 enddo
 print*,'bi_ortho_fock_tc_reigvec_mo'
 do i = 1, mo_num
  write(*,'(100(F16.10,X))')bi_ortho_fock_tc_reigvec_mo(:,i)
 enddo
 do i = 1, mo_num
  do j = 1, mo_num
   do k = 1, mo_num
    do l = 1, mo_num
     tmp (j,i) += bi_ortho_fock_tc_leigvec_mo(l,j) * overlap_fock_tc_eigvec_ao(l,k) * bi_ortho_fock_tc_reigvec_mo(k,i)
    enddo
   enddo
  enddo
 enddo
 print*,'tmp'
do i = 1, mo_num
 write(*,'(100(F16.10,X))')tmp(:,i)
enddo
END_PROVIDER 

 BEGIN_PROVIDER [ double precision, bi_ortho_fock_tc_reigvec_ao, (mo_num, mo_num)]
&BEGIN_PROVIDER [ double precision, bi_ortho_fock_tc_leigvec_ao, (mo_num, mo_num)]
 implicit none
 double precision, allocatable :: reigvec(:,:),leigvec(:,:),overlap(:,:)
 allocate(reigvec(ao_num, mo_num), leigvec(ao_num, mo_num), overlap(mo_num, mo_num))
!  ! MO_R x R
  call dgemm( 'N', 'N', ao_num, mo_num, mo_num, 1.d0          &
            , fock_tc_reigvec_ao, size(fock_tc_reigvec_ao, 1)                   &
            , bi_ortho_fock_tc_reigvec_mo, size(bi_ortho_fock_tc_reigvec_mo, 1) &
            , 0.d0, reigvec, size(reigvec, 1) )
!  ! MO_L x L
  call dgemm( 'N', 'N', ao_num, mo_num, mo_num, 1.d0          &
            , fock_tc_leigvec_ao, size(fock_tc_leigvec_ao, 1)                   &
            , bi_ortho_fock_tc_leigvec_mo, size(bi_ortho_fock_tc_leigvec_mo, 1) &
            , 0.d0, leigvec, size(leigvec, 1) )
  integer :: i,j,p,q
 overlap = 0.d0
 do i = 1, mo_num
  do j = 1, mo_num
   do p = 1, ao_num
    do q = 1, ao_num
     overlap(j,i) += leigvec(q,j) * ao_overlap(q,p) * reigvec(p,i)
    enddo
   enddo
  enddo
 enddo
 print*,'overlap = '
 do i = 1, mo_num
  write(*,'(100(F16.10,X))')overlap(:,i)
 enddo

END_PROVIDER 

