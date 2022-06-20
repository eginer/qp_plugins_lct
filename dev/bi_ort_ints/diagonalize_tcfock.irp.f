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

  !print*,'********'
  !print*,'Diagonal values of Fock_matrix_tc_mo_tot'
  !do i = 1, mo_num
  ! write(*,'(100(F15.10,X))') Fock_matrix_tc_mo_tot(i,i)
  !enddo
  !print*,'********'

!   call non_hrmt_real_diag_new( mo_num, Fock_matrix_tc_mo_tot &
  call non_hrmt_bieig( mo_num, Fock_matrix_tc_mo_tot          &
                     , fock_tc_leigvec_mo, fock_tc_reigvec_mo & 
                     , n_real_tc, eigval_right_tmp )

  eigval_fock_tc_mo = eigval_right_tmp
  !print*,'Eigenvalues of Fock_matrix_tc_mo_tot'
  !do i = 1, mo_num
  !  print*, i, eigval_fock_tc_mo(i)
  !enddo
  deallocate( eigval_right_tmp )

  !overlap_fock_tc_eigvec_mo = 0.d0
  !do i = 1, mo_num
  !  do k = 1, mo_num
  !    do l = 1, mo_num
  !      overlap_fock_tc_eigvec_mo(k,i) += fock_tc_leigvec_mo(l,k) * fock_tc_reigvec_mo(l,i) 
  !    enddo
  !  enddo
  !enddo

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

  if( accu_nd .gt. 1d-7 ) then
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


 BEGIN_PROVIDER [ double precision, overlap_mo_r_mo, (mo_num, mo_num)]
&BEGIN_PROVIDER [ double precision, overlap_mo_l_mo, (mo_num, mo_num)]
 implicit none
 BEGIN_DOC
! overlap_mo_r_mo(j,i) = <MO_i|MO_R_j>
 END_DOC
 integer :: i,j,p,q
 overlap_mo_r_mo = 0.d0
 overlap_mo_l_mo = 0.d0
 do i = 1, mo_num
   do j = 1, mo_num
    do p = 1, ao_num
     do q = 1, ao_num
      overlap_mo_r_mo(j,i) += mo_r_coef(q,i) * mo_coef(p,j) * ao_overlap(q,p) 
      overlap_mo_l_mo(j,i) += mo_l_coef(q,i) * mo_coef(p,j) * ao_overlap(q,p)
     enddo
    enddo
   enddo
 enddo
END_PROVIDER 

 BEGIN_PROVIDER [ double precision, overlap_mo_r, (mo_num, mo_num)]
&BEGIN_PROVIDER [ double precision, overlap_mo_l, (mo_num, mo_num)]
 implicit none
 BEGIN_DOC
! overlap_mo_r_mo(j,i) = <MO_i|MO_R_j>
 END_DOC
 integer :: i,j,p,q
 overlap_mo_r= 0.d0
 overlap_mo_l= 0.d0
 do i = 1, mo_num
   do j = 1, mo_num
    do p = 1, ao_num
     do q = 1, ao_num
      overlap_mo_r(j,i) += mo_r_coef(q,i) * mo_r_coef(p,j) * ao_overlap(q,p) 
      overlap_mo_l(j,i) += mo_l_coef(q,i) * mo_l_coef(p,j) * ao_overlap(q,p)
     enddo
    enddo
   enddo
 enddo
END_PROVIDER 

 BEGIN_PROVIDER [ double precision, fock_tc_reigvec_mo_ortho, (mo_num, mo_num)]
&BEGIN_PROVIDER [ double precision, fock_tc_leigvec_mo_ortho, (mo_num, mo_num)]
 implicit none
 BEGIN_DOC
! Expansion of the right- and left-eigenvectors of the fock matrix on the usual orthonormal MO basis
!
! fock_tc_reigvec_mo_ortho(j,i) = <MO_i|MO_R_j>
 END_DOC
 integer :: i,j,k
 fock_tc_reigvec_mo_ortho = 0.d0
 fock_tc_leigvec_mo_ortho = 0.d0
 do i = 1, mo_num
  do j = 1, mo_num 
   do k = 1, mo_num
    fock_tc_reigvec_mo_ortho(j,i) += fock_tc_reigvec_mo(k,i) * overlap_mo_r_mo(j,k)
    fock_tc_leigvec_mo_ortho(j,i) += fock_tc_leigvec_mo(k,i) * overlap_mo_l_mo(j,k)
   enddo
  enddo
 enddo
END_PROVIDER 

BEGIN_PROVIDER [ double precision, overlap_fock_tc_lreigvec_mo_ortho, (mo_num, mo_num)]
 implicit none
 integer :: i,j,k
 overlap_fock_tc_lreigvec_mo_ortho = 0.D0
 do i = 1, mo_num
  do j = 1, mo_num
   do k = 1, mo_num
    overlap_fock_tc_lreigvec_mo_ortho(j,i) += fock_tc_leigvec_mo_ortho(k,j) * fock_tc_reigvec_mo_ortho(k,i)
   enddo
  enddo
 enddo
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
  double precision              :: accu 
  double precision, allocatable :: tmp(:,:)

!  call rotate_mo_coef(mo_r_coef, fock_tc_reigvec_mo, overlap_mo_r, fock_tc_reigvec_ao)
!  call rotate_mo_coef(mo_l_coef, fock_tc_leigvec_mo, overlap_mo_l, fock_tc_leigvec_ao)


 fock_tc_reigvec_ao = 0.d0
 fock_tc_leigvec_ao = 0.d0
 do i = 1, mo_num
  do p = 1, ao_num
   do k = 1, mo_num
    fock_tc_reigvec_ao(p,i) += mo_coef(p,k) * fock_tc_reigvec_mo_ortho(k,i)
    fock_tc_leigvec_ao(p,i) += mo_coef(p,k) * fock_tc_leigvec_mo_ortho(k,i)
   enddo
  enddo
 enddo
  ! MO_R x R
!  call dgemm( 'N', 'N', ao_num, mo_num, mo_num, 1.d0          &
!            , mo_r_coef, size(mo_r_coef, 1)                   &
!            , fock_tc_reigvec_mo, size(fock_tc_reigvec_mo, 1) &
!            , 0.d0, fock_tc_reigvec_ao, size(fock_tc_reigvec_ao, 1) )

!  ! MO_L x L
!  call dgemm( 'N', 'N', ao_num, mo_num, mo_num, 1.d0          &
!            , mo_l_coef, size(mo_l_coef, 1)                   &
!            , fock_tc_leigvec_mo, size(fock_tc_leigvec_mo, 1) &
!            , 0.d0, fock_tc_leigvec_ao, size(fock_tc_leigvec_ao, 1) )


   overlap_fock_tc_eigvec_ao = 0.d0  
   do i = 1, mo_num 
     do k = 1, mo_num
       do p = 1, ao_num
         do q = 1, ao_num
           overlap_fock_tc_eigvec_ao(k,i) += fock_tc_leigvec_ao(q,k) * ao_overlap(q,p) * fock_tc_reigvec_ao(p,i) 
         enddo
       enddo
     enddo
   enddo

!  allocate( tmp(mo_num,ao_num) )
!
!  ! tmp <-- L.T x S_ao
!  call dgemm( "T", "N", mo_num, ao_num, ao_num, 1.d0                                           &
!            , fock_tc_leigvec_ao, size(fock_tc_leigvec_ao, 1), ao_overlap, size(ao_overlap, 1) &
!            , 0.d0, tmp, size(tmp, 1) )
!
!  ! S <-- tmp x R
!  call dgemm( "N", "N", mo_num, mo_num, ao_num, 1.d0                             &
!            , tmp, size(tmp, 1), fock_tc_reigvec_ao, size(fock_tc_reigvec_ao, 1) &
!            , 0.d0, overlap_fock_tc_eigvec_ao, size(overlap_fock_tc_eigvec_ao, 1) )
!
!  deallocate( tmp )

  ! ---

  accu = 0.d0
  do i = 1, mo_num
    do j = 1, mo_num
      accu += dabs(overlap_fock_tc_eigvec_ao(j,i) - overlap_fock_tc_eigvec_mo(j,i))
    enddo
  enddo
  accu = accu / dble(mo_num**2)

  if(dabs(accu).gt.1.d-10) then
    print*,'Warning !! '
    print*,'overlap_fock_tc_eigvec_ao and overlap_fock_tc_eigvec_mo are different at '
    print*,'an average of ', accu

     print*,'overlap_fock_tc_eigvec_ao = '
     do i = 1, mo_num
       write(*,'(100(F16.10,X))')overlap_fock_tc_eigvec_ao(i,:)
     enddo
     print*,'overlap_fock_tc_eigvec_mo = '
     do i = 1, mo_num
       write(*,'(100(F16.10,X))')overlap_fock_tc_eigvec_mo(i,:)
     enddo
     stop

  endif

END_PROVIDER

! ---

