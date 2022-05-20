
! ---

BEGIN_PROVIDER [ double precision, mo_r_coef, (ao_num, mo_num)]

  BEGIN_DOC
  ! Coef of the RIGHT MOs on the AOS. 
  ! 
  ! Contains normalization factors ensuring BI-ORTHONORMALITY with LEFT MOs
  END_DOC

  implicit none
  integer          :: i, m
  double precision :: norm

  if(bi_ortho)then
    !mo_r_coef = Reig_tcFock_ao
    mo_r_coef = fock_tc_reigvec_ao
    do i = 1, mo_num
      norm = dsqrt(dabs(overlap_fock_tc_eigvec_ao(i,i)))
      do m = 1, ao_num
        mo_r_coef(m,i) *= 1.d0/norm
      enddo
    enddo
  else
    mo_r_coef = mo_coef
  endif

END_PROVIDER 

! ---

BEGIN_PROVIDER [ double precision, mo_l_coef, (ao_num, mo_num)]

  BEGIN_DOC
  ! Coef of the LEFT  MOs on the AOS. 
  ! 
  ! Contains normalization factors ensuring BI-ORTHONORMALITY with RIGHT MOs
  END_DOC

  implicit none
  integer          :: i, m
  double precision :: norm

  if(bi_ortho) then
    ! mo_l_coef = Leig_tcFock_ao
    mo_l_coef = fock_tc_leigvec_ao
    do i = 1, mo_num
      norm = dsqrt(dabs(overlap_fock_tc_eigvec_ao(i,i)))
      do m = 1, ao_num
        mo_l_coef(m,i) *= 1.d0/norm
      enddo
    enddo
  else
    mo_l_coef = mo_coef
  endif

END_PROVIDER 

! ---

 BEGIN_PROVIDER [ double precision, mo_r_coef_transp, (mo_num, ao_num)]
 implicit none
 integer :: j,m
 do j = 1, mo_num
  do m = 1, ao_num
   mo_r_coef_transp(j,m)  = mo_r_coef(m,j)
  enddo
 enddo
 END_PROVIDER 


 BEGIN_PROVIDER [ double precision, mo_l_coef_transp, (mo_num, ao_num)]
 implicit none
 integer :: j,m
 do j = 1, mo_num
  do m = 1, ao_num
   mo_l_coef_transp(j,m)  = mo_l_coef(m,j)
  enddo
 enddo
 END_PROVIDER 

! ---

 BEGIN_PROVIDER [ double precision, overlap_bi_ortho, (mo_num, mo_num)]
&BEGIN_PROVIDER [ double precision, overlap_diag_bi_ortho, (mo_num)]

  BEGIN_DOC
  ! Overlap matrix between the RIGHT and LEFT MOs. Should be the identity matrix 
  END_DOC

  implicit none
  integer                       :: i, k, m, n
  double precision              :: accu_d, accu_nd 
  double precision, allocatable :: tmp(:,:)

  !overlap_bi_ortho = 0.d0
  !do i = 1, mo_num
  !  do k = 1, mo_num
  !    do m = 1, ao_num
  !      do n = 1, ao_num
  !        overlap_bi_ortho(k,i) += ao_overlap(n,m) * mo_l_coef(n,k) * mo_r_coef(m,i)
  !      enddo
  !    enddo
  !  enddo
  !enddo

  allocate( tmp(mo_num,ao_num) )

  ! tmp <-- L.T x S_ao
  call dgemm( "T", "N", mo_num, ao_num, ao_num, 1.d0                         & 
            , mo_l_coef, size(mo_l_coef, 1), ao_overlap, size(ao_overlap, 1) &
            , 0.d0, tmp, size(tmp, 1) )

  ! S <-- tmp x R
  call dgemm( "N", "N", mo_num, mo_num, ao_num, 1.d0           & 
            , tmp, size(tmp, 1), mo_r_coef, size(mo_r_coef, 1) &
            , 0.d0, overlap_bi_ortho, size(overlap_bi_ortho, 1) )

  deallocate( tmp )

  do i = 1, mo_num
    overlap_diag_bi_ortho(i) = overlap_bi_ortho(i,i)
  enddo

  accu_d  = 0.d0
  accu_nd = 0.d0
  print*,'****************'
  print*,'Overlap matrix betwee mo_l_coef and mo_r_coef  '
  do i = 1, mo_num
    write(*,'(100(F16.10,X))')overlap_bi_ortho(i,:)
    do k = 1, mo_num
      if(i==k) then
        accu_d += dabs(overlap_bi_ortho(k,i))
      else
        accu_nd += dabs(overlap_bi_ortho(k,i))
      endif
    enddo  
  enddo
  print*,'accu_d  = ',accu_d
  print*,'accu_nd = ',accu_nd
  print*,'****************'
 
END_PROVIDER 

! ---

