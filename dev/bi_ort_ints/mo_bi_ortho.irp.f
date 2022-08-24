
 BEGIN_PROVIDER [ double precision, overlap_bi_ortho, (mo_num, mo_num)]
&BEGIN_PROVIDER [ double precision, overlap_diag_bi_ortho, (mo_num)]

  BEGIN_DOC
  ! Overlap matrix between the RIGHT and LEFT MOs. Should be the identity matrix 
  END_DOC

  implicit none
  integer                       :: i, k, m, n
  double precision              :: accu_d, accu_nd 
  double precision, allocatable :: tmp(:,:)

  overlap_bi_ortho = 0.d0
  do i = 1, mo_num
    do k = 1, mo_num
      do m = 1, ao_num
        do n = 1, ao_num
          overlap_bi_ortho(k,i) += ao_overlap(n,m) * mo_l_coef(n,k) * mo_r_coef(m,i)
        enddo
      enddo
    enddo
  enddo

!  allocate( tmp(mo_num,ao_num) )
!
!  ! tmp <-- L.T x S_ao
!  call dgemm( "T", "N", mo_num, ao_num, ao_num, 1.d0                         & 
!            , mo_l_coef, size(mo_l_coef, 1), ao_overlap, size(ao_overlap, 1) &
!            , 0.d0, tmp, size(tmp, 1) )
!
!  ! S <-- tmp x R
!  call dgemm( "N", "N", mo_num, mo_num, ao_num, 1.d0           & 
!            , tmp, size(tmp, 1), mo_r_coef, size(mo_r_coef, 1) &
!            , 0.d0, overlap_bi_ortho, size(overlap_bi_ortho, 1) )
!
!  deallocate( tmp )

  do i = 1, mo_num
    overlap_diag_bi_ortho(i) = overlap_bi_ortho(i,i)
  enddo

  accu_d  = 0.d0
  accu_nd = 0.d0
  do i = 1, mo_num
    do k = 1, mo_num
      if(i==k) then
        accu_d += dabs(overlap_bi_ortho(k,i))
      else
        accu_nd += dabs(overlap_bi_ortho(k,i))
      endif
    enddo 
  enddo
  accu_d = accu_d/dble(mo_num)
  accu_nd = accu_nd/dble(mo_num**2-mo_num)
  if(dabs(accu_d-1.d0).gt.1.d-10.or.dabs(accu_nd).gt.1.d-10)then
    print*,'Warning !!!'
    print*,'Average trace of overlap_bi_ortho is different from 1 by ', accu_d
    print*,'And bi orthogonality is off by an average of ',accu_nd
    print*,'****************'
    print*,'Overlap matrix betwee mo_l_coef and mo_r_coef  '
    do i = 1, mo_num
      write(*,'(100(F16.10,X))')overlap_bi_ortho(i,:)
    enddo
  endif
  print*,'Average trace of overlap_bi_ortho (should be 1.)'
  print*,'accu_d  = ',accu_d
  print*,'Sum of off diagonal terms of overlap_bi_ortho (should be zero)'
  print*,'accu_nd = ',accu_nd
  print*,'****************'
 
END_PROVIDER 

! ---

