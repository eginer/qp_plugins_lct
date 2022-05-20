 BEGIN_PROVIDER [ double precision, mo_r_coef, (ao_num, mo_num)]
 implicit none
 BEGIN_DOC
 ! Coef of the RIGHT MOs on the AOS. 
 ! 
 ! Contains normalization factors ensuring BI-ORTHONORMALITY with LEFT MOs
 END_DOC
 if(bi_ortho)then
  mo_r_coef = fock_tc_reigvec_ao
  integer :: i,m
  double precision :: norm
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


 BEGIN_PROVIDER [ double precision, mo_l_coef, (ao_num, mo_num)]
 implicit none
 BEGIN_DOC
 ! Coef of the LEFT  MOs on the AOS. 
 ! 
 ! Contains normalization factors ensuring BI-ORTHONORMALITY with RIGHT MOs
 END_DOC
 if(bi_ortho)then
  mo_l_coef = fock_tc_leigvec_ao
  integer :: i,m
  double precision :: norm
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

 BEGIN_PROVIDER [ double precision, overlap_bi_ortho, (mo_num, mo_num)]
&BEGIN_PROVIDER [ double precision, overlap_diag_bi_ortho, (mo_num)]
 integer :: i,k,m,n
 BEGIN_DOC
! Overlap matrix between the RIGHT and LEFT MOs. Should be the identity matrix 
 END_DOC
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
 do i = 1, mo_num
  overlap_diag_bi_ortho(i) = overlap_bi_ortho(i,i)
 enddo
 double precision :: accu_d, accu_nd 
 accu_d = 0.d0
 accu_nd = 0.d0
 do i = 1, mo_num
  do k = 1, mo_num
   if(i==k)then
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

