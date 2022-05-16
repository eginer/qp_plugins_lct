 BEGIN_PROVIDER [ double precision, mo_r_coef, (ao_num, mo_num)]
 implicit none
 BEGIN_DOC
 ! Coef of the RIGHT MOs on the AOS. 
 ! 
 ! Contains normalization factors ensuring BI-ORTHONORMALITY with LEFT MOs
 END_DOC
 if(bi_ortho)then
!  mo_r_coef = Reig_tcFock_ao
  mo_r_coef = bi_ortho_fock_reigvec_ao
  integer :: i,m
  double precision :: norm
  do i = 1, mo_num
   norm = dsqrt(dabs(overlap_diag_bi_ortho_mo_coef(i)))
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
! mo_l_coef = Leig_tcFock_ao
  mo_l_coef = bi_ortho_fock_leigvec_ao
  integer :: i,m
  double precision :: norm
 do i = 1, mo_num
  norm = dsqrt(dabs(overlap_diag_bi_ortho_mo_coef(i)))
  do m = 1, ao_num
   mo_l_coef(m,i) *= 1.d0/norm
  enddo
 enddo
 else
  mo_l_coef = mo_coef
 endif
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
 print*,'****************'
 print*,'Overlap matrix betwee mo_l_coef and mo_r_coef  '
 do i = 1, mo_num
  write(*,'(100(F16.10,X))')overlap_bi_ortho(i,:)
  do k = 1, mo_num
  if(i==k)then
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

 BEGIN_PROVIDER [ double precision, mo_r_coef_un_norm, (ao_num, mo_num)]
 implicit none
! mo_r_coef_un_norm = mo_coef
! mo_r_coef_un_norm = Reig_tcFock_ao
 mo_r_coef_un_norm = fock_tc_reigvec_ao
 END_PROVIDER 


 BEGIN_PROVIDER [ double precision, mo_l_coef_un_norm, (ao_num, mo_num)]
 implicit none
! mo_l_coef_un_norm = mo_coef
! mo_l_coef_un_norm = Leig_tcFock_ao
 mo_l_coef_un_norm = fock_tc_leigvec_ao
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


 BEGIN_PROVIDER [ double precision, overlap_bi_ortho_un_norm, (mo_num, mo_num)]
&BEGIN_PROVIDER [ double precision, overlap_diag_bi_ortho_un_norm, (mo_num)]
 integer :: i,k,m,n
 overlap_bi_ortho_un_norm = 0.d0
 do i = 1, mo_num
  do k = 1, mo_num
   do m = 1, ao_num
    do n = 1, ao_num
     overlap_bi_ortho_un_norm(k,i) += ao_overlap(n,m) * mo_l_coef_un_norm(n,k) * mo_r_coef_un_norm(m,i)
    enddo
   enddo
  enddo
 enddo
 do i = 1, mo_num
  overlap_diag_bi_ortho_un_norm(i) = overlap_bi_ortho_un_norm(i,i)
 enddo
 double precision :: accu_d, accu_nd 
 accu_d = 0.d0
 accu_nd = 0.d0
 print*,'overlap_bi_ortho_un_norm'
 print*,'****************'
 print*,'Overlap matrix betwee mo_l_coef_un_norm and mo_r_coef_un_norm  '
 do i = 1, mo_num
  write(*,'(100(F16.10,X))')overlap_bi_ortho_un_norm(i,:)
 enddo
 do i = 1, mo_num
  do k = 1, mo_num
  if(i==k)then
   accu_d += dabs(overlap_bi_ortho_un_norm(k,i))
  else
   accu_nd += dabs(overlap_bi_ortho_un_norm(k,i))
  endif
  enddo 
 enddo
 print*,'accu_d  = ',accu_d
 print*,'accu_nd = ',accu_nd
 
 END_PROVIDER 

 BEGIN_PROVIDER [double precision, reigvec_overlap_diag_bi_ortho_un_norm, (mo_num, mo_num)]
&BEGIN_PROVIDER [double precision, leigvec_overlap_diag_bi_ortho_un_norm, (mo_num, mo_num)]
&BEGIN_PROVIDER [double precision, eigval_overlap_diag_bi_ortho_un_norm, (mo_num)]
 implicit none
 integer :: n_real
 call non_hrmt_real_diag(mo_num,overlap_bi_ortho_un_norm,&
         reigvec_overlap_diag_bi_ortho_un_norm,leigvec_overlap_diag_bi_ortho_un_norm,&
         n_real,eigval_overlap_diag_bi_ortho_un_norm)
 integer :: i
 do i = 1, mo_num
  print*,'eigval_overlap_diag_bi_ortho_un_norm',eigval_overlap_diag_bi_ortho_un_norm(i)
 enddo
END_PROVIDER 

 BEGIN_PROVIDER [ double precision, bi_ortho_fock_reigvec_ao, (ao_num, mo_num)]
&BEGIN_PROVIDER [ double precision, bi_ortho_fock_leigvec_ao, (ao_num, mo_num)]
 implicit none
 integer :: i,j,p,q
 bi_ortho_fock_reigvec_ao = 0.d0
 bi_ortho_fock_leigvec_ao = 0.d0
 do i = 1, mo_num
  do j = 1, mo_num
   do p = 1, ao_num
    bi_ortho_fock_reigvec_ao(p,i) += reigvec_overlap_diag_bi_ortho_un_norm(j,i) * fock_tc_reigvec_ao(p,j)
    bi_ortho_fock_leigvec_ao(p,i) += leigvec_overlap_diag_bi_ortho_un_norm(j,i) * fock_tc_leigvec_ao(p,j)
   enddo
  enddo
 enddo
  

 END_PROVIDER 

 BEGIN_PROVIDER [ double precision, overlap_bi_ortho_mo_coef, (mo_num, mo_num)]
&BEGIN_PROVIDER [ double precision, overlap_diag_bi_ortho_mo_coef, (mo_num)]
 implicit none
 integer :: i,j,p,q
 overlap_bi_ortho_mo_coef = 0.d0 
 do i = 1, mo_num
  do j = 1, mo_num
   do p = 1, ao_num
    do q = 1, ao_num
     overlap_bi_ortho_mo_coef(j,i) += ao_overlap(p,q) * bi_ortho_fock_leigvec_ao(p,j) * bi_ortho_fock_reigvec_ao(q,i)
    enddo
   enddo
  enddo
 enddo
 do i = 1, mo_num
  overlap_diag_bi_ortho_mo_coef(i) = overlap_bi_ortho_mo_coef(i,i)
 enddo
 print*,'overlap for bi_ortho_fock_reigvec_ao'
 double precision :: accu_d, accu_nd
 accu_d = 0.d0
 accu_nd = 0.d0
 do i = 1, mo_num
  write(*,'(100(F16.10,X))')overlap_bi_ortho_mo_coef(i,:)
  accu_d += dabs(overlap_bi_ortho_mo_coef(i,i))
  do j = 1, mo_num
   if(i==j)cycle
   accu_nd += dabs(overlap_bi_ortho_mo_coef(j,i))
  enddo
 enddo
 print*,'accu_d  = ',accu_d
 print*,'accu_nd = ',accu_nd
 END_PROVIDER 

 BEGIN_PROVIDER [ double precision, fock_tc_reigvec_mo, (mo_num, mo_num)]
&BEGIN_PROVIDER [ double precision, fock_tc_leigvec_mo, (mo_num, mo_num)]
 implicit none
 integer :: n_real_tc, eigval_right_tmp(mo_num)
 call non_hrmt_real_diag(mo_num,Fock_matrix_tc_mo_tot,fock_tc_reigvec_mo,fock_tc_leigvec_mo,n_real_tc,eigval_right_tmp)
 integer :: i,j,k,l
 double precision :: s_mat(mo_num, mo_num)
 s_mat = 0.d0
 do i = 1, mo_num
  do k = 1, mo_num
   do l = 1, mo_num
    s_mat(k,i) += fock_tc_reigvec_mo(l,i) * fock_tc_leigvec_mo(l,k)
   enddo
  enddo
 enddo
 double precision :: accu_d, accu_nd
 accu_d = 0.D0
 accu_nd = 0.D0
 do i = 1, mo_num
  do k = 1, mo_num
  if(i==k)then
   accu_d += s_mat(k,i)
  else
   accu_nd += s_mat(k,i)
  endif
  enddo 
 enddo
 print*,'accu_d  MO = ',accu_d
 print*,'accu_nd MO = ',accu_nd
 END_PROVIDER 

 BEGIN_PROVIDER [ double precision, fock_tc_reigvec_ao, (ao_num, mo_num)]
&BEGIN_PROVIDER [ double precision, fock_tc_leigvec_ao, (ao_num, mo_num)]
 implicit none
 integer :: i,j,q
 fock_tc_reigvec_ao = 0.d0
 fock_tc_leigvec_ao = 0.d0
 do i = 1, mo_num
  do j = 1, mo_num
   do q = 1, ao_num
    fock_tc_reigvec_ao(q,i) += mo_coef(q,j) * fock_tc_reigvec_mo(j,i)
    fock_tc_leigvec_ao(q,i) += mo_coef(q,j) * fock_tc_leigvec_mo(j,i)
   enddo
  enddo
 enddo

 END_PROVIDER