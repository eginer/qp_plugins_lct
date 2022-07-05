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

 BEGIN_PROVIDER [ double precision, inv_overlap_mo_r_mo, (mo_num, mo_num)]
&BEGIN_PROVIDER [ double precision, inv_overlap_mo_l_mo, (mo_num, mo_num)]
 implicit none
 integer :: i
 call get_inverse(overlap_mo_r_mo,mo_num,mo_num,inv_overlap_mo_r_mo,mo_num)
 call get_inverse(overlap_mo_l_mo,mo_num,mo_num,inv_overlap_mo_l_mo,mo_num)
END_PROVIDER 


 BEGIN_PROVIDER [ double precision, F_mat_tc_mo_tot_ortho_basis, (mo_num, mo_num)]
 implicit none
 integer :: i,j,k,l
!! F_MO,MO = P_MO,mo_l F_mo_l,mo_r P_mo_r,MO
!! |mo_i> = \sum_k,l <mo_r_k|mo_i> | mo_r_l> <mo_r_k|mo_r_l>
!! <mo_l_k|F|mo_i> = \sum_i,j <mo_r_j|mo_i> | mo_r_j> <mo_r_i|mo_r_j>
 allocate(tmp(mo_num, mo_num))
 double precision, allocatable :: tmp(:,:)
 F_mat_tc_mo_tot_ortho_basis = 0.d0
 do i = 1, mo_num
  do j = 1, mo_num
   do k = 1, mo_num
    do l = 1, mo_num
     F_mat_tc_mo_tot_ortho_basis(j,i) += inv_overlap_mo_l_mo(k,i) * inv_overlap_mo_r_mo(l,j) * Fock_matrix_tc_mo_tot(k,l)
    enddo
   enddo
  enddo
 enddo

 END_PROVIDER 
 
 BEGIN_PROVIDER [ double precision, reigvec_f_mat_ortho_basis, (mo_num, mo_num) ]
&BEGIN_PROVIDER [ double precision, leigvec_f_mat_ortho_basis, (mo_num, mo_num) ]
&BEGIN_PROVIDER [ double precision, eigval_f_mat_ortho_basis, (mo_num) ]
&BEGIN_PROVIDER [ double precision, overlap_rleigv_f_mat_ortho_basis, (mo_num, mo_num)]
 implicit none
 integer :: n_real_tc
  call non_hrmt_real_im( mo_num, F_mat_tc_mo_tot_ortho_basis&
                     , leigvec_f_mat_ortho_basis, reigvec_f_mat_ortho_basis &
                     , n_real_tc, eigval_f_mat_ortho_basis)
 integer :: i,j,k,l
 double precision :: accu
 ! Normalizing vectors 
 do i = 1, mo_num
  accu = 0.d0
  do k = 1, mo_num
   accu += leigvec_f_mat_ortho_basis(k,i) * reigvec_f_mat_ortho_basis(k,i)
  enddo
  accu = 1.d0/dsqrt(dabs(accu))
  do k = 1, mo_num
   leigvec_f_mat_ortho_basis(k,i) *= accu
   reigvec_f_mat_ortho_basis(k,i) *= accu
  enddo
 enddo
 overlap_rleigv_f_mat_ortho_basis = 0.d0
 do i = 1, mo_num
  do j = 1, mo_num
   do k = 1, mo_num
    overlap_rleigv_f_mat_ortho_basis(j,i) += leigvec_f_mat_ortho_basis(k,j) * reigvec_f_mat_ortho_basis(k,i)
   enddo
  enddo
 enddo
 print*,'overlap_rleigv_f_mat_ortho_basis'
 accu = 0.d0
 do i = 1,mo_num
  do j = 1, mo_num
   if(i==j)cycle
   accu += dabs(overlap_rleigv_f_mat_ortho_basis(j,i))
  enddo
  write(*,'(100(F16.10,X))')overlap_rleigv_f_mat_ortho_basis(:,i)
 enddo
 print*,'sum non diagonal elements ',accu
 double precision :: mat_tmp(mo_num, mo_num)
 mat_tmp = 0.d0
 do i = 1, mo_num
  do j = 1, mo_num
   do k = 1, mo_num
    do l = 1, mo_num
     mat_tmp(j,i) += leigvec_f_mat_ortho_basis(k,j) * F_mat_tc_mo_tot_ortho_basis(k,l) * reigvec_f_mat_ortho_basis(l,i) 
    enddo
   enddo
  enddo
 enddo
 print*,'mat_tmp'
 accu = 0.d0
 do i = 1, mo_num
  do j = 1, mo_num
   if(i==j)cycle
   accu += dabs(mat_tmp(j,i))
  enddo
  write(*,'(100(F16.10,X))')mat_tmp(:,i)
 enddo
 print*,'sum non diagonal elements ',accu
 END_PROVIDER 

 BEGIN_PROVIDER [ double precision, overlap_mo_r_mo, (mo_num, mo_num)]
&BEGIN_PROVIDER [ double precision, overlap_mo_l_mo, (mo_num, mo_num)]
 implicit none
 BEGIN_DOC
! overlap_mo_r_mo(j,i) = <MO_j|MO_R_i>
 END_DOC
 integer :: i,j,p,q
 overlap_mo_r_mo = 0.d0
 overlap_mo_l_mo = 0.d0
 do i = 1, mo_num
   do j = 1, mo_num
    do p = 1, ao_num
     do q = 1, ao_num
      overlap_mo_r_mo(j,i) += mo_coef(p,j) * mo_r_coef(q,i) * ao_overlap(q,p)
      overlap_mo_l_mo(j,i) += mo_coef(p,j) * mo_l_coef(q,i) * ao_overlap(q,p)
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
! fock_tc_reigvec_mo_ortho(j,i) = <MO_i|F_R_j>
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
