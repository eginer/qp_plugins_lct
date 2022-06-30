
 BEGIN_PROVIDER [ double precision, inv_overlap_mo_r_mo, (mo_num, mo_num)]
&BEGIN_PROVIDER [ double precision, inv_overlap_mo_l_mo, (mo_num, mo_num)]
 implicit none
 integer :: i
 call get_inverse(overlap_mo_r_mo,mo_num,mo_num,inv_overlap_mo_r_mo,mo_num)
 call get_inverse(overlap_mo_l_mo,mo_num,mo_num,inv_overlap_mo_l_mo,mo_num)
END_PROVIDER 

BEGIN_PROVIDER [ double precision, mo_overlap_from_mo_r, (mo_num,mo_num)]
 implicit none
 integer :: i,j,k,l
 mo_overlap_from_mo_r = 0.d0
 do i = 1, mo_num
  do j = 1, mo_num
   do k = 1, mo_num
    do l = 1, mo_num
     mo_overlap_from_mo_r(j,i) += inv_overlap_mo_r_mo(k,i) * inv_overlap_mo_r_mo(l,j) * overlap_mo_r(l,k)
    enddo
   enddo
  enddo
 enddo
END_PROVIDER 

BEGIN_PROVIDER [ double precision, mo_overlap_from_mo_r_l, (mo_num,mo_num)]
 implicit none
 integer :: i,j,k,l
 mo_overlap_from_mo_r_l = 0.d0
 do i = 1, mo_num
  do j = 1, mo_num
   do k = 1, mo_num
    do l = 1, mo_num
     mo_overlap_from_mo_r_l(j,i) += inv_overlap_mo_l_mo(k,i) * inv_overlap_mo_r_mo(l,j) * overlap_bi_ortho(k,l)
    enddo
   enddo
  enddo
 enddo
END_PROVIDER 

 BEGIN_PROVIDER [ double precision, F_mat_tc_mo_tot_ortho_basis, (mo_num, mo_num)]
 implicit none
 integer :: i,j,k,l
!! F_MO,MO = P_MO,mo_l F_mo_l,mo_r P_mo_r,MO
!! |mo_i> = \sum_k,l <mo_r_k|mo_i> | mo_r_l> <mo_r_k|mo_r_l>
!! <mo_l_k|F|mo_i> = \sum_i,j <mo_r_j|mo_i> | mo_r_j> <mo_r_i|mo_r_j>
 allocate(tmp(mo_num, mo_num))
 double precision, allocatable :: tmp(:,:)
! call dgemm( 'N', 'N', mo_num, mo_num, mo_num, 1.d0          &
!           ,overlap_mo_l_mo, size(overlap_mo_l_mo, 1)                   &
!           , Fock_matrix_tc_mo_tot, size(Fock_matrix_tc_mo_tot, 1) &
!           , 0.d0, tmp, size(tmp, 1) )
! call dgemm( 'N', 'N', mo_num, mo_num, mo_num, 1.d0          &
!           , tmp, size(tmp, 1)                   &
!           , inv_overlap_mo_r_mo, size(inv_overlap_mo_r_mo, 1) &
!           , 0.d0, F_mat_tc_mo_tot_ortho_basis, size(F_mat_tc_mo_tot_ortho_basis, 1) )
! tmp = 0.d0
! do i = 1, mo_num
!  do j = 1, mo_num
!   do k = 1, mo_num
!    tmp(j,i) += overlap_mo_l_mo(j,k) * Fock_matrix_tc_mo_tot(k,i)
!   enddo
!  enddo
! enddo
 F_mat_tc_mo_tot_ortho_basis = 0.d0
! do i = 1, mo_num
!  do j = 1, mo_num
!   do k = 1, mo_num
!    F_mat_tc_mo_tot_ortho_basis(j,i) += tmp(j,k) * inv_overlap_mo_r_mo(k,i)
!   enddo
!  enddo
! enddo
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
 implicit none
 integer :: n_real_tc
  call non_hrmt_bieig_real_im( mo_num, F_mat_tc_mo_tot_ortho_basis&
                     , leigvec_f_mat_ortho_basis, reigvec_f_mat_ortho_basis &
                     , n_real_tc, eigval_f_mat_ortho_basis)
 print*,'n_real_tc= ',n_real_tc
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
