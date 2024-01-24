
 BEGIN_PROVIDER [double precision, norm_n2_jastrow_cst_mu_on_grid, (n_points_final_grid,N_states)]
&BEGIN_PROVIDER [double precision, norm_n2_jastrow_cst_mu, (N_states)]
 implicit none
 BEGIN_DOC
! norm_n2_jastrow_cst_mu_on_grid(ipoint) = \int dr2 n2(ipoint,r2)J^2(mu(ipoint)
! 
! norm_n2_jastrow_cst_mu(ipoint) = \int dr1 \int dr2 n2(ipoint,r2)J^2(mu(ipoint)
 END_DOC
 integer :: ipoint,istate,i,j,k,l
 double precision :: weight
 norm_n2_jastrow_cst_mu = 0.d0
 norm_n2_jastrow_cst_mu_on_grid = 0.d0
 do i = 1, N_det
  psi_coef(i,1) = reigvec_trans(i,1)/dsqrt(reigvec_trans_norm(1))
 enddo
 touch psi_coef

 do istate =1, N_states
  do ipoint = 1, n_points_final_grid
   weight = final_weight_at_r_vector(ipoint)
   do i = 1, mo_num ! 1
    do k = 1, mo_num ! 1 
     do j = 1, mo_num ! 2 
      do l = 1, mo_num ! 2
       norm_n2_jastrow_cst_mu_on_grid(ipoint,istate) += mos_in_r_array(i,ipoint) * mos_in_r_array(k,ipoint) * & 
                          ovlp_jastrow2_ij_cst_mu_mo(l,j,ipoint) * full_occ_2_rdm_ab_chemist_mo(l,j,k,i,istate)
      enddo
     enddo
    enddo
   enddo
   norm_n2_jastrow_cst_mu(istate) += norm_n2_jastrow_cst_mu_on_grid(ipoint,istate) * weight
  enddo
 enddo
END_PROVIDER 


BEGIN_PROVIDER [double precision, ovlp_jastrow2_ij_cst_mu_mo, (mo_num,mo_num,n_points_final_grid)]
 implicit none
 BEGIN_DOC
! \int dr2 phi_i(r2) \phi_j(r2) e^{J(r1,r2,\mu(r1))}
 END_DOC
 integer :: i,j,ipoint,n_taylor
 double precision :: r1(3),mu,mo_ints(mo_num,mo_num)
 n_taylor = 4
 double precision :: exponent_exp
 exponent_exp = 1.d0
 do ipoint = 1, n_points_final_grid
  r1(1) = final_grid_points(1,ipoint)
  r1(2) = final_grid_points(2,ipoint)
  r1(3) = final_grid_points(3,ipoint)
  mu = mu_erf 
  call give_jastrow2_ovlp_ints_mo(mu,r1,n_taylor,mo_ints,exponent_exp)
  do i = 1, mo_num
   do j = 1, mo_num
    ovlp_jastrow2_ij_cst_mu_mo(j,i,ipoint) = mo_ints(j,i)
   enddo
  enddo
 enddo

END_PROVIDER 


 BEGIN_PROVIDER [double precision, norm_n2_inv_jastrow_cst_mu_on_grid, (n_points_final_grid,N_states)]
&BEGIN_PROVIDER [double precision, norm_n2_inv_jastrow_cst_mu, (N_states)]
 implicit none
 BEGIN_DOC
! norm_n2_inv_jastrow_cst_mu_on_grid(ipoint) = \int dr2 n2(ipoint,r2)J^2(mu(ipoint)
! 
! norm_n2_inv_jastrow_cst_mu(ipoint) = \int dr1 \int dr2 n2(ipoint,r2)J^2(mu(ipoint)
 END_DOC
 integer :: ipoint,istate,i,j,k,l
 double precision :: weight
 norm_n2_inv_jastrow_cst_mu = 0.d0
 norm_n2_inv_jastrow_cst_mu_on_grid = 0.d0
 do i = 1, N_det
  psi_coef(i,1) = leigvec_trans(i,1)/dsqrt(leigvec_trans_norm(1))
 enddo
 touch psi_coef

 do istate =1, N_states
  do ipoint = 1, n_points_final_grid
   weight = final_weight_at_r_vector(ipoint)
   do i = 1, mo_num ! 1
    do k = 1, mo_num ! 1 
     do j = 1, mo_num ! 2 
      do l = 1, mo_num ! 2
       norm_n2_inv_jastrow_cst_mu_on_grid(ipoint,istate) += mos_in_r_array(i,ipoint) * mos_in_r_array(k,ipoint) * & 
                          ovlp_inv_jastrow2_ij_cst_mu_mo(l,j,ipoint) * full_occ_2_rdm_ab_chemist_mo(l,j,k,i,istate)
      enddo
     enddo
    enddo
   enddo
   norm_n2_inv_jastrow_cst_mu(istate) += norm_n2_inv_jastrow_cst_mu_on_grid(ipoint,istate) * weight
  enddo
 enddo
END_PROVIDER 


BEGIN_PROVIDER [double precision, ovlp_inv_jastrow2_ij_cst_mu_mo, (mo_num,mo_num,n_points_final_grid)]
 implicit none
 BEGIN_DOC
! \int dr2 phi_i(r2) \phi_j(r2) e^{J(r1,r2,\mu(r1))}
 END_DOC
 integer :: i,j,ipoint,n_taylor
 double precision :: r1(3),mu,mo_ints(mo_num,mo_num)
 n_taylor = 4
 double precision :: exponent_exp
 exponent_exp = -1.d0
 do ipoint = 1, n_points_final_grid
  r1(1) = final_grid_points(1,ipoint)
  r1(2) = final_grid_points(2,ipoint)
  r1(3) = final_grid_points(3,ipoint)
  mu = mu_erf 
  call give_jastrow2_ovlp_ints_mo(mu,r1,n_taylor,mo_ints,exponent_exp)
  do i = 1, mo_num
   do j = 1, mo_num
    ovlp_inv_jastrow2_ij_cst_mu_mo(j,i,ipoint) = mo_ints(j,i)
   enddo
  enddo
 enddo

END_PROVIDER 

