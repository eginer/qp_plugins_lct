
BEGIN_PROVIDER [double precision, ovlp_jastrow2_ij_mo, (mo_num,mo_num,n_points_final_grid)]
 implicit none
 BEGIN_DOC
! \int dr2 phi_i(r2) \phi_j(r2) e^{J(r1,r2,\mu(r1))}
 END_DOC
 integer :: i,j,ipoint,n_taylor
 double precision :: r1(3),mu,mo_ints(mo_num,mo_num)
 n_taylor = 4
 do ipoint = 1, n_points_final_grid
  r1(1) = final_grid_points(1,ipoint)
  r1(2) = final_grid_points(2,ipoint)
  r1(3) = final_grid_points(3,ipoint)
  mu = mu_of_r_prov(ipoint,1)
  call give_jastrow2_ovlp_ints_mo(mu,r1,n_taylor,mo_ints)
  do i = 1, mo_num
   do j = 1, mo_num
    ovlp_jastrow2_ij_mo(j,i,ipoint) = mo_ints(j,i)
   enddo
  enddo
 enddo

END_PROVIDER 

BEGIN_PROVIDER [double precision, coulomb_jastrow2_ij_mo, (mo_num,mo_num,n_points_final_grid)]
 implicit none
 BEGIN_DOC
! \int dr2 phi_i(r2) \phi_j(r2) e^{J(r1,r2,\mu(r1))} 1/r12
 END_DOC
 integer :: i,j,ipoint,n_taylor
 double precision :: r1(3),muj,muc,mo_ints(mo_num,mo_num)
 n_taylor = 4
 muc = 1500.d0
 do ipoint = 1, n_points_final_grid
  r1(1) = final_grid_points(1,ipoint)
  r1(2) = final_grid_points(2,ipoint)
  r1(3) = final_grid_points(3,ipoint)
  muj = mu_of_r_prov(ipoint,1)
  call give_jastrow2_erf_ints_mo(muj,muc,r1,n_taylor,mo_ints)
  do i = 1, mo_num
   do j = 1, mo_num
    coulomb_jastrow2_ij_mo(j,i,ipoint) = mo_ints(j,i)
   enddo
  enddo
 enddo

END_PROVIDER 

BEGIN_PROVIDER [double precision, erf_mu_of_r_ij_mo, (mo_num,mo_num,n_points_final_grid)]
 implicit none
 BEGIN_DOC
! \int dr2 phi_i(r2) \phi_j(r2)  erf(\mu(r1) r12) / r12
 END_DOC
 integer :: i,j,ipoint,n_taylor
 double precision :: r1(3),muj,muc,mo_ints(mo_num,mo_num)
 n_taylor = 4
 muj = 15000.d0
 do ipoint = 1, n_points_final_grid
  r1(1) = final_grid_points(1,ipoint)
  r1(2) = final_grid_points(2,ipoint)
  r1(3) = final_grid_points(3,ipoint)
  muc = mu_of_r_prov(ipoint,1)
  call give_jastrow2_erf_ints_mo(muj,muc,r1,n_taylor,mo_ints)
  do i = 1, mo_num
   do j = 1, mo_num
    erf_mu_of_r_ij_mo(j,i,ipoint) = mo_ints(j,i)
   enddo
  enddo
 enddo

END_PROVIDER 

BEGIN_PROVIDER [double precision, erf_mu_of_r_jastrow2_ij_mo, (mo_num,mo_num,n_points_final_grid)]
 implicit none
 BEGIN_DOC
! \int dr2 phi_i(r2) \phi_j(r2) J(r1,r2,\mu(r1)) erf(\mu(r1) r12) / r12
 END_DOC
 integer :: i,j,ipoint,n_taylor
 double precision :: r1(3),muj,muc,mo_ints(mo_num,mo_num)
 n_taylor = 4
 do ipoint = 1, n_points_final_grid
  r1(1) = final_grid_points(1,ipoint)
  r1(2) = final_grid_points(2,ipoint)
  r1(3) = final_grid_points(3,ipoint)
  muc = mu_of_r_prov(ipoint,1)
  muj = muc
  call give_jastrow2_erf_ints_mo(muj,muc,r1,n_taylor,mo_ints)
  do i = 1, mo_num
   do j = 1, mo_num
    erf_mu_of_r_jastrow2_ij_mo(j,i,ipoint) = mo_ints(j,i)
   enddo
  enddo
 enddo

END_PROVIDER 

 BEGIN_PROVIDER [double precision, norm_n2_jastrow_on_grid, (n_points_final_grid,N_states)]
&BEGIN_PROVIDER [double precision, norm_n2_jastrow, (N_states)]
 implicit none
 BEGIN_DOC
! norm_n2_jastrow_on_grid(ipoint) = \int dr2 n2(ipoint,r2)J^2(mu(ipoint)
! 
! norm_n2_jastrow(ipoint) = \int dr1 \int dr2 n2(ipoint,r2)J^2(mu(ipoint)
 END_DOC
 integer :: ipoint,istate,i,j,k,l
 double precision :: weight
 norm_n2_jastrow = 0.d0
 norm_n2_jastrow_on_grid = 0.d0
 do istate =1, N_states
  do ipoint = 1, n_points_final_grid
   weight = final_weight_at_r_vector(ipoint)
   do i = 1, mo_num ! 1
    do k = 1, mo_num ! 1 
     do j = 1, mo_num ! 2 
      do l = 1, mo_num ! 2
       norm_n2_jastrow_on_grid(ipoint,istate) += mos_in_r_array(i,ipoint) * mos_in_r_array(k,ipoint) * & 
                          ovlp_jastrow2_ij_mo(l,j,ipoint) * full_occ_2_rdm_ab_chemist_mo(l,j,k,i,istate)
      enddo
     enddo
    enddo
   enddo
   norm_n2_jastrow(istate) += norm_n2_jastrow_on_grid(ipoint,istate) * weight
  enddo
 enddo
END_PROVIDER 

 BEGIN_PROVIDER [double precision, coulomb_n2_jastrow_on_grid, (n_points_final_grid,N_states)]
&BEGIN_PROVIDER [double precision, coulomb_n2_jastrow, (N_states)]
&BEGIN_PROVIDER [double precision, coulomb_n2_jastrow_renorm, (N_states)]
 implicit none
 BEGIN_DOC
! coulomb_n2_jastrow_on_grid(ipoint) = \int dr2 n2(ipoint,r2)J^2(mu(ipoint) 1/r12
!
! coulomb_n2_jastrow = \int \int dr1 dr2 n2(r1,r2)e^{J(r1,r2,\mu(r1)} 1/r12
!
! coulomb_n2_jastrow_renorm = \int dr1 \int dr2 n2(r1,r2)e^{J(r1,r2,\mu(r1)} 1/r12 / \int dr2 n2(r1,r2)e^{J(r1,r2,\mu(r1)} 
 END_DOC
 integer :: ipoint,istate,i,j,k,l
 double precision :: weight,dens,dm_a,dm_b
 coulomb_n2_jastrow_on_grid = 0.d0
 coulomb_n2_jastrow = 0.d0
 coulomb_n2_jastrow_renorm = 0.d0
 do istate =1, N_states
  do ipoint = 1, n_points_final_grid
   weight = final_weight_at_r_vector(ipoint)
   dm_a = one_e_dm_and_grad_alpha_in_r(4,ipoint,istate)
   dm_b =  one_e_dm_and_grad_beta_in_r(4,ipoint,istate)
   dens = 0.5d0 * (dm_a + dm_b)
   do i = 1, mo_num ! 1
    do k = 1, mo_num ! 1 
     do j = 1, mo_num ! 2 
      do l = 1, mo_num ! 2
       coulomb_n2_jastrow_on_grid(ipoint,istate) += mos_in_r_array(i,ipoint) * mos_in_r_array(k,ipoint) * & 
                          coulomb_jastrow2_ij_mo(l,j,ipoint) * full_occ_2_rdm_ab_chemist_mo(l,j,k,i,istate)
      enddo
     enddo
    enddo
   enddo
   coulomb_n2_jastrow(istate) += coulomb_n2_jastrow_on_grid(ipoint,istate) * weight
   coulomb_n2_jastrow_renorm(istate) += coulomb_n2_jastrow_on_grid(ipoint,istate) * dens / norm_n2_jastrow_on_grid(ipoint,istate) * weight 
  enddo
 enddo
END_PROVIDER 

 BEGIN_PROVIDER [double precision, wee_mu_of_r_n2_jastrow_on_grid, (n_points_final_grid,N_states)]
&BEGIN_PROVIDER [double precision, wee_mu_of_r_n2_jastrow, (N_states)]
&BEGIN_PROVIDER [double precision, wee_mu_of_r_n2_jastrow_renorm, (N_states)]
 implicit none
 BEGIN_DOC
! wee_mu_of_r_n2_jastrow_on_grid(ipoint) = \int dr2 n2(ipoint,r2)J^2(mu(ipoint) erf(\mu(r1) r12)/r12
!
! wee_mu_of_r_n2_jastrow = \int \int dr1 dr2 n2(r1,r2)e^{J(r1,r2,\mu(r1)} erf(\mu(r1) r12)/r12
!
! wee_mu_of_r_n2_jastrow_renorm = \int dr1 \int dr2 n2(r1,r2)e^{J(r1,r2,\mu(r1)} erf(\mu(r1) r12)/r12 / \int dr2 n2(r1,r2)e^{J(r1,r2,\mu(r1)} 
 END_DOC
 integer :: ipoint,istate,i,j,k,l
 double precision :: weight,dens,dm_a,dm_b
 wee_mu_of_r_n2_jastrow_on_grid = 0.d0
 wee_mu_of_r_n2_jastrow = 0.d0
 wee_mu_of_r_n2_jastrow_renorm = 0.d0
 do istate =1, N_states
  do ipoint = 1, n_points_final_grid
   weight = final_weight_at_r_vector(ipoint)
   dm_a = one_e_dm_and_grad_alpha_in_r(4,ipoint,istate)
   dm_b =  one_e_dm_and_grad_beta_in_r(4,ipoint,istate)
   dens = 0.5d0 * (dm_a + dm_b)
   do i = 1, mo_num ! 1
    do k = 1, mo_num ! 1 
     do j = 1, mo_num ! 2 
      do l = 1, mo_num ! 2
       wee_mu_of_r_n2_jastrow_on_grid(ipoint,istate) += mos_in_r_array(i,ipoint) * mos_in_r_array(k,ipoint) * & 
!                          erf_mu_of_r_ij_mo(l,j,ipoint) * full_occ_2_rdm_ab_chemist_mo(l,j,k,i,istate)
                          erf_mu_of_r_jastrow2_ij_mo(l,j,ipoint) * full_occ_2_rdm_ab_chemist_mo(l,j,k,i,istate)
      enddo
     enddo
    enddo
   enddo
   wee_mu_of_r_n2_jastrow(istate) += wee_mu_of_r_n2_jastrow_on_grid(ipoint,istate) * weight
   wee_mu_of_r_n2_jastrow_renorm(istate) += wee_mu_of_r_n2_jastrow_on_grid(ipoint,istate) * dens / norm_n2_jastrow_on_grid(ipoint,istate) * weight 
  enddo
 enddo
END_PROVIDER 

BEGIN_PROVIDER [double precision, psi_wee_mu_of_r, (N_states)]
 implicit none
 BEGIN_DOC
! psi_wee_mu_of_r = <Psi | W_ee^{\mu(r)} | Psi>
 END_DOC
 integer :: ipoint,istate,i,j,k,l
 double precision :: muj, muc,weight
 muj = 100000.d0
 psi_wee_mu_of_r = 0.d0
 do istate =1, N_states
  do ipoint = 1, n_points_final_grid
   weight = final_weight_at_r_vector(ipoint)
   muc = mu_of_r_prov(ipoint,istate)
   do i = 1, mo_num ! 1
    do k = 1, mo_num ! 1 
     do j = 1, mo_num ! 2 
      do l = 1, mo_num ! 2
       psi_wee_mu_of_r += weight * mos_in_r_array(i,ipoint) * mos_in_r_array(k,ipoint) * & 
                          erf_mu_of_r_ij_mo(l,j,ipoint) * full_occ_2_rdm_ab_chemist_mo(l,j,k,i,istate)
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo

END_PROVIDER 

BEGIN_PROVIDER [double precision, psi_wee_mu_of_r_sr, (N_states)]
 implicit none
 BEGIN_DOC
! psi_wee_mu_of_r_sr = <Psi| W_ee^{\mu(r),sr} | Psi>
 END_DOC
 psi_wee_mu_of_r_sr = psi_energy_two_e - psi_wee_mu_of_r

END_PROVIDER 

 BEGIN_PROVIDER [double precision, wee_mu_of_r_sr_n2_jastrow, (N_states)]
&BEGIN_PROVIDER [double precision, wee_mu_of_r_sr_n2_jastrow_renorm, (N_states)]
 implicit none
 BEGIN_DOC
! wee_mu_of_r_sr_n2_jastrow = <Psi| e^J w_ee^{\mu(r),sr} e^J |Psi>
!
! wee_mu_of_r_sr_n2_jastrow_renorm = <Psi| e^J w_ee^{\mu(r),sr} e^J |Psi>/<Psi| e^J e^J |Psi>
 END_DOC
 wee_mu_of_r_sr_n2_jastrow = coulomb_n2_jastrow - wee_mu_of_r_n2_jastrow
 wee_mu_of_r_sr_n2_jastrow_renorm = coulomb_n2_jastrow_renorm- wee_mu_of_r_n2_jastrow_renorm

END_PROVIDER 
