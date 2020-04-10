
BEGIN_PROVIDER[double precision, eff_two_e, (mo_num,mo_num,mo_num,mo_num,N_states) ]
 implicit none
 BEGIN_DOC
 !                = \int dr de_cmd/dn2(r) \hat{n}_2(r,r)
 !
 ! \hat{n}_2(r,r) = \sum_{i,j,k,l} a^{\dagger}_k  a^{\dagger}_l  a_j a_i \phi_k(r)  \phi_l(r)  \phi_j(r)  \phi_i(r) 
 !
 ! \sum_{i,j,k,l} a^{\dagger}_k  a^{\dagger}_l  a_j a_i   \int dr de_cmd/dn2(r) * \phi_k(r)  \phi_l(r)  \phi_j(r)  \phi_i(r) 
 !
 ! \sum_{i,j,k,l} a^{\dagger}_k  a^{\dagger}_l  a_j a_i  <ij|kl>_{eff}
 !
 ! <ij|kl>_{eff} =  \int dr de_cmd/dn2(r) * \phi_k(r)  \phi_l(r)  \phi_j(r)  \phi_i(r)
 END_DOC 
 integer :: istate,ipoint,m,i,j,k,l
 double precision :: two_dm_in_r_exact
 double precision :: weight, r(3)
 double precision :: ec_srmuPBE,mu
 double precision :: rho2, rho_a,rho_b,grad_rho_a(3),grad_rho_b(3),grad_rho_a_2,grad_rho_b_2,grad_rho_a_b
 double precision :: decdrho_a, decdrho_b, decdrho, decdrho2
 double precision :: decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b, decdgrad_rho_2

 eff_two_e = 0.d0
 do istate = 1, N_states
  do ipoint = 1, n_points_final_grid
   r(1) = final_grid_points(1,ipoint)
   r(2) = final_grid_points(2,ipoint)
   r(3) = final_grid_points(3,ipoint)

   weight = final_weight_at_r_vector(ipoint)
   rho_a =  one_e_dm_and_grad_alpha_in_r(4,ipoint,istate)
   rho_b =  one_e_dm_and_grad_beta_in_r(4,ipoint,istate)
   !rho2 = two_dm_in_r_exacti(r,r,istate)
   call give_on_top_in_r_one_state_local(r,istate,rho2)
   grad_rho_a(1:3) =  one_e_dm_and_grad_alpha_in_r(1:3,ipoint,istate)
   grad_rho_b(1:3) =  one_e_dm_and_grad_beta_in_r(1:3,ipoint,istate)
   grad_rho_a_2 = 0.d0
   grad_rho_b_2 = 0.d0
   grad_rho_a_b = 0.d0
   do m = 1, 3
    grad_rho_a_2 += grad_rho_a(m) * grad_rho_a(m)
    grad_rho_b_2 += grad_rho_b(m) * grad_rho_b(m)
    grad_rho_a_b += grad_rho_a(m) * grad_rho_b(m)
   enddo

   rho2 = rho2*2.d0 ! normalization 
   ! mu = mu_b(r)
   mu = mu_of_r_prov(ipoint,istate)
   call ecmdsrPBEn2(mu,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,rho2,ec_srmuPBE,decdrho_a,decdrho_b, decdrho, decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b, decdgrad_rho_2,decdrho2)
   decdrho2 = 2.d0*decdrho2 ! normalization
   do i = 1, mo_num
    do j = 1, mo_num
     do k = 1, mo_num
      do l = 1, mo_num
        eff_two_e(l,k,j,i,istate) += weight * decdrho2 * mos_in_r_array(i,ipoint) * mos_in_r_array(j,ipoint) * mos_in_r_array(k,ipoint) * mos_in_r_array(l,ipoint)
      enddo
     enddo
    enddo
   enddo

  enddo
 enddo

END_PROVIDER


