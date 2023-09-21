
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
 double precision :: weight
 double precision :: decdrho2

 eff_two_e = 0.d0
 do istate = 1, N_states
  do ipoint = 1, n_points_final_grid

   weight = final_weight_at_r_vector(ipoint)
   decdrho2 = d_dn2_e_cmd_basis(ipoint,istate)
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


BEGIN_PROVIDER[double precision, tot_eff_two_e, (mo_num,mo_num,mo_num,mo_num,N_states) ]
 implicit none
 integer :: istate,ipoint,m,i,j,k,l
 double precision :: weight
 double precision :: decdrho2,get_two_e_integral,integral_regular

 tot_eff_two_e = 0.d0
 do istate = 1, N_states
   do i = 1, mo_num
    do j = 1, mo_num
     do k = 1, mo_num
      do l = 1, mo_num
        integral_regular = get_two_e_integral(i,j,k,l,mo_integrals_map)
        tot_eff_two_e(l,k,j,i,istate) =  eff_two_e(l,k,j,i,istate) + integral_regular
      enddo
     enddo
    enddo
   enddo
 enddo

END_PROVIDER

BEGIN_PROVIDER [double precision, phi_ij_eff_pot_in_r, (n_points_final_grid, mo_num, mo_num)]
 implicit none
 integer :: istate,ipoint,m,i,j,k,l
 double precision :: weight
 double precision :: decdrho2

  istate = 1
   do i = 1, mo_num
    do j = 1, mo_num
     do ipoint = 1, n_points_final_grid
      weight = final_weight_at_r_vector(ipoint)
      decdrho2 = d_dn2_e_cmd_basis(ipoint,istate)
      phi_ij_eff_pot_in_r(ipoint,j,i)  = weight * decdrho2 * mos_in_r_array(i,ipoint) * mos_in_r_array(j,ipoint) 
     enddo
    enddo
   enddo
END_PROVIDER 

BEGIN_PROVIDER[double precision, eff_two_e_lr, (mo_num,mo_num,mo_num,mo_num,N_states) ]
 implicit none
 BEGIN_DOC
 END_DOC 
 integer :: istate,ipoint,m,i,j,k,l
 double precision :: weight
 double precision :: mu, rho_a, rho_b, ec_srmuLDAn,decdrho_a,decdrho_b,d2ecdrho_a2,d2ecdrho_b2
 double precision, allocatable :: d2ecdrho2(:)
 double precision :: rho, delta_rho, delta_rho_a, delta_rho_b, rho_a_plus_delta_rho, rho_b_plus_delta_rho
 double precision :: ec_srmuLDAn_delta_rho
 double precision :: decdrho_a_delta_rho_a,decdrho_b_delta_rho_a
 double precision :: d2ecdrho_a2_delta_rho_a,d2ecdrho_b2_delta_rho_a
 double precision :: decdrho_a_delta_rho_b,decdrho_b_delta_rho_b
 double precision :: d2ecdrho_a2_delta_rho_b,d2ecdrho_b2_delta_rho_b

 eff_two_e_lr = 0.d0

 allocate(d2ecdrho2(n_points_final_grid))
 call fc_LDAUEG(d2ecdrho2) 

 do istate = 1, N_states
  do ipoint = 1, n_points_final_grid

   weight = final_weight_at_r_vector(ipoint)

   do i = 1, mo_num
    do j = 1, mo_num
     do k = 1, mo_num
      do l = 1, mo_num
        eff_two_e_lr(l,k,j,i,istate) += weight * d2ecdrho2(ipoint) * mos_in_r_array(i,ipoint) * mos_in_r_array(j,ipoint) * mos_in_r_array(k,ipoint) * mos_in_r_array(l,ipoint)
      enddo
     enddo
    enddo
   enddo


  enddo
 enddo

 deallocate(d2ecdrho2)
END_PROVIDER

