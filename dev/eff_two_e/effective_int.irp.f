
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
   decdrho2 = d_dn2_e_cmd_sr_pbe_n2(ipoint,istate)
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
      decdrho2 = d_dn2_e_cmd_sr_pbe_n2(ipoint,istate)
      phi_ij_eff_pot_in_r(ipoint,j,i)  = weight * decdrho2 * mos_in_r_array(i,ipoint) * mos_in_r_array(j,ipoint) 
     enddo
    enddo
   enddo
END_PROVIDER 

