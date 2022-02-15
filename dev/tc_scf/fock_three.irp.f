BEGIN_PROVIDER [ double precision, fock_3_mat, (mo_num, mo_num)] 
 implicit none
  integer :: i,j
  double precision :: contrib
  call give_fock_ia_real_space_prov(1,1,contrib)
  fock_3_mat = 0.d0
!  !$OMP PARALLEL                  &
!  !$OMP DEFAULT (NONE)            &
!  !$OMP PRIVATE (i,j,m,integral) & 
!  !$OMP SHARED (mo_num,three_body_3_index)
!  !$OMP DO SCHEDULE (guided) COLLAPSE(3)
  do i = 1, mo_num
   do j = 1, mo_num
    call give_fock_ia_three_e_total(j,i,contrib)
    fock_3_mat(j,i) = -contrib
   enddo
  enddo
!  !$OMP END DO
!  !$OMP END PARALLEL
!  do i = 1, mo_num
!   do j = 1, i-1
!    mat_three(j,i) = mat_three(i,j)
!   enddo
!  enddo

END_PROVIDER 


subroutine give_fock_ia_three_e_total(i,a,contrib)
 implicit none
 BEGIN_DOC
! contrib is the TOTAL (same spins / opposite spins) contribution from the three body term to the Fock operator 
!
 END_DOC
 integer, intent(in) :: i,a
 double precision, intent(out) :: contrib
 double precision :: int_1, int_2, int_3
 double precision :: mos_i, mos_a, w_ia
 double precision :: mos_ia, weight

 integer :: mm, ipoint,k,l

 int_1 = 0.d0
 int_2 = 0.d0
 int_3 = 0.d0
 do mm = 1, 3
  do ipoint = 1, n_points_final_grid
   weight = final_weight_at_r_vector(ipoint)                                                                          
   mos_i  = mos_in_r_array_transp(ipoint,i) 
   mos_a  = mos_in_r_array_transp(ipoint,a) 
   mos_ia = mos_a * mos_i
   w_ia   = x_W_ij_erf_rk(ipoint,mm,i,a) 
     
   int_1  += weight * fock_3_w_kk_sum(ipoint,mm) * (4.d0 * fock_3_rho_beta(ipoint) * w_ia               & 
                                                  + 2.0d0 * mos_ia * fock_3_w_kk_sum(ipoint,mm)         & 
                                                  - 2.0d0 * fock_3_w_ki_mos_k(ipoint,mm,i) * mos_a      & 
                                                  - 2.0d0 * fock_3_w_ki_mos_k(ipoint,mm,a) * mos_i      )
   int_2  += weight * (-1.d0) * ( 2.0d0 * fock_3_w_kl_mo_k_mo_l(ipoint,mm) * w_ia                     & 
                                + 2.0d0 * fock_3_rho_beta(ipoint) * fock_3_w_ki_wk_a(ipoint,mm,i,a)   & 
                                + 1.0d0 * mos_ia * fock_3_trace_w_tilde(ipoint,mm)                    )

   int_3  += weight *   1.d0  * (fock_3_w_kl_wla_phi_k(ipoint,mm,i) * mos_a + fock_3_w_kl_wla_phi_k(ipoint,mm,a) * mos_i & 
                                +fock_3_w_ki_mos_k(ipoint,mm,i)     * fock_3_w_ki_mos_k(ipoint,mm,a)                     )
  enddo
 enddo
 contrib = int_1 + int_2 + int_3

end

BEGIN_PROVIDER [double precision, diag_three_elem_hf]
 implicit none
 integer :: i,j,k,ipoint,mm
 double precision :: contrib,weight,four_third,one_third,two_third,exchange_int_231
 one_third = 1.d0/3.d0
 two_third = 2.d0/3.d0
 four_third = 4.d0/3.d0
 diag_three_elem_hf = 0.d0
 do i = 1, elec_beta_num
  do j = 1, elec_beta_num
   do k = 1, elec_beta_num
    call  give_integrals_3_body(k,j,i,j,i,k,exchange_int_231)   
    diag_three_elem_hf += two_third * exchange_int_231
   enddo
  enddo
 enddo
 do mm = 1, 3
  do ipoint = 1, n_points_final_grid
   weight = final_weight_at_r_vector(ipoint)                                                                          
   contrib   = 3.d0 * fock_3_w_kk_sum(ipoint,mm) * fock_3_rho_beta(ipoint) * fock_3_w_kk_sum(ipoint,mm)  & 
              -2.d0 * fock_3_w_kl_mo_k_mo_l(ipoint,mm) * fock_3_w_kk_sum(ipoint,mm)                                 & 
              -1.d0 * fock_3_rho_beta(ipoint) * fock_3_w_kl_w_kl(ipoint,mm)
   contrib  *= four_third
   contrib  += -two_third  * fock_3_rho_beta(ipoint)     * fock_3_w_kl_w_kl(ipoint,mm) & 
              - four_third * fock_3_w_kk_sum(ipoint,mm)  * fock_3_w_kl_mo_k_mo_l(ipoint,mm)
   diag_three_elem_hf += weight * contrib
  enddo
 enddo
 diag_three_elem_hf = - diag_three_elem_hf
END_PROVIDER 

