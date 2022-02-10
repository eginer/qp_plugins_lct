BEGIN_PROVIDER [ double precision, ref_fock_three, (mo_num, mo_num)]
 implicit none
 integer :: i,a,k,l
 double precision :: direct_int, exch_1, exch_2
 ref_fock_three = 0.d0
 do i = 1, mo_num
  do a = 1, mo_num
   do k = 1, elec_beta_num
    do l = k+1, elec_beta_num
     call  give_integrals_3_body(k,l,i,k,l,a,direct_int)   
     call  give_integrals_3_body(l,k,i,k,l,a,exch_1)   
     ref_fock_three(i,a) += 1.d0 * direct_int - exch_1 
    enddo
    do l = 1, elec_beta_num
     call  give_integrals_3_body(k,l,i,k,l,a,direct_int)   
     call  give_integrals_3_body(i,l,k,k,l,a,exch_1)   
     ref_fock_three(i,a) += 1.d0 * direct_int - exch_1 
    enddo
   enddo
  enddo
 enddo
END_PROVIDER 

BEGIN_PROVIDER [ double precision, ref_fock_three_new, (mo_num, mo_num)]
 implicit none
 integer :: i,a,k,l
 double precision :: direct_int, contrib
 ref_fock_three_new = 0.d0
 do i = 1, mo_num
  do a = 1, mo_num
   call give_contrib_three_fock(i,a,contrib)
   ref_fock_three_new(i,a) += contrib 
  enddo
 enddo
END_PROVIDER 


subroutine give_contrib_three_fock(i,a,contrib)
 implicit none
 integer, intent(in) :: i,a
 double precision, intent(out) :: contrib
 double precision :: direct_int, exch_1, exch_2 
 integer :: k,l
 contrib = 0.d0
 h = i
 p = a
 do k = 1, elec_beta_num
  do l = 1, elec_beta_num
   call  give_integrals_3_body(h,k,l,p,k,l,direct_int)  ! < h k l | p k l >
   call  give_integrals_3_body(h,k,l,k,p,l,exch_1)      ! < h k l | k p l >
   call  give_integrals_3_body(h,k,l,p,l,k,exch_2)      ! < h k l | p l k >
   !   anti parallel spin : 1.5 direct - 1 exch_1 - 0.5 exch_2
   contrib += 2.0d0 * direct_int - 1.5d0 * exch_1 - 1.0d0 * exch_2  
!  !!contrib += 1.5d0 * direct_int - 1.0d0 * exch_1 - 0.5d0 * exch_2  
  enddo
 enddo
 double precision :: exchange_int_13,exchange_int_12, exchange_int_23
 double precision :: exchange_int_123, exchange_int_321,same_spin
 integer :: h,p
 same_spin = 0.d0
 do k = 1, elec_beta_num
  do l = 1, elec_beta_num
!   call  give_integrals_3_body(h,k,l,p,k,l,direct_int)         ! <h k l | p k l >
!   call  give_integrals_3_body(h,k,l,p,l,k,exchange_int_23)    ! <h k l | p l k >
!   call  give_integrals_3_body(h,k,l,k,p,l,exchange_int_12)    ! <h k l | k p l >
   call  give_integrals_3_body(h,k,l,k,l,p,exchange_int_123)   ! <h k l | k l p > 
   call  give_integrals_3_body(h,k,l,l,p,k,exchange_int_321)   ! <h k l | l p k > 
   call  give_integrals_3_body(h,k,l,l,k,p,exchange_int_13)    ! <h k l | l k p >
   same_spin += exchange_int_123 + exchange_int_321 - exchange_int_13
            
  enddo
 enddo
 contrib += 0.5d0 * same_spin
 
end


BEGIN_PROVIDER [ double precision, fock_3_w_kk_sum, (n_points_final_grid,3)]
 implicit none
 integer :: mm, ipoint,k
 double precision :: w_kk
 fock_3_w_kk_sum = 0.d0
 do k = 1, elec_beta_num
  do mm = 1, 3
   do ipoint = 1, n_points_final_grid
    w_kk   = x_W_ij_erf_rk(ipoint,mm,k,k) 
    fock_3_w_kk_sum(ipoint,mm) += w_kk
   enddo
  enddo
 enddo
END_PROVIDER 

BEGIN_PROVIDER [ double precision, fock_3_w_ki_mos_k, (n_points_final_grid,3,mo_num)]
 implicit none
 integer :: mm, ipoint,k,i
 double precision :: w_ki, mo_k
 fock_3_w_ki_mos_k = 0.d0
 do i = 1, mo_num
  do k = 1, elec_beta_num
   do mm = 1, 3
    do ipoint = 1, n_points_final_grid
     w_ki   = x_W_ij_erf_rk(ipoint,mm,k,i) 
     mo_k = mos_in_r_array(k,ipoint)
     fock_3_w_ki_mos_k(ipoint,mm,i) += w_ki * mo_k
    enddo
   enddo
  enddo
 enddo

END_PROVIDER 

BEGIN_PROVIDER [ double precision, fock_3_w_kl_w_kl, (n_points_final_grid,3)]
 implicit none
 integer :: k,j,ipoint,mm
 double precision :: w_kj
 fock_3_w_kl_w_kl = 0.d0
 do j = 1, elec_beta_num
  do k = 1, elec_beta_num
   do mm = 1, 3
    do ipoint = 1, n_points_final_grid
     w_kj   = x_W_ij_erf_rk(ipoint,mm,k,j) 
     fock_3_w_kl_w_kl(ipoint,mm) += w_kj * w_kj
    enddo
   enddo
  enddo
 enddo


END_PROVIDER 

BEGIN_PROVIDER [ double precision, fock_3_rho_beta, (n_points_final_grid)]
 implicit none
 integer :: ipoint,k
 fock_3_rho_beta = 0.d0
 do ipoint = 1, n_points_final_grid
  do k = 1, elec_beta_num
   fock_3_rho_beta(ipoint) += mos_in_r_array(k,ipoint) * mos_in_r_array(k,ipoint)
  enddo
 enddo
END_PROVIDER 

BEGIN_PROVIDER [ double precision, fock_3_w_kl_mo_k_mo_l, (n_points_final_grid,3)]
 implicit none
 integer :: ipoint,k,l,mm
 double precision :: mos_k, mos_l, w_kl
 fock_3_w_kl_mo_k_mo_l = 0.d0
 do k = 1, elec_beta_num
  do l = 1, elec_beta_num
   do mm = 1, 3
    do ipoint = 1, n_points_final_grid
     mos_k  = mos_in_r_array_transp(ipoint,k) 
     mos_l  = mos_in_r_array_transp(ipoint,l) 
     w_kl   = x_W_ij_erf_rk(ipoint,mm,l,k)
     fock_3_w_kl_mo_k_mo_l(ipoint,mm) += w_kl * mos_k * mos_l 
    enddo
   enddo
  enddo
 enddo

END_PROVIDER 

BEGIN_PROVIDER [ double precision, fock_3_w_ki_wk_a, (n_points_final_grid,3,mo_num, mo_num)]
 implicit none
 integer :: ipoint,i,a,k,mm
 double precision :: w_ki,w_ka
 fock_3_w_ki_wk_a = 0.d0
 do i = 1, mo_num
  do a = 1, mo_num
   do mm = 1, 3
    do ipoint = 1, n_points_final_grid
     do k = 1, elec_beta_num
      w_ki   = x_W_ij_erf_rk(ipoint,mm,k,i)
      w_ka   = x_W_ij_erf_rk(ipoint,mm,k,a)
      fock_3_w_ki_wk_a(ipoint,mm,a,i) += w_ki * w_ka
     enddo
    enddo
   enddo
  enddo
 enddo
END_PROVIDER 

BEGIN_PROVIDER [ double precision, fock_3_trace_w_tilde, (n_points_final_grid,3)]
 implicit none
 integer :: ipoint,k,mm
 fock_3_trace_w_tilde = 0.d0
 do k = 1, elec_beta_num
   do mm = 1, 3
    do ipoint = 1, n_points_final_grid
     fock_3_trace_w_tilde(ipoint,mm) += fock_3_w_ki_wk_a(ipoint,mm,k,k)
    enddo
   enddo
 enddo

END_PROVIDER 

BEGIN_PROVIDER [ double precision, fock_3_w_kl_wla_phi_k, (n_points_final_grid,3,mo_num)]
 implicit none
 integer :: ipoint,a,k,mm,l
 double precision :: w_kl,w_la, mo_k
 fock_3_w_kl_wla_phi_k = 0.d0
 do a = 1, mo_num
  do k = 1, elec_beta_num 
   do l = 1, elec_beta_num
    do mm = 1, 3
     do ipoint = 1, n_points_final_grid
      w_kl   = x_W_ij_erf_rk(ipoint,mm,l,k)
      w_la   = x_W_ij_erf_rk(ipoint,mm,l,a)
      mo_k  = mos_in_r_array_transp(ipoint,k) 
      fock_3_w_kl_wla_phi_k(ipoint,mm,a) += w_kl * w_la * mo_k
     enddo
    enddo
   enddo
  enddo
 enddo
END_PROVIDER 

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
    call give_fock_ia_real_space_prov(j,i,contrib)
    fock_3_mat(j,i) = contrib
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

subroutine give_fock_ia_real_space_prov(i,a,contrib)
 implicit none
 BEGIN_DOC
! contrib is the contribution from the three body term to the Fock operator 
!
! WARNING : It takes into account only the closed shell part, and neglect the purely parallel spins contributions.
 END_DOC
 integer, intent(in) :: i,a
 double precision, intent(out) :: contrib
 double precision :: big_v_r
 double precision :: int_1, int_2
 double precision :: mos_i, mos_a, w_ia
 double precision :: mos_ia, weight

 integer :: mm, ipoint,k,l

 int_1 = 0.d0
 int_2 = 0.d0
 do mm = 1, 3
  do ipoint = 1, n_points_final_grid
   weight = final_weight_at_r_vector(ipoint)                                                                          
   mos_i  = mos_in_r_array_transp(ipoint,i) 
   mos_a  = mos_in_r_array_transp(ipoint,a) 
   mos_ia = mos_a * mos_i
   w_ia   = x_W_ij_erf_rk(ipoint,mm,i,a) 
   big_v_r  = fock_3_trace_w_tilde(ipoint,mm)
     
   int_1  += weight * fock_3_w_kk_sum(ipoint,mm) * (3.d0 * fock_3_rho_beta(ipoint) * w_ia       & 
                                                  + 1.5d0 * mos_ia * fock_3_w_kk_sum(ipoint,mm) & 
                                                  - 1.0d0 * fock_3_w_ki_mos_k(ipoint,mm,i) * mos_a      & 
                                                  - 1.0d0 * fock_3_w_ki_mos_k(ipoint,mm,a) * mos_i )
   int_2  += weight * (-1.d0) * ( 1.0d0 * fock_3_w_kl_mo_k_mo_l(ipoint,mm) * w_ia + 1.0d0 * fock_3_rho_beta(ipoint) * fock_3_w_ki_wk_a(ipoint,mm,i,a)  + 0.5d0 * mos_ia * big_v_r)
  enddo
 enddo
 contrib = int_1 + int_2 

end

subroutine give_fock_ia_scaled_op_spin(i,a,contrib)
 implicit none
 BEGIN_DOC
! contrib is the contribution from the three body term to the Fock operator 
!
! WARNING : It takes into the scaled anti-parrallel part to contain a part of direct and exchange of parrallel
 END_DOC
 integer, intent(in) :: i,a
 double precision, intent(out) :: contrib
 double precision :: big_v_r
 double precision :: int_1, int_2
 double precision :: mos_i, mos_a, w_ia
 double precision :: mos_ia, weight

 integer :: mm, ipoint,k,l

 int_1 = 0.d0
 int_2 = 0.d0
 do mm = 1, 3
  do ipoint = 1, n_points_final_grid
   weight = final_weight_at_r_vector(ipoint)                                                                          
   mos_i  = mos_in_r_array_transp(ipoint,i) 
   mos_a  = mos_in_r_array_transp(ipoint,a) 
   mos_ia = mos_a * mos_i
   w_ia   = x_W_ij_erf_rk(ipoint,mm,i,a) 
   big_v_r  = fock_3_trace_w_tilde(ipoint,mm)
     
   int_1  += weight * fock_3_w_kk_sum(ipoint,mm) * (4.d0 * fock_3_rho_beta(ipoint) * w_ia       & 
                                                  + 2.0d0 * mos_ia * fock_3_w_kk_sum(ipoint,mm) & 
                                                  - 1.5d0 * fock_3_w_ki_mos_k(ipoint,mm,i) * mos_a      & 
                                                  - 1.5d0 * fock_3_w_ki_mos_k(ipoint,mm,a) * mos_i )
   int_2  += weight * (-1.d0) * ( 2.0d0 * fock_3_w_kl_mo_k_mo_l(ipoint,mm) * w_ia + 1.5d0 * fock_3_rho_beta(ipoint) * fock_3_w_ki_wk_a(ipoint,mm,i,a)  + 1.0d0 * mos_ia * big_v_r)
  enddo
 enddo
 contrib = int_1 + int_2 

end


subroutine give_fock_ia_same_spin(i,a,contrib)
 implicit none
 BEGIN_DOC
! contrib is the contribution from the three body term to the Fock operator 
!
! WARNING : It takes into parrallel part of the three body term
 END_DOC
 integer, intent(in) :: i,a
 double precision, intent(out) :: contrib
 double precision :: big_v_r
 double precision :: int_1, int_2
 double precision :: mos_i, mos_a, w_ia
 double precision :: mos_ia, weight

 integer :: mm, ipoint,k,l

 int_1 = 0.d0
 int_2 = 0.d0
 do mm = 1, 3
  do ipoint = 1, n_points_final_grid
   weight = final_weight_at_r_vector(ipoint)                                                                          
   mos_i  = mos_in_r_array_transp(ipoint,i) 
   mos_a  = mos_in_r_array_transp(ipoint,a) 
   mos_ia = mos_a * mos_i
   w_ia   = x_W_ij_erf_rk(ipoint,mm,i,a) 
     
   int_1  += weight *   2.d0  * (fock_3_w_kl_wla_phi_k(ipoint,mm,i) * mos_a + fock_3_w_kl_wla_phi_k(ipoint,mm,a) * mos_i & 
             +fock_3_w_ki_mos_k(ipoint,mm,i) * fock_3_w_ki_mos_k(ipoint,mm,a))                          
   int_2  += weight *  ( fock_3_rho_beta(ipoint) * fock_3_w_ki_wk_a(ipoint,mm,i,a) & 
             + fock_3_w_kk_sum(ipoint,mm) * ( fock_3_w_ki_mos_k(ipoint,mm,a) * mos_i  + fock_3_w_ki_mos_k(ipoint,mm,i) * mos_a) )
  enddo
 enddo
 contrib = 0.5d0 * (int_1 - int_2 )

end

BEGIN_PROVIDER [double precision, diag_three_elem_hf_old]
 implicit none
 double precision :: direct_int,exch_int,four_third,one_third,two_third
 double precision :: exchange_int_231 , exchange_int_312 , exchange_int_23 , exchange_int_12 , exchange_int_13
 integer :: i,j,k
 one_third = 1.d0/3.d0
 two_third = 2.d0/3.d0
 four_third = 4.d0/3.d0
 diag_three_elem_hf_old = 0.d0

 do i = 1, elec_beta_num
  do j = 1, elec_beta_num
   do k = 1, elec_beta_num
    call  give_integrals_3_body(i,j,k,i,j,k,direct_int)   
    call  give_integrals_3_body(k,j,i,k,i,j,exchange_int_23)   
    diag_three_elem_hf_old += four_third * (direct_int - exchange_int_23) 
    call  give_integrals_3_body(k,j,i,j,i,k,exchange_int_231)   
    call  give_integrals_3_body(k,j,i,j,k,i,exchange_int_12)   
    call  give_integrals_3_body(k,j,i,i,j,k,exchange_int_13)   
    diag_three_elem_hf_old += two_third *  exchange_int_231 & 
                            - one_third * (exchange_int_12 + exchange_int_13)
                            
   enddo
  enddo
 enddo
 double precision :: same_spin
 same_spin = 0.d0
! do i = 1, elec_beta_num
!  do j = 1, elec_beta_num
!   do k = 1, elec_beta_num
!    call  give_integrals_3_body(k,j,i,k,j,i,direct_int)   
!    call  give_integrals_3_body(k,j,i,j,i,k,exchange_int_231)   
!    call  give_integrals_3_body(k,j,i,j,k,i,exchange_int_12)   
!    call  give_integrals_3_body(k,j,i,i,j,k,exchange_int_13)   
!    call  give_integrals_3_body(k,j,i,k,i,j,exchange_int_23)   
!    same_spin  += direct_int + 2.d0 * exchange_int_231 & 
!                - exchange_int_12 - exchange_int_13 - exchange_int_23
!   enddo
!  enddo
! enddo
 same_spin *= 1.d0/3.d0

 diag_three_elem_hf_old = - diag_three_elem_hf_old - same_spin
END_PROVIDER 

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

