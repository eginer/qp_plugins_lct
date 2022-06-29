subroutine give_integrals_3_body_real_space(i,a,k,l,contrib)
 implicit none
 integer, intent(in) :: i,a,k,l
 double precision, intent(out) :: contrib
 double precision :: direct_int, exch_2, exch_1
 double precision :: mos_k, mos_l, mos_i, mos_a, w_ll, w_kk, w_ia
 double precision :: mos_kk, mos_ll, mos_ia, weight
 double precision :: mos_kl, mos_ki, mos_ka, mos_li, mos_la
 double precision :: w_ka,w_ki,w_la,w_li, w_kl

 integer :: mm, ipoint

 direct_int = 0.d0
 exch_1     = 0.d0
 exch_2     = 0.d0
 do mm = 1, 3
  do ipoint = 1, n_points_final_grid
   weight = final_weight_at_r_vector(ipoint)                                                                          
   ! MOs at r 
   mos_i  = mos_in_r_array_transp(ipoint,i) 
   mos_a  = mos_in_r_array_transp(ipoint,a) 
   mos_ia = mos_a * mos_i
   mos_k  = mos_in_r_array_transp(ipoint,k) 
   mos_kk = mos_k * mos_k 
   mos_l  = mos_in_r_array_transp(ipoint,l) 
   mos_ll = mos_l * mos_l
   mos_kl = mos_k * mos_l
   mos_ki = mos_k * mos_i
   mos_ka = mos_k * mos_a
   mos_la = mos_l * mos_a
   mos_li = mos_l * mos_i
   ! Integral at r
   w_ll   = x_W_ij_erf_rk(ipoint,mm,l,l) 
   w_kk   = x_W_ij_erf_rk(ipoint,mm,k,k) 
   w_ia   = x_W_ij_erf_rk(ipoint,mm,i,a) 
   w_ka   = x_W_ij_erf_rk(ipoint,mm,k,a)
   w_ki   = x_W_ij_erf_rk(ipoint,mm,k,i)
   w_kl   = x_W_ij_erf_rk(ipoint,mm,k,l)
   w_la   = x_W_ij_erf_rk(ipoint,mm,l,a)
   w_li   = x_W_ij_erf_rk(ipoint,mm,l,i)
   ! Direct integral 
   direct_int += weight * ( mos_kk * w_ll * w_ia + mos_ll * w_kk * w_ia + mos_ia * w_kk * w_ll )
   exch_1     += weight * ( mos_ki * w_ll * w_ka + mos_ll * w_ki * w_ka + mos_ka * w_ki * w_ll )
   exch_2     += weight * ( mos_kl * w_kl * w_ia + mos_kl * w_kl * w_ia + mos_ia * w_kl * w_kl ) 
  enddo
 enddo
 contrib = 1.5d0 * direct_int - exch_1 - 0.5d0 * exch_2  

end

subroutine give_fock_ia_real_space(i,a,contrib)
 implicit none
 integer, intent(in) :: i,a
 double precision, intent(out) :: contrib
 double precision :: direct_int, exch_2, exch_1
 double precision :: mos_k, mos_l, mos_i, mos_a, w_ll, w_kk, w_ia
 double precision :: mos_kk, mos_ll, mos_ia, weight
 double precision :: mos_kl, mos_ki, mos_ka, mos_li, mos_la
 double precision :: w_ka,w_ki,w_la,w_li, w_kl

 integer :: mm, ipoint,k,l

 direct_int = 0.d0
 exch_1     = 0.d0
 exch_2     = 0.d0
 do mm = 1, 3
  do ipoint = 1, n_points_final_grid
   do k = 1, elec_beta_num
    do l = 1, elec_beta_num
     weight = final_weight_at_r_vector(ipoint)                                                                          
     ! MOs at r 
     mos_i  = mos_in_r_array_transp(ipoint,i) 
     mos_a  = mos_in_r_array_transp(ipoint,a) 
     mos_ia = mos_a * mos_i
     mos_k  = mos_in_r_array_transp(ipoint,k) 
     mos_kk = mos_k * mos_k 
     mos_l  = mos_in_r_array_transp(ipoint,l) 
     mos_ll = mos_l * mos_l
     mos_kl = mos_k * mos_l
     mos_ki = mos_k * mos_i
     mos_ka = mos_k * mos_a
     mos_la = mos_l * mos_a
     mos_li = mos_l * mos_i
     ! Integral at r
     w_ll   = x_W_ij_erf_rk(ipoint,mm,l,l) 
     w_kk   = x_W_ij_erf_rk(ipoint,mm,k,k) 
     w_ia   = x_W_ij_erf_rk(ipoint,mm,i,a) 
     w_ka   = x_W_ij_erf_rk(ipoint,mm,k,a)
     w_ki   = x_W_ij_erf_rk(ipoint,mm,k,i)
     w_kl   = x_W_ij_erf_rk(ipoint,mm,k,l)
     w_la   = x_W_ij_erf_rk(ipoint,mm,l,a)
     w_li   = x_W_ij_erf_rk(ipoint,mm,l,i)
     ! Direct integral 
     direct_int += weight * ( mos_kk * w_ll * w_ia + mos_ll * w_kk * w_ia + mos_ia * w_kk * w_ll )
     exch_1     += weight * ( mos_ki * w_ll * w_ka + mos_ll * w_ki * w_ka + mos_ka * w_ki * w_ll )
     exch_2     += weight * ( mos_kl * w_kl * w_ia + mos_kl * w_kl * w_ia + mos_ia * w_kl * w_kl ) 
    enddo
   enddo
  enddo
 enddo
 contrib = 1.5d0 * direct_int - exch_1 - 0.5d0 * exch_2  

end

subroutine give_fock_ia_real_space_old(i,a,contrib)
 implicit none
 integer, intent(in) :: i,a
 double precision, intent(out) :: contrib
 double precision :: int_1, int_2, int_3, int_4
 double precision :: mos_k, mos_l, mos_i, mos_a, w_ll, w_kk, w_ia
 double precision :: mos_kk, mos_ll, mos_ia, weight
 double precision :: mos_kl, mos_ki, mos_ka, mos_li, mos_la
 double precision :: w_ka,w_ki,w_la,w_li, w_kl

 integer :: mm, ipoint,k,l

 int_1 = 0.d0
 int_2 = 0.d0
 int_3 = 0.d0
 int_4 = 0.d0
 do mm = 1, 3
  do ipoint = 1, n_points_final_grid
   do k = 1, elec_beta_num
    do l = 1, elec_beta_num
     weight = final_weight_at_r_vector(ipoint)                                                                          
     ! MOs at r 
     mos_i  = mos_in_r_array_transp(ipoint,i) 
     mos_a  = mos_in_r_array_transp(ipoint,a) 
     mos_ia = mos_a * mos_i
     mos_k  = mos_in_r_array_transp(ipoint,k) 
     mos_kk = mos_k * mos_k 
     mos_l  = mos_in_r_array_transp(ipoint,l) 
     mos_ll = mos_l * mos_l
     mos_kl = mos_k * mos_l
     mos_ki = mos_k * mos_i
     mos_ka = mos_k * mos_a
     mos_la = mos_l * mos_a
     mos_li = mos_l * mos_i
     ! Integral at r
     w_ll   = x_W_ij_erf_rk(ipoint,mm,l,l) 
     w_kk   = x_W_ij_erf_rk(ipoint,mm,k,k) 
     w_ia   = x_W_ij_erf_rk(ipoint,mm,i,a) 
     w_ka   = x_W_ij_erf_rk(ipoint,mm,k,a)
     w_ki   = x_W_ij_erf_rk(ipoint,mm,k,i)
     w_kl   = x_W_ij_erf_rk(ipoint,mm,k,l)
     w_la   = x_W_ij_erf_rk(ipoint,mm,l,a)
     w_li   = x_W_ij_erf_rk(ipoint,mm,l,i)
     ! Direct integral 
     int_1  += weight * w_kk * 1.5d0 * (mos_ll * w_ia + mos_ia * w_ll) 
     int_2  += weight * w_ll * (1.5d0  * mos_kk * w_ia - mos_ka * w_ki - mos_ki * w_ka )
     int_3  += weight * (-1.d0  * w_ia * mos_kl * w_kl -1.d0 * mos_ll * w_ki * w_ka )
     int_4  += weight * (-0.5d0 * mos_ia * w_kl * w_kl )
    enddo
   enddo
  enddo
 enddo
 contrib = int_1 + int_2 + int_3 + int_4 

end

subroutine give_fock_ia_real_space_old_bis(i,a,contrib)
 implicit none
 integer, intent(in) :: i,a
 double precision, intent(out) :: contrib
 double precision :: int_1, int_2, int_3, int_4
 double precision :: mos_k, mos_l, mos_i, mos_a, w_ll, w_kk, w_ia
 double precision :: mos_kk, mos_ll, mos_ia, weight
 double precision :: mos_kl, mos_ki, mos_ka, mos_li, mos_la
 double precision :: w_ka,w_ki,w_la,w_li, w_kl

 integer :: mm, ipoint,k,l

 int_1 = 0.d0
 int_2 = 0.d0
 int_3 = 0.d0
 int_4 = 0.d0
 do mm = 1, 3
  do ipoint = 1, n_points_final_grid
   weight = final_weight_at_r_vector(ipoint)                                                                          
   ! MOs at r 
   mos_i  = mos_in_r_array_transp(ipoint,i) 
   mos_a  = mos_in_r_array_transp(ipoint,a) 
   mos_ia = mos_a * mos_i
   w_ia   = x_W_ij_erf_rk(ipoint,mm,i,a) 
   double precision :: w_kk_sum,psi_k_i,psi_k_a,w_ki_mos_k,w_ka_mos_k
   w_kk_sum = 0.d0
   psi_k_a = 0.d0
   psi_k_i = 0.d0
   w_ki_mos_k = 0.d0
   w_ka_mos_k = 0.d0
   do k = 1, elec_beta_num
     mos_k  = mos_in_r_array_transp(ipoint,k) 
     mos_kk = mos_k * mos_k 
     w_kk   = x_W_ij_erf_rk(ipoint,mm,k,k) 
     psi_k_i += mos_k * mos_i
     psi_k_a += mos_k * mos_a
     w_kk_sum += w_kk
     w_ka   = x_W_ij_erf_rk(ipoint,mm,k,a)
     w_ki   = x_W_ij_erf_rk(ipoint,mm,k,i)
     w_ki_mos_k += w_ki * mos_k 
     w_ka_mos_k += w_ka * mos_k 
   enddo
   double precision :: rho_beta
   rho_beta = 0.d0
   do l = 1, elec_beta_num
    mos_l  = mos_in_r_array_transp(ipoint,l) 
    mos_ll = mos_l * mos_l
    rho_beta += mos_ll
   enddo
!    w_ll   = x_W_ij_erf_rk(ipoint,mm,l,l) 
    int_1  += weight * w_kk_sum * 1.5d0 * (rho_beta * w_ia + mos_ia * w_kk_sum)
    int_2  += weight * w_kk_sum * ( 1.5d0 *  rho_beta * w_ia - w_ki_mos_k * mos_a - w_ka_mos_k * mos_i)
  enddo
 enddo
 contrib = int_1 + int_2 + int_3 + int_4 

end

subroutine give_fock_ia_real_space_new(i,a,contrib)
 implicit none
 integer, intent(in) :: i,a
 double precision, intent(out) :: contrib
 double precision :: int_1, int_2, int_3, int_4
 double precision :: mos_k, mos_l, mos_i, mos_a, w_ll, w_kk, w_ia
 double precision :: mos_kk, mos_ll, mos_ia, weight
 double precision :: mos_kl, mos_ki, mos_ka, mos_li, mos_la
 double precision :: w_ka,w_ki,w_la,w_li, w_kl

 integer :: mm, ipoint,k,l

 int_1 = 0.d0
 int_2 = 0.d0
 int_3 = 0.d0
 int_4 = 0.d0
 do mm = 1, 3
  do ipoint = 1, n_points_final_grid
   weight = final_weight_at_r_vector(ipoint)                                                                          
   mos_i  = mos_in_r_array_transp(ipoint,i) 
   mos_a  = mos_in_r_array_transp(ipoint,a) 
   mos_ia = mos_a * mos_i
   w_ia   = x_W_ij_erf_rk(ipoint,mm,i,a) 
   

   integer :: m 
   double precision :: big_v_r,w_ia_tilde(mo_num, mo_num)
   w_ia_tilde = 0.d0
   do l = 1, mo_num
    do m = 1, mo_num
     do k = 1, elec_beta_num
      w_ka   = x_W_ij_erf_rk(ipoint,mm,k,m)
      w_ki   = x_W_ij_erf_rk(ipoint,mm,k,l)
      w_ia_tilde(m,l) += w_ki * w_ka 
     enddo
    enddo
   enddo
   big_v_r = 0.d0
   do l = 1, elec_beta_num
    big_v_r += w_ia_tilde(l,l)
   enddo
   
   double precision :: v_r
   v_r = 0.d0
   do k = 1, elec_beta_num
    mos_k  = mos_in_r_array_transp(ipoint,k) 
    do l = 1, elec_beta_num
     mos_l  = mos_in_r_array_transp(ipoint,l) 
     w_kl   = x_W_ij_erf_rk(ipoint,mm,l,k)
     v_r += w_kl * mos_k * mos_l 
    enddo
   enddo

   rho_k    = 0.d0
   v_i      = 0.d0
   big_w_kk = 0.d0
   big_v_r  = 0.d0
   w_ki_mos_k = 0.d0
   w_ka_mos_k = 0.d0
   do k = 1, elec_beta_num
    double precision :: big_w_kk,rho_k,v_i
    double precision :: w_ki_mos_k,w_ka_mos_k
    mos_k  = mos_in_r_array_transp(ipoint,k) 
    mos_kk = mos_k * mos_k 
    mos_ki = mos_k * mos_i
    mos_ka = mos_k * mos_a
    w_kk   = x_W_ij_erf_rk(ipoint,mm,k,k) 
    w_ka   = x_W_ij_erf_rk(ipoint,mm,k,a)
    w_ki   = x_W_ij_erf_rk(ipoint,mm,k,i)
    rho_k    += mos_kk 
    v_i      += mos_k * w_ki 
    big_w_kk += w_kk
    big_v_r  += w_ia_tilde(k,k)
    w_ki_mos_k += w_ki * mos_k 
    w_ka_mos_k += w_ka * mos_k 
   enddo
     
   int_2  += weight * big_w_kk * (3.d0 * rho_k * w_ia + 1.5d0 * mos_ia * big_w_kk - w_ki_mos_k * mos_a - w_ka_mos_k * mos_i )
   int_3  += weight * (-1.d0) * ( v_r * w_ia + rho_k * w_ia_tilde(i,a)  + 0.5d0 * mos_ia * big_v_r)
  enddo
 enddo
 contrib = int_1 + int_2 + int_3 

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

 diag_three_elem_hf_old = - diag_three_elem_hf_old 
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


