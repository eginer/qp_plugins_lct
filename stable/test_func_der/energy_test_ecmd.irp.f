subroutine energy_cmd_supbeot_test (rho_a, rho_b, grad_rho_a, grad_rho_b, ec_supbeot_test)
 implicit none
 BEGIN_DOC
! exchange energy with the lda functional
 END_DOC
 integer :: i,j,m
 double precision, intent(in) :: rho_a(n_points_final_grid), rho_b(n_points_final_grid), grad_rho_a(3,n_points_final_grid), grad_rho_b(3,n_points_final_grid)
 double precision, intent(out) :: ec_supbeot_test
 double precision :: mu,weight
 double precision :: ecmd_ot
 double precision :: grad_rho_a_2,grad_rho_b_2,grad_rho_a_b

 ec_supbeot_test = 0.d0
 mu = 999 !mu_of_r

  do i = 1, n_points_final_grid
   weight = final_weight_at_r_vector(i)
   grad_rho_a_2 = 0.d0
   grad_rho_b_2 = 0.d0
   grad_rho_a_b = 0.d0
   do m = 1, 3
    grad_rho_a_2 += grad_rho_a(m,i) * grad_rho_a(m,i)
    grad_rho_b_2 += grad_rho_b(m,i) * grad_rho_b(m,i)
    grad_rho_a_b += grad_rho_a(m,i) * grad_rho_b(m,i)
   enddo

   ecmd_ot = 999! call ecmd_ot ...

   ec_supbeot_test += ecmd_ot * weight
 
  enddo
end

subroutine int_potential_c_supbeot_test (delta_rho_11, int_vc_supbeot_test)
 implicit none
 BEGIN_DOC
! delta_Ex[n] = int_dr delta_n(r) dEx/dn(r) 
! delta_n(r) <- subroutine delta_density_for_energy_test
! d_e_d_n    <- subroutine de_dn_at_r
 END_DOC
 integer :: i
 double precision, intent(in) :: delta_rho_11
 double precision, intent(out):: int_vc_supbeot_test
 double precision :: vc_i_j, vc_i_j_n2
 double precision :: delta_gamma_i_j(mo_num, mo_num)
 integer :: j,k,l,m
 
 int_vc_supbeot_test = 0.d0
 call delta_gamma_i_j_for_energy_test (delta_rho_11,delta_gamma_i_j)
 do k=1, mo_num
  do l=1, mo_num
   call potential_c_supbeot(j,k,l,m,vc_i_j, vc_i_j_n2)
   int_vc_supbeot_test += delta_gamma_i_j(k,l)*vc_i_j
  enddo
 enddo
end

subroutine int_potential_c_supbeot_n2_test (delta_rho_11, int_vc_supbeot_n2_test)
 implicit none
 BEGIN_DOC
! delta_Ex[n] = int_dr delta_n(r) dEx/dn(r) 
! delta_n(r) <- subroutine delta_density_for_energy_test
! d_e_d_n    <- subroutine de_dn_at_r
 END_DOC
 integer :: i
 double precision, intent(in) :: delta_rho_11
 double precision, intent(out):: int_vc_supbeot_n2_test
 double precision :: vc_i_j, vc_i_j_n2
 double precision :: delta_gamma_i_j_k_l(mo_num, mo_num,mo_num,mo_num)
 integer :: j,k,l,m
 
 int_vc_supbeot_n2_test = 0.d0
 call delta_gamma_i_j_k_l_for_energy_test (delta_rho_11,delta_gamma_i_j_k_l)
 do j=1, mo_num
  do k=1, mo_num
   do l=1, mo_num
    do m=1, mo_num
     call potential_c_supbeot(j,k,l,m,vc_i_j, vc_i_j_n2)
     int_vc_supbeot_n2_test += delta_gamma_i_j_k_l(j,k,l,m)*vc_i_j_n2
    enddo
   enddo
  enddo
 enddo
end



subroutine potential_c_supbeot(j,k,l,m,vc_i_j, vc_i_j_n2)
 implicit none
 BEGIN_DOC
! dEx/dn(r) = sum(k,l) {(vx_lda*mo_k*mo_l)}
! vx_lda = vx_lda_alpha + vx_lda_beta
 END_DOC
 double precision, intent(out) :: vc_i_j, vc_i_j_n2
 integer, intent(in) :: j,k,l,m
 integer :: i
 double precision  :: r(3), weight
 double precision :: mo_j(mo_num), mo_k(mo_num), mo_l(mo_num), mo_m(mo_num)
 double precision :: potential_c_alpha_mo_supbeot(mo_num,mo_num),potential_c_beta_mo_supbeot(mo_num,mo_num),potential_c_alpha_mo_supbeot_n2(mo_num,mo_num,mo_num,mo_num), potential_c_beta_mo_supbeot_n2(mo_num,mo_num,mo_num,mo_num)
 vc_i_j = 0.d0
 vc_i_j_n2 = 0.d0
 
 ! Où sont les providers suivants ?
 potential_c_alpha_mo_supbeot = 0.d0
 potential_c_beta_mo_supbeot = 0.d0
 potential_c_alpha_mo_supbeot_n2 = 0.d0
 potential_c_beta_mo_supbeot_n2 = 0.d0

 do i = 1, n_points_final_grid
   r(1) = final_grid_points(1,i)   
   r(2) = final_grid_points(2,i) 
   r(3) = final_grid_points(3,i)
   weight = final_weight_at_r_vector(i)

   call give_all_mos_at_r(r, mo_j)
   call give_all_mos_at_r(r, mo_k)
   call give_all_mos_at_r(r, mo_l)
   call give_all_mos_at_r(r, mo_m)

   vc_i_j += weight*(potential_c_alpha_mo_supbeot(k,l) + potential_c_beta_mo_supbeot(k,l)) * mo_k(k) * mo_l(l)
   vc_i_j_n2 += weight*(potential_c_alpha_mo_supbeot_n2(j,k,l,m) + potential_c_beta_mo_supbeot_n2(j,k,l,m)) * mo_j(j) * mo_k(k)* mo_l(l) * mo_m(m)
 enddo
end
 
!BEGIN_PROVIDER [double precision, potential_c_alpha_mo_supbeot, (mo_num,mo_num)]
!implicit none
!BEGIN_DOC
!!
!END_DOC
!! call ao_to_mo(potential_c_alpha_ao_supbeot,ao_num,potential_c_alpha_mo_supbeot,mo_num)
!potential_c_alpha_mo_supbeot(
!END_PROVIDER 

!BEGIN_PROVIDER [double precision, potential_c_beta_mo_supbeot, (mo_num,mo_num)]
!implicit none
!BEGIN_DOC
!!
!END_DOC
!call ao_to_mo(potential_c_beta_ao_supbeot,ao_num,potential_c_beta_mo_supbeot,mo_num)
!END_PROVIDER

!BEGIN_PROVIDER [double precision, potential_c_alpha_mo_supbeot_n2, (mo_num,mo_num)]
!implicit none
!BEGIN_DOC
!!
!END_DOC
!call ao_to_mo(potential_c_alpha_ao_supbeot,ao_num,potential_c_alpha_mo_supbeot,mo_num)
!END_PROVIDER 

!BEGIN_PROVIDER [double precision, potential_c_alpha_mo_supbeot_n2, (mo_num,mo_num)]
!implicit none
!BEGIN_DOC
!!
!END_DOC
!call ao_to_mo(potential_c_alpha_ao_supbeot,ao_num,potential_c_alpha_mo_supbeot,mo_num)
!END_PROVIDER 

