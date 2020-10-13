program my_func
 implicit none
 read_wf = .True.
 touch read_wf 
 call test_x_lda
end


BEGIN_PROVIDER [double precision, delta_gamma, (mo_num, mo_num)]
 implicit none
 double precision :: eps
 integer :: i,j
 delta_gamma = 0.d0
!eps = delta_dm_input/dble(mo_num*mo_num)
!do i = 1, mo_num
! do j = 1, mo_num
!  delta_gamma(j,i) = eps
! enddo
!enddo
 eps = 1.d-5
 delta_gamma(1,1) = eps
END_PROVIDER 

BEGIN_PROVIDER [double precision, delta_gamm_pot_x_lda]
 implicit none
 integer :: i,j
 delta_gamm_pot_x_lda = 0.d0
 do i = 1, mo_num
  do j = 1, mo_num
   delta_gamm_pot_x_lda += delta_gamma(j,i) * potential_x_mo_lda(j,i)
  enddo
 enddo
END_PROVIDER 

BEGIN_PROVIDER [double precision, delta_gamm_pot_x_pbe]
 implicit none
 integer :: i,j
 delta_gamm_pot_x_pbe = 0.d0
 do i = 1, mo_num
  do j = 1, mo_num
   delta_gamm_pot_x_pbe += delta_gamma(j,i) * (potential_x_alpha_mo_pbe(j,i) + potential_x_beta_mo_pbe(j,i))
  enddo
 enddo
END_PROVIDER 

double precision function give_operator_at_r(r,op_mat_mo)
 implicit none
 double precision, intent(in) :: r(3),op_mat_mo(mo_num,mo_num)
 integer :: i,j
 double precision :: mos_array(mo_num)
 call give_all_mos_at_r(r, mos_array)
 give_operator_at_r = 0.d0
 do i = 1, mo_num
  do j = 1, mo_num
   give_operator_at_r += mos_array(j) * mos_array(i) * op_mat_mo(j,i)
  enddo
 enddo
end

BEGIN_PROVIDER [double precision, potential_x_mo_lda, (mo_num, mo_num) ]
 implicit none
 potential_x_mo_lda = potential_x_alpha_mo_lda  + potential_x_beta_mo_lda 
END_PROVIDER 

subroutine test_x_lda
 implicit none
 integer :: i
 double precision :: n_delta_dm_a(n_points_final_grid), n_delta_dm_b(n_points_final_grid)
 double precision :: delta_dm_a(n_points_final_grid), delta_dm_b(n_points_final_grid)
 double precision :: delta_dm(n_points_final_grid)
 double precision :: dm_a(n_points_final_grid), dm_b(n_points_final_grid)
 double precision :: d_exlda_dn(n_points_final_grid)
 double precision :: d_exlda_dn_a(n_points_final_grid)
 double precision :: d_exlda_dn_b(n_points_final_grid)
 double precision :: give_operator_at_r,r(3),delta_n, d_ex,delta
 double precision :: ex_lda_test,ex_lda_test_plus_delta, delta_num
 do i = 1, n_points_final_grid
  r(1) = final_grid_points(1,i) 
  r(2) = final_grid_points(2,i) 
  r(3) = final_grid_points(3,i)
  call dm_dft_alpha_beta_at_r(r,dm_a(i),dm_b(i))
  delta_n = give_operator_at_r(r,delta_gamma)
  n_delta_dm_a(i) =  delta_n + dm_a(i)
  n_delta_dm_b(i) =  delta_n + dm_b(i)
  delta_dm(i) = 1.d0 * delta_n  
  delta_dm_a(i) =  delta_n 
  delta_dm_b(i) =  delta_n 
  d_ex = give_operator_at_r(r,potential_x_mo_lda)
  d_exlda_dn(i) = d_ex
  d_ex = give_operator_at_r(r,potential_x_alpha_mo_lda)
  d_exlda_dn_a(i) = d_ex
  d_ex = give_operator_at_r(r,potential_x_beta_mo_lda)
  d_exlda_dn_b(i) = d_ex
 enddo
 call energy_x_lda_test(dm_a, dm_b, ex_lda_test)
 call energy_x_lda_test(n_delta_dm_a, n_delta_dm_b, ex_lda_test_plus_delta)
 delta = ex_lda_test_plus_delta - ex_lda_test
 print*,'ex_lda_test            = ',ex_lda_test
 print*,'ex_lda_test_plus_delta = ',ex_lda_test_plus_delta
 print*,'delta                  = ',delta
 print*,'delta_gamm_pot_x_lda   = ',delta_gamm_pot_x_lda
 print*,'delta/delta_gamma_x_lda= ',delta/delta_gamm_pot_x_lda
 print*,'delta_gamm_pot_x_pbe   = ',delta_gamm_pot_x_pbe
 stop
 call integrate_prods(delta_dm, d_exlda_dn,delta_num)
 print*,'delta_num              = ',delta_num
 call integrate_prods(delta_dm_a, d_exlda_dn_a,delta_num)
 print*,'delta_num              = ',delta_num
 call integrate_prods(delta_dm_b, d_exlda_dn_b,delta_num)
 print*,'delta_num              = ',delta_num
 print*,'2 * delta_num          = ',2.d0 * delta_num

end


subroutine integrate_prods(array_a, array_b,ints)
 implicit none
 double precision, intent(in)  :: array_a(n_points_final_grid), array_b(n_points_final_grid)
 double precision, intent(out) :: ints
 integer :: i
 double precision :: weight
 ints = 0.d0
 do i = 1, n_points_final_grid
  weight = final_weight_at_r_vector(i)
  ints += array_a(i) * array_b(i) * weight
 enddo
end
