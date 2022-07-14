BEGIN_PROVIDER [ double precision, tc_scf_dm_in_r, (n_points_final_grid) ]
 implicit none
 integer :: i,j
 tc_scf_dm_in_r = 0.d0
  do i = 1, n_points_final_grid
   do j = 1, elec_beta_num
    tc_scf_dm_in_r(i) += mos_r_in_r_array(j,i) * mos_l_in_r_array(j,i)
   enddo
  enddo
END_PROVIDER 

BEGIN_PROVIDER [ double precision, w_sum_in_r, (n_points_final_grid,3)]
 implicit none
 integer :: i,j,xi
 w_sum_in_r = 0.d0
 do j = 1, elec_beta_num
  do xi = 1, 3
   do i = 1, n_points_final_grid
    w_sum_in_r(i,xi) += x_W_ki_bi_ortho_erf_rk(i,xi,j,j)
   enddo
  enddo
 enddo
END_PROVIDER 

subroutine direct_term_bi_ortho(a,i,integral)
 implicit none
 double precision, intent(out) :: integral
 integer, intent(in) :: i,a
 BEGIN_DOC
! computes sum_(j,m =1, elec_beta_num) <a m j | i m j > with bi ortho mos
 END_DOC
 integer :: ipoint,xi
 double precision :: weight
 integral = 0.d0
 do xi = 1, 3
  do ipoint = 1, n_points_final_grid
   weight = final_weight_at_r_vector(ipoint)
   integral += (mos_l_in_r_array(a,ipoint) * mos_r_in_r_array(i,ipoint) * w_sum_in_r(ipoint,xi) * w_sum_in_r(ipoint,xi) & 
             + 2.d0 * tc_scf_dm_in_r(ipoint) * w_sum_in_r(ipoint,xi) * x_W_ki_bi_ortho_erf_rk(ipoint,xi,a,i) ) * weight
  enddo
 enddo
end
