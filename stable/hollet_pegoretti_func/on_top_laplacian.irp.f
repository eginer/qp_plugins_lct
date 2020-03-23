double precision function on_top_aa_lapl(i_point)
 implicit none
 BEGIN_DOC
 ! Laplacian of the 
 END_DOC
 integer, intent(in) :: i_point
 double precision :: r(3),nabla_2_at_r_mo(n_act_orb,n_act_orb)
 r(1)= final_grid_points(1,i_point) 
 r(2)= final_grid_points(2,i_point) 
 r(3)= final_grid_points(3,i_point)
 call give_nabla_2_at_r_act_mo(r,nabla_2_at_r_mo)
 
 double precision :: direct,exchange,lapl_jj,lapl_ij
 integer :: i,j 
 direct   = 0.d0
 exchange = 0.d0
 do i = 1, elec_alpha_num
  do j = 1, elec_alpha_num
   lapl_jj   = nabla_2_at_r_mo(j,j)
   lapl_ij   = nabla_2_at_r_mo(j,i)
   direct   += mos_in_r_array(i,i_point)**2 * lapl_jj
   exchange += mos_in_r_array(i,i_point) * mos_in_r_array(j,i_point) * lapl_ij
  enddo
 enddo
 on_top_aa_lapl = 0.5d0 * (direct - exchange)
end

double precision function on_top_bb_lapl(i_point)
 implicit none
 integer, intent(in) :: i_point
 double precision :: r(3),nabla_2_at_r_mo(n_act_orb,n_act_orb)
 r(1)= final_grid_points(1,i_point) 
 r(2)= final_grid_points(2,i_point) 
 r(3)= final_grid_points(3,i_point)
 call give_nabla_2_at_r_act_mo(r,nabla_2_at_r_mo)
 
 double precision :: direct,exchange,lapl_jj,lapl_ij
 integer :: i,j 
 direct   = 0.d0
 exchange = 0.d0
 do i = 1, elec_beta_num
  do j = 1, elec_beta_num
   lapl_jj   = nabla_2_at_r_mo(j,j)
   lapl_ij   = nabla_2_at_r_mo(j,i)
   direct   += mos_in_r_array(i,i_point)**2 * lapl_jj
   exchange += mos_in_r_array(i,i_point) * mos_in_r_array(j,i_point) * lapl_ij
  enddo
 enddo
 on_top_bb_lapl = 0.5d0 * (direct - exchange)
end
