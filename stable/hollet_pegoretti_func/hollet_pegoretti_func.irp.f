program hollet_pegoretti_func
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  integer :: i_point
  double precision :: weight
  integer :: i,j
  double precision :: accu_1,accu_2,on_top,rho,v,V_ab_holl_peg
  accu_1 = 0.d0
  accu_2 = 0.d0
  do i_point = 1, n_points_final_grid
   weight = final_weight_at_r_vector(i_point)
   on_top = core_inact_act_on_top_of_r(i_point,1)
   rho    = one_e_dm_alpha_at_r(i_point,1) + one_e_dm_beta_at_r(i_point,1)
   v      = V_ab_holl_peg(rho,on_top)
   if(v.gt.0.d0)then
    accu_1 += weight * v
   else
    accu_2 += weight * v
   endif
  enddo
  print*,'accu_1 = ',accu_1
  print*,'accu_2 = ',accu_2
end
