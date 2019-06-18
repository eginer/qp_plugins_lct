program hollet_pegoretti_func
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  integer :: i_point
  double precision :: weight
  integer :: i,j
  double precision :: accu_ab,accu_aa,on_top,rho,v,V_ab_holl_peg
  double precision :: lapl,on_top_aa_lapl,V_aa_holl_peg
  double precision :: accu_bb,on_top_bb_lapl
  accu_ab = 0.d0
  accu_aa = 0.d0
  accu_bb = 0.d0
  do i_point = 1, n_points_final_grid
   weight = final_weight_at_r_vector(i_point)
   on_top = core_inact_act_on_top_of_r(i_point,1)
   rho    = one_e_dm_alpha_at_r(i_point,1) + one_e_dm_beta_at_r(i_point,1)
   v      = V_ab_holl_peg(rho,on_top)
   accu_ab += weight * v
   
   lapl   = on_top_aa_lapl(i_point) 
   v      = V_aa_holl_peg(lapl)
   accu_aa += weight * v

   lapl   = on_top_bb_lapl(i_point) 
   v      = V_aa_holl_peg(lapl)
   accu_bb += weight * v
  enddo
  print*,'accu_ab = ',accu_ab
  print*,'accu_aa = ',accu_aa
  print*,'accu_bb = ',accu_bb
end

