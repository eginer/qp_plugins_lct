program mu_of_r_print
 implicit none
 read_wf = .True.
 touch read_wf
 no_core_density = "no_core_dm"
 touch no_core_density
 call routine 

end

subroutine routine 
  implicit none
  double precision :: z,ontop,ec_pbe_ontop,ec_pbe,dm_a,dm_b,mu
  integer :: i
  provide cas_full_mu_of_r_grid_mur_psi_coal_vector
  provide E_c_md_mur_grid_PBE
  provide E_c_md_mur_grid_n_and_on_top_PBE
  character*(128) :: output
  integer :: i_unit_output,getUnitAndOpen
 
  output=trim(ezfio_filename)//'.'//trim("mu_of_r")
  output=trim(output)
  print*,'output = ',trim(output)
  i_unit_output = getUnitAndOpen(output,'w')

  do i = 1, n_points_print_mur
   z = grid_points_mur(3,i) 
   call dm_dft_alpha_beta_at_r(grid_points_mur(1,i),dm_a,dm_b)
   ec_pbe_ontop = E_c_md_mur_grid_n_and_on_top_PBE(1,i)
   ec_pbe       = E_c_md_mur_grid_PBE(1,i)
   ontop        = core_inact_act_on_top_of_r_grid_mur(i,1)
   mu           = cas_full_mu_of_r_grid_mur_psi_coal_vector(i)
   write(i_unit_output,'(100(F16.10,X))')z,ec_pbe_ontop,ec_pbe,mu,ontop,(dm_a+dm_b)
  enddo
end
