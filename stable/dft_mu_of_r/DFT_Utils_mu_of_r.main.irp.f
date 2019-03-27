program DFT_Utils_mu_of_r
  implicit none
  BEGIN_DOC
! TODO
  END_DOC
  read_wf = .True.
  touch read_wf 
  call routine 
end

subroutine routine
 implicit none

 io_mo_one_e_integrals = "None"
 touch io_mo_two_e_integrals
 io_ao_two_e_integrals = "None"
 touch io_ao_two_e_integrals



 print*,'**********************'
 print*,'**********************'
 print*,'LDA / HF coallescence with FULL DENSITY'
 mu_of_r_functional ="basis_set_LDA"
 touch mu_of_r_functional
 mu_of_r_potential = "hf_coallescence"
 touch mu_of_r_potential 
 no_core_density = "pouet"
 touch no_core_density
 call print_contribution_dft_mu_of_r
 print*,'**********************'
 print*,'**********************'
 print*,'**********************'
 print*,'LDA / HF VALENCE coallescence with no_core_density'
 mu_of_r_functional ="basis_set_LDA"
 touch mu_of_r_functional
 mu_of_r_potential = "hf_valence_coallescence"
 touch mu_of_r_potential 
 no_core_density = "no_core_dm"
 touch no_core_density
 call print_contribution_dft_mu_of_r
 print*,'**********************'
 print*,'**********************'


end


subroutine routine2
 implicit none
 integer :: i_point
 if(dabs(mu_of_r_hf_coal_vv_vector(i_point) - mu_of_r_hf_coal_vector(i_point)).gt.1.d-10)then
  print*,i_point
  print*,final_grid_points(:,i_point)
  print*,mu_of_r_hf_coal_vv_vector(i_point),mu_of_r_hf_coal_vector(i_point)
 endif
 
end
