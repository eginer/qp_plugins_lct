program tc_scf

  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC

  implicit none

  my_grid_becke  = .True.
  my_n_pt_r_grid = 30
  my_n_pt_a_grid = 50
  touch  my_grid_becke my_n_pt_r_grid my_n_pt_a_grid
 call routine_scf
end

subroutine routine_scf
 implicit none
 integer :: i
 i = 0
  print*,'iteration = ',i
  print*,'grad_good_hermit_tc_fock_mat = ',grad_good_hermit_tc_fock_mat
  print*,'***'
   print*,'TC HF total energy = ',TC_right_HF_energy
   print*,'TC HF 1 e   energy = ',TC_right_HF_one_electron_energy
   print*,'TC HF 2 e hermit   = ',TC_right_HF_two_e_hermit_energy
   print*,'TC HF 2 non hermit = ',TC_right_HF_two_e_n_hermit_energy
   print*,'TC HF 3 body       = ',diag_three_elem_hf
  print*,'***'
 do while(grad_good_hermit_tc_fock_mat.gt.thresh_scf)
  i += 1
! do i = 1, 10
  print*,'iteration = ',i
  print*,'grad_good_hermit_tc_fock_mat = ',grad_good_hermit_tc_fock_mat
  print*,'***'
   print*,'TC HF total energy = ',TC_right_HF_energy
   print*,'TC HF 1 e   energy = ',TC_right_HF_one_electron_energy
   print*,'TC HF 2 e hermit   = ',TC_right_HF_two_e_hermit_energy
   print*,'TC HF 2 non hermit = ',TC_right_HF_two_e_n_hermit_energy
   print*,'TC HF 3 body       = ',diag_three_elem_hf
  print*,'***'
  call save_good_hermit_tc_eigvectors
  touch mo_coef 
  call save_mos
 enddo
end

