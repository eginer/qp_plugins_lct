program tc_scf

  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC

  implicit none

  print *, 'starting ...'

  my_grid_becke  = .True.
  my_n_pt_r_grid = 30
  my_n_pt_a_grid = 50
  touch my_grid_becke my_n_pt_r_grid my_n_pt_a_grid

  call routine_scf()

end

subroutine routine_scf()

  implicit none
  integer          :: i, j, it
  double precision :: e_save, e_delta

  mo_l_coef = fock_tc_leigvec_ao
  mo_r_coef = fock_tc_reigvec_ao
  call ezfio_set_bi_ortho_mos_mo_l_coef(mo_l_coef)
  call ezfio_set_bi_ortho_mos_mo_r_coef(mo_r_coef)
  
end subroutine routine_scf

! ---

