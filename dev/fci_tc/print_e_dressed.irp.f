
program fci
  implicit none
  BEGIN_DOC
  ! Selected Full Configuration Interaction with stochastic selection
  ! and PT2.
  !
  ! This program performs a |CIPSI|-like selected |CI| using a
  ! stochastic scheme for both the selection of the important Slater
  ! determinants and the computation of the |PT2| correction. This
  ! |CIPSI|-like algorithm will be performed for the lowest states of
  ! the variational space (see :option:`determinants n_states`). The
  ! |FCI| program will stop when reaching at least one the two following
  ! conditions:
  !
  ! * number of Slater determinants > :option:`determinants n_det_max`
  ! * abs(|PT2|) less than :option:`perturbation pt2_max`
  !
  ! The following other options can be of interest:
  !
  ! :option:`determinants read_wf`
  !   When set to |false|, the program starts with a ROHF-like Slater
  !   determinant as a guess wave function. When set to |true|, the
  !   program starts with the wave function(s) stored in the |EZFIO|
  !   directory as guess wave function(s).
  !
  ! :option:`determinants s2_eig`
  !   When set to |true|, the selection will systematically add all the
  !   necessary Slater determinants in order to have a pure spin wave
  !   function with an |S^2| value corresponding to
  !   :option:`determinants expected_s2`.
  !
  ! For excited states calculations, it is recommended to start with
  ! :ref:`cis` or :ref:`cisd` guess wave functions, eventually in
  ! a restricted set of |MOs|, and to set :option:`determinants s2_eig`
  ! to |true|.
  !
  END_DOC
  my_grid_becke = .True.
  my_n_pt_r_grid = 30
  my_n_pt_a_grid = 170
  touch  my_grid_becke my_n_pt_r_grid my_n_pt_a_grid
  extra_grid_type_sgn = 1
  touch extra_grid_type_sgn
  my_extra_grid_becke = .False.
  touch my_extra_grid_becke
  print*,'Warning : the Becke grid parameters are automatically set to '
  print*,'my_n_pt_a_grid = ',my_n_pt_a_grid
  print*,'my_n_pt_r_grid = ',my_n_pt_r_grid
  print*,'If you want to modify them, you have to modify the following file '
  print*,'qp2/plugins/qp_plugins_lct/dev/transcorr_h/transcorr_general.irp.f'
  print*,'and recompile doing ninja'
  if(linear_tc)then
   three_body_h_tc = .False.
   touch three_body_h_tc
   grad_squared = .False.
   touch grad_squared
  endif
 read_wf = .true.
 touch read_wf
 call routine
end


subroutine routine
 implicit none
! print*,'ci_energy_dressed = ',ci_energy_dressed
 call get_e_dressed_scf
! provide ci_energy_dressed_scf
! call test_tc

end
