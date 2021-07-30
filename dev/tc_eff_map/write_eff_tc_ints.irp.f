program write_eff_tc_ints
 implicit none
 io_mo_two_e_integrals = "None"
 touch io_mo_two_e_integrals 
 read_wf = .True.
 touch read_wf
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

 call save_eff_tc_pot_in_coulomb
end

subroutine save_eff_tc_pot_in_coulomb
 implicit none
 integer :: i,j,k,l
 PROVIDE mo_two_e_integrals_tc_int_in_map
 call ezfio_set_work_empty(.False.)
 call map_save_to_disk(trim(ezfio_filename)//'/work/mo_ints',mo_integrals_tc_int_map)
 call ezfio_set_mo_two_e_ints_io_mo_two_e_integrals('Read')
end

