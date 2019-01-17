program write_erf_and_regular_ints
 implicit none
 io_mo_two_e_integrals = "None"
 touch io_mo_two_e_integrals
 io_ao_two_e_integrals = "None"
 touch io_ao_two_e_integrals
 io_mo_two_e_integrals_erf = "None" 
 touch io_mo_two_e_integrals_erf
 io_ao_two_e_integrals_erf = "None" 
 touch io_ao_two_e_integrals_erf
 call ezfio_set_work_empty(.False.)
 call routine_provide 
 call routine_write
end

subroutine routine_provide
 implicit none
 provide mo_two_e_integrals_in_map
 provide mo_two_e_integrals_erf_in_map
end

subroutine routine_write
 implicit none
 call map_save_to_disk(trim(ezfio_filename)//'/work/mo_ints_regular',mo_integrals_map)
 call ezfio_set_mo_two_e_ints_io_mo_two_e_integrals('Read')

 call map_save_to_disk(trim(ezfio_filename)//'/work/mo_ints_erf',mo_integrals_erf_map)
 call ezfio_set_mo_two_e_erf_ints_io_mo_two_e_integrals_erf("Read")
end
