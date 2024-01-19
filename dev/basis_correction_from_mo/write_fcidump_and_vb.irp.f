program fcidump_and_vbarb
 implicit none
 read_wf = .true.
 touch read_wf
 ! total one-e integrals 
 io_mo_one_e_integrals = "None"
 touch io_mo_one_e_integrals  
 ! Vne integrals on the MO basis 
 io_mo_integrals_n_e = "None"
 touch io_mo_integrals_n_e
 ! kinetic integrals on the MO basis 
 io_mo_integrals_kinetic = "None"
 touch io_mo_integrals_kinetic 
 ! Vne integrals on the AO basis 
 io_ao_integrals_n_e = "None"
 touch io_ao_integrals_n_e 
 ! kinetic integrals on the AO basis 
 io_ao_integrals_kinetic = "None"
 touch io_ao_integrals_kinetic 

 ! regular 1/r12 integrals  on the MO basis
 io_mo_two_e_integrals = "None"
 touch io_mo_two_e_integrals
 ! regular 1/r12 integrals  on the AO basis
 io_ao_two_e_integrals = "None"
 touch io_ao_two_e_integrals
 ! integral of the effective potential 
! io_mo_int_mu_of_r = "None" 
! touch io_mo_int_mu_of_r

 no_core_density = .True.
 touch no_core_density

 call write_fcidump_vbarb
 call write_only_vbarb
end

subroutine write_only_vbarb
 implicit none
 integer :: i,j
 !          set directory             variable 
 double precision, allocatable :: pot_a(:,:)
 allocate(pot_a(mo_num, mo_num))
 do i = 1, mo_num
  do j = 1, mo_num
   pot_a(j,i) = pot_basis_alpha_mo(j,i,1)
  enddo
 enddo
 call ezfio_set_basis_correction_from_mo_v_one_e_a_mo( pot_a)

 do i = 1, mo_num
  do j = 1, mo_num
   pot_a(j,i) = pot_basis_beta_mo(j,i,1)
  enddo
 enddo
 call ezfio_set_basis_correction_from_mo_v_one_e_b_mo( pot_a)

 
end

