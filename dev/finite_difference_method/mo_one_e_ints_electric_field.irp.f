BEGIN_PROVIDER [ double precision, mo_one_e_integrals_electric_field,(mo_num,mo_num)]
  implicit none
  integer                        :: i,j,n,l
  BEGIN_DOC
  ! array of the one-electron Hamiltonian on the |MO| basis :
  ! sum of the kinetic and nuclear electronic potentials (and pseudo potential if needed)
  END_DOC
  print*,'Providing the one-electron integrals'

  IF (read_mo_one_e_integrals) THEN
        call ezfio_get_mo_one_e_ints_mo_one_e_integrals(mo_one_e_integrals_electric_field)
  ELSE
      mo_one_e_integrals_electric_field  = mo_integrals_n_e + mo_kinetic_integrals + field_strenght*mo_dipole_z

  ENDIF

  IF (write_mo_one_e_integrals) THEN
        call ezfio_set_mo_one_e_ints_mo_one_e_integrals(mo_one_e_integrals_electric_field)
       print *,  'MO one-e integrals written to disk'
  ENDIF

END_PROVIDER
