program yuan_plugins
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  print *, 'Hello world'
  call save_new_mo_coefs
  call read_two_rdm_and_write_to_ezfio(dim_two_rdm,n_act_orb)
end
