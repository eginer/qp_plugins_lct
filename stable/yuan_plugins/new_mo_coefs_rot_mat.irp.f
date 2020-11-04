program yuan_plugins
  implicit none
  BEGIN_DOC
! Assumptions : 
!
! There is a file called "rotation_matrix" 
  END_DOC
  print *, 'Hello world'
  call save_new_mo_coefs
  call read_one_rdm_sp_tr_uniq_and_write_to_ezfio(n_act_orb,dim_one_rdm)
!  call read_two_rdm_yuan_and_write_to_ezfio(dim_two_rdm,n_act_orb)
  call read_two_rdm_yuan_uniq_and_write_to_ezfio(dim_two_rdm,n_act_orb)

!  call read_two_rdm_and_write_to_ezfio(dim_two_rdm,n_act_orb)
!  call read_one_rdm_sp_tr_and_write_to_ezfio(n_act_orb)
end
