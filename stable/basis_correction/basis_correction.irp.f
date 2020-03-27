program basis_correction
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  read_wf = .True.
  touch read_wf
  call print_e_b
end

subroutine print_e_b
  print *, 'Hello world'
  print*,'ecmd_lda_mu_of_r = ',ecmd_lda_mu_of_r
  print*,'psi_energy + E^B = ',psi_energy + ecmd_lda_mu_of_r
  print*,'mu_average_prov  = ',mu_average_prov
end
