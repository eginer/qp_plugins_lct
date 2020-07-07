program mo_two_e_eff_pot
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  print *, 'Hello world'
! call test_gauss_ints_mos
! call test_gauss_ints_mos_exchange
 call test_coulomb_exchange
 call test_fit_pouet
end
