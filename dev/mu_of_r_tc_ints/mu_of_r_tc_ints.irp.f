program mu_of_r_tc_ints
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  constant_mu = .False.
  touch constant_mu
!  call test_gauss_ij_rk
 call test_erf_mu_squared_ij_rk
end
