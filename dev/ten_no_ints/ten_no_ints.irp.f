program ten_no_ints
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  print *, 'Hello world'
! call test_overlap_x_gauss_r12_ao
! call test_overlap_x_gauss_ten_no
! call test_deriv_ints
! call test_dgemm 
! print*,'n_max_fit_ten_no_slat = ',n_max_fit_ten_no_slat
! call full_num_deriv_ao
 call full_num_deriv_mo
! call test_dgemm_square
! call test_dgemm_lapl
! call full_num_lapl_ao
! call full_num_lapl_mo
! call full_num_square_ao
! call full_num_square_mo
 call test_mo_final
 call full_num_square_lapl_mo
end
