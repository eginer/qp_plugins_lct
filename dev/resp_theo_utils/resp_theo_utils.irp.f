program resp_theo_utils
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  print *, 'Hello world'
  integer :: i
  do i = 1, n_singles_resp
   print*,'eigval_a_resp_h_mat(i) = ',eigval_a_resp_h_mat(i)
  enddo
  print*,''
  print*,''
  print*,''
  print*,''
  provide eigval_ab_resp_hf_mat
end
