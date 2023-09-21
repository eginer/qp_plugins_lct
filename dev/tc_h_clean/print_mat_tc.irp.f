program print_mat_tc
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC

  implicit none

  my_grid_becke = .True.
  my_n_pt_r_grid = 30
  my_n_pt_a_grid = 50
  touch  my_grid_becke my_n_pt_r_grid my_n_pt_a_grid

  read_wf = .True.
  touch read_wf 

  call routine
  call routine_print_brutal
end

subroutine routine_print_brutal
 integer :: i,j
 do i = 1,N_det
  do j = 1,N_det
   print*,j,i,htilde_matrix_elmt(j,i)
  enddo
 enddo
 print*,'eigval_right_tc'
 do i = 1,N_det
  print*,i,eigval_right_tc(i) 
 enddo
end

subroutine routine
 implicit none
 integer :: i
 print*,''
 print*,'Total H_TC matrix '
 do i= 1, N_det
  write(*,'(100(F18.12,X))')htilde_matrix_elmt(i,:)
 enddo
 print*,''
 print*,'3-e   H_TC matrix '
 do i= 1, N_det
  write(*,'(100(F18.12,X))')htilde_matrix_elmt_hthree(i,:)
 enddo
end

subroutine routine_bis
 implicit none
 use bitmasks
 integer :: i,j
 double precision :: hmono, heff, hderiv, hthree, htot,ref
 ref = -0.41604835979406465D-01
 i = 2
 j = 4
 call htilde_mu_mat(psi_det(1,1,i), psi_det(1,1,j), N_int, hmono, heff, hderiv, hthree, htot)
 print*,'hmono+heff+hderiv = ', hmono+heff+hderiv
 print*,'hthree            = ',hthree
 print*,'***'
 print*,'htot              = ',htot
 print*,'ref               = ',ref
end
