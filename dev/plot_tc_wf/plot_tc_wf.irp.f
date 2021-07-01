program plot_tc_wf
  implicit none
  read_wf = .True.
  touch read_wf
  call routine

end

subroutine routine
implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  print *, 'Hello world'
  print*,'norm_n2_jastrow_cst_mu     = ',norm_n2_jastrow_cst_mu
  print*,'norm_n2_inv_jastrow_cst_mu = ',norm_n2_inv_jastrow_cst_mu
  call print_psi_right_psi_trans
  call print_psi_left_psi_trans
end
