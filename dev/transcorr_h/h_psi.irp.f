program transcorr_h
 implicit none
 read_wf = .True.
 touch read_wf
 my_grid_becke = .True. 
 my_n_pt_r_grid = 30
 my_n_pt_a_grid = 170
 touch  my_grid_becke my_n_pt_r_grid my_n_pt_a_grid 
 extra_grid_type_sgn = 1 
 touch extra_grid_type_sgn 
 my_extra_grid_becke = .False.
 touch my_extra_grid_becke 
 print*,'Warning : the Becke grid parameters are automatically set to '
 print*,'my_n_pt_a_grid = ',my_n_pt_a_grid
 print*,'my_n_pt_r_grid = ',my_n_pt_r_grid
 print*,'If you want to modify them, you have to modify the following file '
 print*,'qp2/plugins/qp_plugins_lct/dev/transcorr_h/transcorr_general.irp.f'
 print*,'and recompile doing ninja'
 if(linear_tc)then
  three_body_h_tc = .False. 
  touch three_body_h_tc
  grad_squared = .False. 
  touch grad_squared 
 endif
 if(read_tc_ints)then
  call read_fcidump_1_tc
 endif

 call h_psi_tc
end

subroutine h_psi_tc
 use bitmasks
 implicit none
 double precision, allocatable :: htilde_psi(:)
 double precision :: a
 a = 1.d0
 allocate(htilde_psi(N_det))
 call h_non_hermite(htilde_psi,psi_coef(1,1),htilde_matrix_elmt,a,1,N_det)   
 print*,'htilde_psi/c0',htilde_psi(1)/psi_coef(1,1)
 print*,'htilde_psi'
 integer :: i
 do i = 1, N_det
  print*,htilde_psi(i)
 enddo
end

