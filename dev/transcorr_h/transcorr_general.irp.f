program transcorr_h
 implicit none
 read_wf = .True.
 touch read_wf
 my_grid_becke = .True. 
 my_n_pt_r_grid = 30
 my_n_pt_a_grid = 50
 touch  my_grid_becke my_n_pt_r_grid my_n_pt_a_grid 
 print*,'Warning : the Becke grid parameters are automatically set to '
 print*,'my_n_pt_a_grid = 50'
 print*,'my_n_pt_r_grid = 30'
 print*,'If you want to modify them, you have to modify the following file '
 print*,'qp2/plugins/qp_plugins_lct/dev/transcorr_h/transcorr_general.irp.f'
 print*,'and recompile doing ninja'
 if(linear_tc)then
  three_body_h_tc = .False. 
  touch three_body_h_tc
  grad_squared = .False. 
  touch grad_squared 
 endif

 call provide_all
 call print_energy_tc
 call print_e_comp_transcorr
 call print_eigv
 call print_pert
 call write_left_right
end

subroutine provide_all
 use bitmasks
 integer(bit_kind) :: key_i(N_int,2), key_j(N_int,2)
 integer :: i,j,degree
 double precision :: hij,s2,hmono,herf,heff,hderiv,htot,hthree
 double precision :: accu
 accu = 0.d0
 key_i(:,:) = psi_det(:,:,1)
 call htilde_mat(key_i,key_i,hmono,herf,heff,hderiv,hthree,htot)
 provide eigval_trans
end

