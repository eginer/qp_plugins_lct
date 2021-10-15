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

 call provide_all
 call print_energy_tc
 call print_e_comp_transcorr
 call print_eigv
 call print_pert
 call write_left_right
 call routine_save
end

subroutine provide_all
 use bitmasks
 integer(bit_kind) :: key_i(N_int,2), key_j(N_int,2)
 integer :: i,j,degree
 double precision :: hij,s2,hmono,heff,hderiv,htot,hthree
 double precision :: accu
 accu = 0.d0
 key_i(:,:) = psi_det(:,:,1)
 call htilde_mat(key_i,key_i,hmono,heff,hderiv,hthree,htot)
 provide eigval_trans
end


subroutine routine_save
implicit none
 double precision, allocatable :: coef_tmp(:,:)
 N_states = 1

 allocate(coef_tmp(N_det, N_states))
 integer :: i
 do i = 1, N_det
  coef_tmp(i,1) = reigvec_trans(i,1)
!  coef_tmp(i,2) = leigvec_trans(i,1)
 enddo
 call save_wavefunction_general(N_det,N_states,psi_det,size(coef_tmp,1),coef_tmp(1,1))
end
