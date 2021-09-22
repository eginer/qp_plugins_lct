program print_delta
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

 call routine
end

subroutine routine
 implicit none
 integer :: i
 double precision, allocatable :: delta_u0(:), h_u0(:), delta_mat(:,:)
 allocate(delta_u0(N_det),h_u0(N_det))
 delta_mat = htilde_matrix_elmt - h_matrix_all_dets ! Delta = Htilde - H
 delta_u0 = 0.d0
 double precision :: a
 a = 1.d0
 call h_non_hermite(delta_u0,psi_coef,delta_mat,a,1,N_det)  ! delta_u0 = Delta |u0> 
 print*,''
 print*,'printing Delta'
 print*,''
 do i = 1, N_det
  print*,delta_u0(i)
 enddo

 h_u0 = 0.d0
 call h_non_hermite(h_u0,psi_coef,h_matrix_all_dets,a,1,N_det)  ! delta_u0 = Delta |u0> 
 do i = 1, N_det
  print*,'delta_u0(i) + h_u0(i) / c_0(i) = ',(delta_u0(i) + h_u0(i))/psi_coef(i,1)
 enddo
end

