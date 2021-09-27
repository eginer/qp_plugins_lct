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
 double precision, allocatable ::  delta_u0_new(:)

 allocate(delta_u0_new(N_det))
 call get_delta_no_store(psi_det,psi_coef,n_det,delta_u0_new) ! delta_u0_new = Delta |u0> 

 print*,''
 print*,'printing Delta'
 print*,''
 call dset_order(delta_u0_new,psi_bilinear_matrix_order,N_det)

! do i = 1, N_det
!  print*,delta_u0(i)
! enddo

 call ezfio_set_dmc_dress_dmc_delta_h(delta_u0_new)

end

