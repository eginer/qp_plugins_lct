program diag_h_t_iter
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

 call routine_diagonalize_htilde
end

subroutine routine_diagonalize_htilde
 implicit none
 integer :: j,i,dressing_state,idress
 double precision :: res,ebefore,thr,hij
 double precision, allocatable :: u1(:),u0(:,:),H_jj(:),dressing_vec(:),H_jj_tmp(:),delta(:)
 double precision, allocatable :: energies(:)
 logical :: converged 
 external hcalc_template
 call provide_integrals_for_tc
 thr = threshold_davidson
 allocate(u1(N_det),u0(N_det,N_states_diag),H_jj(N_det),dressing_vec(N_det),H_jj_tmp(N_det),delta(N_det))
 allocate(energies(N_states_diag))
 do i = 1, N_det
  call  i_H_j(psi_det(1,1,i),psi_det(1,1,i),N_int,hij)
  H_jj(i) = hij
 enddo
 dressing_state = 1
 res = 1.d0
 j = 0
 u0 = 0.d0
 u0(1:N_det,1) = psi_coef(1:N_det,1)
 idress = 1
 energies = 0.d0
 do while(res .gt. thr) 
  j += 1
  print*,''
  print*,'**************'
  print*,'**************'
  print*,'**************'
  print*,'Iteration j ',j
  H_jj_tmp = H_jj
  ! delta(i) = <I|(H_t - H)|u0>
  call get_delta_no_store_no_provide(psi_det,u0(1,1),n_det,delta)
  ! dressing_vec(i) = dressing vector for the davidson dress
  call delta_to_dressing_vector(u0(1,1),delta,n_det,dressing_vec,idress)
  ! get the ground state energy with Davidson including a dressing 
  call davidson_general_ext_rout_dressed(u0,H_jj_tmp,energies,N_det,N_states,N_states_diag,dressing_state,dressing_vec,idress,converged,hcalc_template)
  print*,'**************'
  if(j.gt.1)then
   print*,'energies,De = ',energies(1),dabs(ebefore-energies(1))
  else
   print*,'energies    = ',energies(1)
  print*,'**************'
  print*,'**************'
  print*,'**************'
  endif
  u1 = 0.d0
  call htilde_psi_no_store(psi_det,u0(1,1),n_det,u1)
  u1(:) -= energies(1) * u0(:,1)
  res = 0.d0
  do i = 1, N_det
   res += u1(i)*u1(i)
  enddo
  res = dsqrt(res)
  print*,'Norm of the residual vector ', res 
  ebefore = energies(1)
 enddo
 call routine_save(u0)
end

subroutine routine_save(u0)
implicit none
 double precision, intent(in)  :: u0(N_det)
 double precision, allocatable :: coef_tmp(:,:)
 N_states = 1

 allocate(coef_tmp(N_det, N_states))
 integer :: i
 do i = 1, N_det
  coef_tmp(i,1) = u0(i)
 enddo
 call save_wavefunction_general(N_det,N_states,psi_det,size(coef_tmp,1),coef_tmp(1,1))
end

subroutine provide_integrals_for_tc
 implicit none
 PROVIDE scalar_mu_r_pot_physicist_mo deriv_mu_r_pot_physicist_mo
 PROVIDE three_body_3_index three_body_3_index_exch_12 three_body_3_index_exch_13 three_body_3_index_exch_23
 PROVIDE three_body_5_index three_body_5_index_exch_13 three_body_5_index_exch_32
 PROVIDE three_body_4_index three_body_4_index_exch_12 three_body_4_index_exch_12_part
end
