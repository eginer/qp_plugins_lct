program test_det
 implicit none
 read_wf = .True.
 touch read_wf 
! call routine
! call test_conf
! call test_prov_s_p
 call test_prov_s_m
end

subroutine test_prov_s_p
 implicit none
 integer :: i,j,k,istate
 print*,'Initial wave function '
 do istate = 1,N_states
  print*,'********************'
  print*,'state ',istate
  print*,' S values ',s_values(istate)
  do i = 1, N_configuration
   do j =  psi_configuration_to_psi_det(1,i), psi_configuration_to_psi_det(2,i)
    k = psi_configuration_to_psi_det_data(j) ! index of the determinant "j" in the configuration i 
    print*,'coef = ',psi_coef(k,istate)
    call debug_det(psi_det(1,1,k),N_int)
   enddo 
  enddo
  print*,'Applying S^+ on the state'
  do i = 1, N_configuration
   do j = 1, n_det_per_conf_s_p(i) 
    call debug_det(psi_det_s_p(1,1,j,i),N_int)
    print*,'coef = ',psi_coef_s_p(j,i,istate)
   enddo
  enddo
  print*,'psi_norm_s_p',psi_norm_s_p(istate)
  print*,'s2_values',s2_values
  print*,'s_values ',s_values
 enddo
end

subroutine test_prov_s_m
 implicit none
 integer :: i,j,k,istate
 print*,'Initial wave function '
 do istate = 1,N_states
  print*,'********************'
  print*,'state ',istate
  print*,' S values ',s_values(istate)
  do i = 1, N_configuration
   do j =  psi_configuration_to_psi_det(1,i), psi_configuration_to_psi_det(2,i)
    k = psi_configuration_to_psi_det_data(j) ! index of the determinant "j" in the configuration i 
    print*,'coef = ',psi_coef(k,istate)
    call debug_det(psi_det(1,1,k),N_int)
   enddo 
  enddo
  print*,'Applying S^- on the state'
  do i = 1, N_configuration
   do j = 1, n_det_per_conf_s_m(i) 
    call debug_det(psi_det_s_m(1,1,j,i),N_int)
    print*,'coef = ',psi_coef_s_m(j,i,istate)
   enddo
  enddo
  print*,'psi_norm_s_m',psi_norm_s_m(istate)
 enddo
end

subroutine test_conf
 implicit none
  use bitmasks ! you need to include the bitmasks_module.f90 features
 integer                        :: i,j,k,l,m,ndet_conf,ndet_conf_sp,n_alpha_p,ndet_out,degree
 integer(bit_kind), allocatable :: det_tmp(:,:),dets_sp(:,:,:)
 integer(bit_kind), allocatable :: det_out(:,:,:)
 double precision, allocatable  :: phase(:),coef_sp(:)
 allocate(det_tmp(N_int,2))
 allocate(phase(elec_beta_num),det_out(N_int,2,elec_beta_num))
 do i = 1, N_configuration
  det_tmp(:,:) = psi_configuration(:,:,i)
  call debug_det(det_tmp,N_int)
  ndet_conf = psi_configuration_n_det(i)
  n_alpha_p = n_elec_alpha_for_psi_configuration(i)
  n_alpha_p = elec_alpha_num + 1
  call configuration_to_dets_size(psi_configuration(1,1,i),ndet_conf_sp,n_alpha_p,N_int)
  allocate(dets_sp(N_int,2,ndet_conf_sp),coef_sp(ndet_conf_sp))
  coef_sp = 0.d0
  call configuration_to_dets(psi_configuration(1,1,i),dets_sp,ndet_conf_sp,n_alpha_p,N_int)
  do j =  psi_configuration_to_psi_det(1,i), psi_configuration_to_psi_det(2,i)
   k = psi_configuration_to_psi_det_data(j)
   call s_plus_det(psi_det(1,1,k),det_out,phase,ndet_out)
   do m = 1, ndet_conf_sp
    do l = 1, ndet_out
     call get_excitation_degree(dets_sp(1,1,m),det_out(1,1,l),degree,N_int)
     if(degree == 0)then
      coef_sp(m) += psi_coef(k,1) * phase(m)
     endif
    enddo
   enddo
  enddo
 enddo
 print*,'coef_sp = ',coef_sp
end

subroutine routine
 implicit none
  use bitmasks ! you need to include the bitmasks_module.f90 features
 integer                        :: exc(0:2,2,2)
 integer                        :: h2,p2,s2,degree,i,j,ndet_out
 integer(bit_kind), allocatable :: det_out(:,:,:)
 double precision, allocatable  :: phase(:)
 double precision :: accu
 allocate(phase(elec_beta_num),det_out(N_int,2,elec_beta_num))
 accu = 0.d0
 do i = 1, N_det
  call s_plus_det(psi_det(1,1,i),det_out,phase,ndet_out)
  print*,''
  print*,'Applying S+ on '
  call debug_det(psi_det(1,1,i),N_int)
  do j = 1, ndet_out
   print*,'phase = ',phase
   call debug_det(det_out(1,1,j),N_int)
   accu += psi_coef(i,1) * phase(j)
  enddo
 enddo
 print*,'accu = ',accu
 print*,'sqrt(3)',dsqrt(3.d0)
end
