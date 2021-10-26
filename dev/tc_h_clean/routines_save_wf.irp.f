subroutine routine_save_left_right
 implicit none
 double precision, allocatable :: coef_tmp(:,:)
 integer :: i
 if(comp_left_eigv)then
  N_states = 2
  allocate(coef_tmp(N_det, N_states))
  do i = 1, N_det
   coef_tmp(i,1) = reigvec_tc(i,1)
   coef_tmp(i,1) = leigvec_tc(i,1)
  enddo
 else
  N_states = 1
  allocate(coef_tmp(N_det, N_states))
  do i = 1, N_det
   coef_tmp(i,1) = reigvec_tc(i,1)
  enddo
 endif
 call save_wavefunction_general(N_det,N_states,psi_det,size(coef_tmp,1),coef_tmp(1,1))
end

subroutine routine_save_right
 implicit none
 double precision, allocatable :: coef_tmp(:,:)
 integer :: i
 N_states = 1
 allocate(coef_tmp(N_det, N_states))
 do i = 1, N_det
  coef_tmp(i,1) = reigvec_tc(i,1)
 enddo
 call save_wavefunction_general(N_det,N_states,psi_det,size(coef_tmp,1),coef_tmp(1,1))
end
