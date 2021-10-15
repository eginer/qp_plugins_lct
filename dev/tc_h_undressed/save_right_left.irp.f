program save_right_eigv
 implicit none
 read_wf = .True.
 touch read_wf 
 read_rl_eigv = .True.
 touch read_rl_eigv
 print*,'Saving right/left eigenvectors in EZFIO as the two first states'
 call routine_save
end

subroutine routine_save
implicit none
 double precision, allocatable :: coef_tmp(:,:)
 N_states = 2

 allocate(coef_tmp(N_det, 2))
 integer :: i
 do i = 1, N_det
  coef_tmp(i,1) = reigvec_trans(i,1)
  coef_tmp(i,2) = leigvec_trans(i,1)
 enddo
 call save_wavefunction_general(N_det,N_states,psi_det,size(coef_tmp,1),coef_tmp(1,1))
end
