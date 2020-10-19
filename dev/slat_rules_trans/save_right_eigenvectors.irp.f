program save_right_eigv
 implicit none
 read_wf = .True.
 touch read_wf 
 read_rl_eigv = .True.
 touch read_rl_eigv

 call routine_save
end

subroutine routine_save
implicit none
 call save_wavefunction_general(N_det,min(N_states,N_det),psi_det,size(reigvec_trans,1),reigvec_trans(1,1))
end
