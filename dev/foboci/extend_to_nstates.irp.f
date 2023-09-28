program extend_to_nstates
 implicit none
 read_wf = .True.
 touch read_wf
 call routine

end

subroutine routine
 implicit none
 print*,'**************************'
 print*,'Extending the wave function to this number of states'
 print*,'n_states_extend   === ',n_states_extend
 use bitmasks ! you need to include the bitmasks_module.f90 features
 double precision, allocatable :: psi_coef_new(:,:)
 integer :: dim_psicoef
 dim_psicoef = N_det
 allocate(psi_coef_new(dim_psicoef,n_states_extend))
 integer :: i,j
 psi_coef_new = 0.d0
 do i = 1, N_det
  psi_coef_new(i,1) = psi_coef(i,1)
 enddo
 do i = 2, n_states_extend
  psi_coef_new(i,i) = 1.d0
 enddo
 call save_wavefunction_general(n_det,n_states_extend,psi_det,dim_psicoef,psi_coef_new)
 

end
