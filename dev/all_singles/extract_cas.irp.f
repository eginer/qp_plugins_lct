program get_cas_wf
 implicit none
 read_wf = .True.
 touch read_wf
 call routine_get_cas
end

subroutine routine_get_cas
 implicit none
 integer :: i,j
 double precision :: accu(N_states)
 accu = 0.d0
 do j = 1, N_states
  do i = 1, N_det_cas
   accu(j) += psi_cas_coef(i,j) **2
  enddo
  accu(j) = 1.d0/dsqrt(accu(j))
  do i = 1, N_det_cas
   psi_cas_coef(i,j) *= accu(j) 
  enddo
 enddo
 call save_wavefunction_general(N_det_cas,min(N_states,N_det_cas),psi_cas,size(psi_cas_coef,1),psi_cas_coef)

end
