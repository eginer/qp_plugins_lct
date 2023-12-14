program cisd_sc2
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  read_wf = .True.
  touch read_wf 
  call routine
end

subroutine routine
 implicit none
 integer :: it,sze,i,istate
 double precision, allocatable :: u_in(:,:), energies(:), psi_coef_save(:,:)
 double precision :: ecorr,ebefore, eafter, thresh 
 logical :: converged
 external :: H_u_0_nstates_openmp
 sze = N_det
 it = 0
 thresh = 1.d-6
 converged = .False.
 allocate(u_in(sze,N_states_diag),energies(N_states_diag),psi_coef_save(N_det, N_states))
 u_in = 0.d0
 do istate = 1, N_states
   psi_coef_save(1:N_det,istate) = psi_coef(1:N_det,istate)
 enddo
 do i = 1, N_det
  u_in(i,1) = psi_coef_save(i,1)
 enddo

 energies = 0.d0
 ebefore = 0.d0
 print*,'Doing CISDSC2 method: converging ground state'
 do while(.not.converged)
  it += 1
  ecorr = ecorr_tot
  print*,'------------------'
  print*,'it = ',it
  print*,'ecorr = ',ecorr
  call davidson_general_ext_rout_diag_dressed(u_in,h_matrix_diag_all_dets,diag_sc2_dressing,energies,sze,1,N_states_diag,converged,H_u_0_nstates_openmp)           
  print*,'energies = ',energies(1)
  eafter = energies(1)
  converged = dabs(ebefore-eafter).lt.thresh
  ebefore = eafter
  do istate = 1, 1
   do i = 1, N_det
    psi_coef(i,istate)= u_in(i,istate) 
   enddo
  enddo
  touch psi_coef
 enddo
 if(N_states.gt.1)then
  do istate = 2, N_states
    u_in(1:N_det,istate) = psi_coef_save(1:N_det,istate)
  enddo
  call davidson_general_ext_rout_diag_dressed(u_in,h_matrix_diag_all_dets,diag_sc2_dressing,energies,sze,N_states,N_states_diag,converged,H_u_0_nstates_openmp)           
  do i = 2, N_states
   print*,'delta E = ',energies(i) - energies(1)
  enddo
 endif
 do i = 1, N_states
  psi_coef(1:N_det, 1:N_states) = u_in(1:N_det, 1:N_states)
 enddo
 soft_touch psi_coef 
 call save_wavefunction_general(N_det,N_states,psi_det,size(psi_coef,1),psi_coef)
 

end
