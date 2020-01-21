program test
 implicit none
  read_wf = .True.
  touch read_wf
  call routine
end

subroutine routine
 implicit none
 integer :: i,j,degree,nsingles,idx_hf
 integer, allocatable :: idx(:)
 double precision, allocatable :: coef(:)
 double precision :: hij,hjj,eref,oij
 nsingles = 0
 eref = ref_bitmask_energy
 do i = 1, N_det
  call get_excitation_degree(ref_bitmask,psi_det(1,1,i),degree,N_int)
  if(degree == 0)then
   psi_coef(i,1) = 1.d0
  else 
   call i_H_j(psi_det(1,1,i),psi_det(1,1,i),N_int,hjj)
   call i_H_j(ref_bitmask,psi_det(1,1,i),N_int,hij)
   psi_coef(i,1) = hij/(eref - hjj)
  endif
 enddo
 call save_wavefunction_general(N_det,N_states,psi_det,size(psi_coef,1),psi_coef)
end
