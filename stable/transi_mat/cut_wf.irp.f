program cut_wf
 implicit none
 read_wf = .True.
 touch read_wf
 call routine
end

subroutine routine
 implicit none
 double precision :: thr, accu
 integer :: i,j
 thr = 1.d-2
 do i = 1, N_states
  accu = 0.d0
  do j = 1, N_det
   if(dabs(psi_coef(j,i)).lt.thr)then
    psi_coef(j,i) = 0.d0
   endif
   accu += psi_coef(j,i)**2.d0
  enddo
  accu = 1.d0/dsqrt(accu)
  do j = 1, N_det
   psi_coef(j,i) = psi_coef(j,i) * accu 
  enddo
 enddo
 touch psi_coef
 call save_wavefunction_general(N_det,N_states,psi_det,size(psi_coef,1),psi_coef)
end
