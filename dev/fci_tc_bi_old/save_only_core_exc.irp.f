program only_core_exc
 implicit none
 read_wf= .True.
 touch read_wf
 call routine

end

subroutine routine
use bitmasks                                                                                                              
 implicit none
 integer :: number_of_holes,i,igood,ngood
 integer(bit_kind), allocatable :: det_good(:,:,:)
 double precision, allocatable :: coef_good(:,:)
 integer,allocatable :: list_good(:)
 ngood=0
 do i = 1,N_det
  if(number_of_holes(psi_det(1,1,i))==0.and. i.ne.index_hf_psi_det)cycle
  ngood += 1
 enddo
 allocate(det_good(N_int,2,ngood),list_good(ngood),coef_good(ngood,N_states))
 ngood=0
 do i = 1,N_det
  if(number_of_holes(psi_det(1,1,i))==0.and. i.ne.index_hf_psi_det)cycle
  ngood += 1
  det_good(:,:,ngood) = psi_det(:,:,i)
  coef_good(ngood,:) = psi_coef(i,:)
 enddo
 call save_wavefunction_general(ngood,N_states,det_good,ngood,coef_good)

end
