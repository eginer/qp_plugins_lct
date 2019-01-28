program remove_dets
 implicit none
 read_wf = .True.
 touch read_wf
 call routine

end

subroutine routine 
 implicit none
 use bitmasks ! you need to include the bitmasks_module.f90 features
 integer :: i,j
 integer :: n_det_remove,n_det_new,n_det_tmp
 integer, allocatable :: index_det_remove(:)
 integer(bit_kind), allocatable :: psi_new(:,:,:)
 logical, allocatable :: psi_ok(:)
 double precision, allocatable :: psi_coef_new(:,:)
 n_det_remove = 2 
 n_det_new = N_det - n_det_remove 
 allocate(index_det_remove(n_det_remove))

 index_det_remove(1) = 7
 index_det_remove(2) = 8

 allocate(psi_ok(N_det),psi_new(N_int,2,n_det_new),psi_coef_new(n_det_new,N_states))

 psi_ok = .True.
 do i = 1, n_det_remove
  psi_ok(index_det_remove(i)) = .False. 
 enddo

 n_det_tmp = 0
 do i = 1, N_det 
  if(psi_ok(i))then
    n_det_tmp += 1
    do j = 1, N_int
     psi_new(j,1,n_det_tmp) = psi_det(j,1,i)
     psi_new(j,2,n_det_tmp) = psi_det(j,2,i)
    enddo
    do j = 1, N_states
     psi_coef_new(n_det_tmp,j) = psi_coef(i,j)
    enddo
  endif
 enddo

 call save_wavefunction_general(n_det_new,N_states,psi_new,size(psi_coef_new,1),psi_coef_new)

end
