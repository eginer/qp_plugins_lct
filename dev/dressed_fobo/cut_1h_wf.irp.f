program remove_small_1h
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  read_wf = .true.
  touch read_wf
  call routine_print_1h_and_cas
  pause
  call routine_cut_1h_and_cas
end

subroutine routine_cut_1h
 implicit none
 integer :: i
 use bitmasks ! you need to include the bitmasks_module.f90 features
 integer(bit_kind), allocatable :: psi_new(:,:,:)
 double precision, allocatable :: psi_coef_new(:,:)
 integer, allocatable :: index_good(:)
 allocate(index_good(N_det))
 integer :: n_det_new
 logical :: is_a_1h
 thresh_coef_lmct = 1.d-2
 integer :: n_1h_tot,n_1h_small
 n_1h_tot = 0
 n_1h_small = 0
 n_det_new = 0
 do i = 1, N_det
  if(is_a_1h(psi_det(1,1,i)))then
   n_1h_tot += 1
   if(dabs(psi_coef(i,1)/psi_coef(1,1)).lt.thresh_coef_lmct)then
    n_1h_small += 1
    cycle
   else 
    print*,''
    print*,''
    call debug_det(psi_det(1,1,i),N_int)
    print*,'dabs(psi_coef(i,1)/psi_coef(1,1))',dabs(psi_coef(i,1)/psi_coef(1,1))
   endif
  endif
  n_det_new += 1
  index_good(n_det_new) = i
 enddo
 print*,'*******************'
 print*,'*******************'
 print*,'N_det     = ',N_det
 print*,'n_det_new = ',n_det_new
 print*,'n_1h_tot  = ',n_1h_tot
 print*,'n_1h_small= ',n_1h_small
 print*,'*******************'
 print*,'*******************'
 print*,'n_1h_big  = ',n_1h_tot - n_1h_small
 allocate(psi_new(N_int,2,n_det_new),psi_coef_new(n_det_new,n_states))
 integer :: istate 
 integer :: k
 do i = 1, n_det_new
  do k = 1, N_int 
   psi_new(k,1,i) = psi_det(k,1,index_good(i))
   psi_new(k,2,i) = psi_det(k,2,index_good(i))
  enddo
 enddo

 do istate = 1, N_states
  do i = 1, n_det_new
   psi_coef_new(i,istate) = psi_coef(index_good(i),istate)
  enddo
 enddo

 call save_wavefunction_general(N_det_new,N_states,psi_new,size(psi_coef_new,1),psi_coef_new)

end

subroutine routine_cut_1h_and_cas
 implicit none
 integer :: i
 use bitmasks ! you need to include the bitmasks_module.f90 features
 integer(bit_kind), allocatable :: psi_new(:,:,:)
 double precision, allocatable :: psi_coef_new(:,:)
 integer, allocatable :: index_good(:)
 allocate(index_good(N_det))
 integer :: n_det_new
 logical :: is_a_1h
 integer :: number_of_holes,nh
 integer :: number_of_particles,np
 n_det_new = 0
 do i = 1, N_det
  nh = number_of_holes(psi_det(1,1,i))
  np = number_of_particles(psi_det(1,1,i))
  if(is_a_1h(psi_det(1,1,i)))then
   if(dabs(psi_coef(i,1)/psi_coef(1,1)).gt.thresh_coef_lmct)then
    n_det_new += 1
    index_good(n_det_new) = i
   endif
  endif
  if(nh==0.and.np==0)then
   n_det_new += 1
   index_good(n_det_new) = i
  endif
 enddo
 print*,'*******************'
 print*,'*******************'
 print*,'N_det     = ',N_det
 print*,'n_det_new = ',n_det_new
 print*,'*******************'
 print*,'*******************'
 allocate(psi_new(N_int,2,n_det_new),psi_coef_new(n_det_new,n_states))
 integer :: istate 
 integer :: k
 do i = 1, n_det_new
  do k = 1, N_int 
   psi_new(k,1,i) = psi_det(k,1,index_good(i))
   psi_new(k,2,i) = psi_det(k,2,index_good(i))
  enddo
 enddo

 do istate = 1, N_states
  do i = 1, n_det_new
   psi_coef_new(i,istate) = psi_coef(index_good(i),istate)
  enddo
 enddo

 call save_wavefunction_general(N_det_new,N_states,psi_new,size(psi_coef_new,1),psi_coef_new)

end

subroutine routine_print_1h_and_cas
 implicit none
 integer :: i
 use bitmasks ! you need to include the bitmasks_module.f90 features
 integer(bit_kind), allocatable :: psi_new(:,:,:)
 double precision, allocatable :: psi_coef_new(:,:)
 integer, allocatable :: index_good(:)
 allocate(index_good(N_det))
 integer :: n_det_new
 logical :: is_a_1h
 integer :: number_of_holes,nh
 integer :: number_of_particles,np
 n_det_new = 0
 do i = 1, N_det
  nh = number_of_holes(psi_det(1,1,i))
  np = number_of_particles(psi_det(1,1,i))
  if(is_a_1h(psi_det(1,1,i)))then
   if(dabs(psi_coef(i,1)/psi_coef(1,1)).gt.thresh_coef_lmct)then
    n_det_new += 1
    index_good(n_det_new) = i
   endif
  endif
  if(nh==0.and.np==0)then
   n_det_new += 1
   index_good(n_det_new) = i
  endif
 enddo
 print*,'*******************'
 print*,'*******************'
 print*,'N_det     = ',N_det
 print*,'n_det_new = ',n_det_new
 print*,'*******************'
 print*,'*******************'

end


