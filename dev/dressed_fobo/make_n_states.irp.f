program make_n_states
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  read_wf = .true.
  touch read_wf
  call routine_make_n_states
end

subroutine routine_make_n_states
 implicit none
 use bitmasks ! you need to include the bitmasks_module.f90 features
 double precision, allocatable :: psi_coef_new(:,:)
 integer :: n_states_new,istate,i
 write(*,*)'how many states do you want ?'
 read(5,*)n_states_new
 allocate(psi_coef_new(N_det,n_states_new))

 istate = 1
  do i = 1, n_det
   psi_coef_new(i,istate) = psi_coef(i,istate)
  enddo
 double precision :: dnorm 
 dnorm = 1.d0/dsqrt(dble(N_det))
 do istate = 2, n_states_new
  do i = 1, n_det
   psi_coef_new(i,istate) = dnorm
  enddo
 enddo

 call save_wavefunction_general(N_det,n_states_new,psi_det,size(psi_coef_new,1),psi_coef_new)

end
