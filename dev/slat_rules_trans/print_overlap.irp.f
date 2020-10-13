program print_ovrlp
 implicit none
 read_wf = .True.
 touch read_wf
 read_rl_eigv = .True.
 touch read_rl_eigv
 call routine

end

subroutine routine
 implicit none
 integer :: i,j
 double precision :: accu1,accu2
 double precision, allocatable :: coefs(:,:)
 allocate(coefs(N_det,n_states))
 do i = 1, N_states
  do j = 1, N_det
   coefs(j,i) = psi_coef(j,i)
  enddo
 enddo
 do i = 1, N_states
  accu1 = 0.d0
  do j = 1, N_det
   print*,CI_eigenvectors(j,i),reigvec_trans(j,i),coefs(j,i)
   accu1 += CI_eigenvectors(j,i) * coefs(j,i)
   accu2 += reigvec_trans(j,i)   * coefs(j,i)
  enddo
  print*,'***'
  print*,'***'
  print*,'***'
  print*,'***'
  print*,'accu1 = ',accu1
  print*,'accu2 = ',accu2
 enddo


end
