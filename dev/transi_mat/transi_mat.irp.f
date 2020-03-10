program transi_mat
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  print *, 'Hello world'
  read_wf = .True.
  touch read_wf
  call routine

end

subroutine routine
  implicit none
  integer :: i,j
  do i = 1, N_states
   do j = i, N_states
    print*,i,j,dabs(trans_dipole_z_bourrin(j,i)-trans_dipole_z(j,i))
   enddo
  enddo
end
