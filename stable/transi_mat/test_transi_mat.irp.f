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
  integer :: i,j,istate,jstate
  double precision :: accu
  do istate = 1, N_states
   do jstate = 1, N_states
    accu = 0.d0
    do i = 1, mo_num
     do j = 1, mo_num
      accu += transition_matrix(j,i,jstate,istate)  * mo_dipole_z(j,i)
     enddo
    enddo
    print*,'istate,jstate',istate,jstate
    print*,accu,trans_dipole_z(jstate,istate),dabs(accu - trans_dipole_z(jstate,istate))
   enddo
  enddo
end
