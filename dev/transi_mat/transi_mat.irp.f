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
  provide trans_dipole_z
print*,'Matrix of z dipole '
  do i = 1, N_states
   write(*,'(100(F10.5,X))')trans_dipole_z(i,:)
!   do j = i, N_states
!    print*,i,j,dabs(trans_dipole_z_bourrin(j,i)-trans_dipole_z(j,i))
!   enddo
  enddo



print*,''
print*,''
print*,''
print*,''
print*,'                    Transition Dipole Moments and            '
print*,'            GMH Couplings Between Ground and Singlet Excited States       '
print*,' -------------------------------------------------------------------------------- '
print*,'    States   X          Y          Z(a.u.)     Coupling(eV)  '
print*,' -------------------------------------------------------------------------------- '
  do i = 1, 1
   do j = 2, N_states
    write(*,'(2(I3,X),10(F16.10,X))')i-1,j-1,trans_dipole_x(i,j),trans_dipole_y(i,j),trans_dipole_z(i,j) , oscillator_strength(j)
   enddo
  enddo

print*,''

print*,'                    Transition Dipole Moments and            '
print*,'            GMH Couplings Between Excited and Singlet Excited States       '

  print*,' --------------------------------------------------------------------------------'
  print*,'States   X          Y          Z(a.u.)     Coupling(eV)'
  print*,' --------------------------------------------------------------------------------'
  do i = 2, N_states
   do j = 2, N_states
    write(*,'(2(I3,X),3(F16.10,X))')i-1,j-1,trans_dipole_x(i,j),trans_dipole_y(i,j),trans_dipole_z(i,j)   
   enddo
  enddo

  print*,' --------------------------------------------------------------------------------'
  print*,' --------------------------------------------------------------------------------'
  print*,' --------------------------------------------------------------------------------'
  print*,' --------------------------------------------------------------------------------'
  print*,'States   X          Y          Z(a.u.)     Coupling(eV)'
  print*,' --------------------------------------------------------------------------------'
  do i = 2, N_states
   do j = 2, N_states
    write(*,'(2(I3,X),3(F16.10,X))')i-1,j-1,trans_dipole_x_bourrin(i,j),trans_dipole_y_bourrin(i,j),trans_dipole_z_bourrin(i,j)   
   enddo
  enddo
end
