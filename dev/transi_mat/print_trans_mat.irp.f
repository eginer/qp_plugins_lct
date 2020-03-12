program print_leo
 implicit none
 read_wf = .True.
 touch read_wf 
 call routine

end

subroutine routine
 implicit none
 integer :: i,j,k
 OPEN(unit=77,file='energies.dat',status='unknown',form='formatted')
 DO i=1, N_states
    WRITE(77,*) i,  CI_energy(i)
 ENDDO
 CLOSE(77)

 OPEN(unit=72,file='transmom_GStoEXC.dat',status='unknown',form='formatted')
 DO i=2, N_states
   WRITE(72,'(I3,X,3(F16.10,X))') i, trans_dipole_x(1,i), trans_dipole_y(1,i), trans_dipole_z(1,i)
 ENDDO
 CLOSE(72)

 OPEN(unit=73,file='transmom_EXCntoEXCn.dat',status='unknown',form='formatted')
 DO i=2, N_states
   WRITE(73,'(i7,2x,f19.12,2x,f19.12,2x,f19.12)') i, trans_dipole_x(i,i), trans_dipole_y(i,i), trans_dipole_z(i,i)
 ENDDO
 CLOSE(73)

 OPEN(unit=74,file='transmom_GStoGS.dat',status='unknown',form='formatted')
    WRITE(74,'(i7,2x,f19.12,2x,f19.12,2x,f19.12)') 1, trans_dipole_x(1,1), trans_dipole_y(1,1), trans_dipole_z(1,1)
 CLOSE(74)

 open(unit=75,file='transmom_EXCntoEXCm.dat',status='unknown',form='formatted')
 k = 0
 DO i=2, N_states
   do j = i+1, N_states
    k +=1
    WRITE(75,'(i20,2x,f19.12,2x,f19.12,2x,f19.12)') k, trans_dipole_x(i,j), trans_dipole_y(i,j), trans_dipole_z(i,j)
   enddo
 ENDDO
 CLOSE(75)

end
