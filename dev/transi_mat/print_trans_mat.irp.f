program print_leo
 implicit none
 read_wf = .True.
 touch read_wf 
 call routine

end

subroutine routine
 implicit none
 integer :: i,j,k,NO,NV,istate
 OPEN(unit=77,file='energies.dat',status='unknown',form='formatted')
 do i=1, N_states
    WRITE(77,*) i,  CI_energy(i)
 enddo
 CLOSE(77)

 OPEN(unit=72,file='transmom_GStoEXC.dat',status='unknown',form='formatted')
 do i=2, N_states
   WRITE(72,'(I3,X,3(F16.10,X))') i, trans_dipole_x(1,i), trans_dipole_y(1,i), trans_dipole_z(1,i)
 enddo
 CLOSE(72)

 OPEN(unit=73,file='transmom_EXCntoEXCn.dat',status='unknown',form='formatted')
 do i=2, N_states
   WRITE(73,'(i7,2x,f19.12,2x,f19.12,2x,f19.12)') i, trans_dipole_x(i,i), trans_dipole_y(i,i), trans_dipole_z(i,i)
 enddo
 CLOSE(73)

 OPEN(unit=74,file='transmom_GStoGS.dat',status='unknown',form='formatted')
    WRITE(74,'(i7,2x,f19.12,2x,f19.12,2x,f19.12)') 1, trans_dipole_x(1,1), trans_dipole_y(1,1), trans_dipole_z(1,1)
 CLOSE(74)

 open(unit=75,file='transmom_EXCntoEXCm.dat',status='unknown',form='formatted')
 k = 0
 do i=2, N_states
   do j = i+1, N_states
    k +=1
    WRITE(75,'(i20,2x,f19.12,2x,f19.12,2x,f19.12)') k, trans_dipole_x(i,j), trans_dipole_y(i,j), trans_dipole_z(i,j)
   enddo
 enddo
 CLOSE(75)

 OPEN(unit=91,file='amplitudes.dat',status='unknown',form='formatted')
 NO = elec_alpha_num - n_core_orb
 NV = n_act_orb - NO
 WRITE(91,'(i5)') NO
 ! Original heuristic model (only one escape length value)
 do istate = 2, N_states
  do i = 1, NO
   do j = 1, NV
    write(91,'(f7.4)') cis_amplitudes(i,j,istate)
   enddo  
  enddo
 enddo

 OPEN(unit=50,file='virtual.dat',status='unknown',form='formatted')
 write(50,'(i5)')    NV
 write(50,'(F7.3)')  (Fock_matrix_diag_mo(i),i=elec_alpha_num+1,elec_alpha_num + NV)
 close(50)
 
end
