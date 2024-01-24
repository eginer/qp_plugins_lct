
 BEGIN_PROVIDER [double precision, tc_bi_ortho_dipole, (3,N_states)]
 implicit none
 integer :: i,j,istate,m
 double precision :: nuclei_part(3)
 tc_bi_ortho_dipole = 0.d0
 do istate = 1, N_states
  do i = 1, mo_num
   do j = 1, mo_num
    tc_bi_ortho_dipole(1,istate) += -(tc_transition_matrix(j,i,istate,istate)) *  mo_bi_orth_bipole_x(j,i)
    tc_bi_ortho_dipole(2,istate) += -(tc_transition_matrix(j,i,istate,istate)) *  mo_bi_orth_bipole_y(j,i)
    tc_bi_ortho_dipole(3,istate) += -(tc_transition_matrix(j,i,istate,istate)) *  mo_bi_orth_bipole_z(j,i)
   enddo
  enddo
 enddo

 nuclei_part = 0.d0
 do m = 1, 3
  do i = 1,nucl_num
   nuclei_part(m) += nucl_charge(i) * nucl_coord(i,m)
  enddo
 enddo
!
 do istate = 1, N_states
  do m = 1, 3
    tc_bi_ortho_dipole(m,istate) += nuclei_part(m)
  enddo
 enddo
 END_PROVIDER

