program test_linear_response
  implicit none
  BEGIN_DOC
  END_DOC

  read_wf=.true.
  touch read_wf

  call routine
end

subroutine routine
  implicit none
  BEGIN_DOC
  END_DOC
 integer :: i,j,i_state, det_I, det_J
 integer :: det_II
 integer :: k,l
 double precision :: contrib
 double precision, allocatable :: gamma_i_j_test(:,:,:)

print*, "Test gamma_i_j VS one e dm mo"
do i_state=1, N_states
 do i = 1, mo_num
  do j = 1, mo_num
   contrib = dabs(gamma_i_j(i,j,i_state) - (one_e_dm_mo_alpha(i,j,i_state) + one_e_dm_mo_beta(i,j,i_state)))
   if (contrib.gt.1E-10) then
    print*, "j=", j," ;i=",i, " ;istate=", i_state 
    print*, gamma_i_j(i,j,i_state), (one_e_dm_mo_alpha(i,j,i_state) + one_e_dm_mo_beta(i,j,i_state)), contrib
   endif
  enddo
 enddo
enddo


print*, "Test gamma_I_i_j"
allocate(gamma_i_j_test(mo_num, mo_num,N_states))
gamma_i_j_test =0.d0
do i_state=1, N_states
 do i = 1, mo_num
  do j = 1, mo_num

   do det_I=1, N_det
    gamma_i_j_test(j, i, i_state) += psi_coef(det_I,i_state)*gamma_I_i_j(det_I, j, i, i_state)
   enddo

  contrib = gamma_i_j_test(j, i, i_state) - gamma_i_j(j,i,i_state)
   if (contrib.gt.1E-10) then
    print*, "j=", j," ;i=",i 
    print*, gamma_i_j_test(j, i, i_state), gamma_i_j(j,i,i_state), contrib
   endif
  enddo
 enddo
enddo

print*, "Test kernel cis"
print*, "Print kernel_lr prov for cis:"
do det_II=1, N_det-1
 write(*,'(100(f10.5,x))') K_test_cis(1:N_det-1,det_II)
enddo
print*, "Print kernel_lr prov:"
do det_I=1, N_det
 write(*,'(100(f10.5,x))') kernel_lr(1:N_det,det_I,1)
enddo

!do det_I=1, N_det
! do det_J=1, N_det
!  print*, "A(", det_I, ",", det_J, ")=", A_IJ(det_I, det_J,1)
!  print*, "B(", det_I, ",", det_J, ")=", B_IJ(det_I, det_J,1)
!  print*, "S(", det_I, ",", det_J, ")=", S_IJ(det_I, det_J,1)
! enddo
!enddo
end subroutine
