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
 integer :: i,j,i_state, det_I
 double precision :: contrib
 double precision, allocatable :: gamma_i_j_test(:,:,:)

do i_state=1, N_states
 do i = 1, mo_num
  do j = 1, mo_num
   contrib = dabs(gamma_i_j(j,i,i_state) - (one_e_dm_mo_alpha(j,i,i_state) + one_e_dm_mo_beta(j,i,i_state)))
   if (contrib.gt.1E-10) then
    print*, "j=", j," ;i=",i 
    print*, gamma_i_j(j,i,i_state), (one_e_dm_mo_alpha(j,i,i_state) + one_e_dm_mo_beta(j,i,i_state)), contrib
   endif
  enddo
 enddo
enddo

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

end subroutine
