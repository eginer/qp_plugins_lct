BEGIN_PROVIDER[double precision, gamma_i_j, (mo_num, mo_num, N_states)]
 implicit none
 BEGIN_DOC
 ! sum_I(1 -> n_det) sum_J(1 -> n_det) c_I c_J <phi_I | E_ij | phi_J>
 END_DOC 
 integer :: i, j, i_state
 integer :: ii, jj
 integer :: det_I, det_J
 double precision :: phase
 
 gamma_i_j = 0.d0

 do i_state=1, N_states
  do i = 1, n_core_orb 
   ii=list_core(i)
   gamma_i_j(ii,ii,i_state) = 2.d0
  enddo
 enddo


do i_state=1, N_states
 do i = 1, n_act_orb 
  ii = list_act(i)
  do j = 1, n_act_orb 
   jj = list_act(j)
   do det_I = 1, N_det
     det_J = Aindx_mpq(det_I,i,j,1) ! <det_J | E_ij | det_I>
     if (det_J.eq.-1) cycle 
     phase=Aval_mpq(det_I,i,j,1)
     gamma_i_j(ii,jj,i_state) = gamma_i_j(ii,jj,i_state) + psi_coef(det_I,i_state)*psi_coef(det_J,i_state)*phase 
     if (secnd_Aval_mpq(det_I,i,j)) then
      det_J=Aindx_mpq(det_I,i,j,2)
      phase=Aval_mpq(det_I,i,j,2)
      gamma_i_j(ii,jj,i_state) = gamma_i_j(ii,jj,i_state) + psi_coef(det_I,i_state)*psi_coef(det_J,i_state)*phase 
     end if
   enddo
  enddo
 enddo
enddo

!do i=1, mo_num
! print*, gamma_i_j(i,:,1)
!enddo
END_PROVIDER

 BEGIN_PROVIDER[double precision, gamma_I_i_j, (N_det, mo_num, mo_num, N_states)]
 implicit none
 BEGIN_DOC
 ! sum_J(1 -> n_det) c_J <phi_I | E_ij | phi_J>
 END_DOC 
 integer :: i, j, i_state
 integer :: ii, jj
 integer :: det_I, det_J
 double precision :: phase
 
 gamma_I_i_j = 0.d0

 do i_state=1, N_states
   do i = 1, n_core_orb 
    ii=list_core(i)
    do det_I=1, N_det
     det_J = Aindx_mpq(det_I,i,i,1) ! <det_J | E_ij | det_I>
     if (det_J.eq.-1) cycle
     gamma_I_i_j(det_I,ii,ii,i_state) = gamma_I_i_j(det_I,ii,ii,i_state) + 2.d0*psi_coef(det_J,i_state)
   enddo
  enddo
 enddo

 do i_state=1, N_states
  do i = 1, n_act_orb
   ii = list_act(i)
   do j = 1, n_act_orb 
    jj = list_act(j)
    do det_I = 1, N_det
     det_J = Aindx_mpq(det_I,i,j,1) ! <det_J | E_ij | det_I>
     if (det_J.eq.-1) cycle 
     phase=Aval_mpq(det_I,i,j,1)
     gamma_I_i_j(det_I,jj,ii,i_state) = gamma_I_i_j(det_I,jj,ii,i_state) +  psi_coef(det_J,i_state)* phase 
     if (secnd_Aval_mpq(det_I,i,j)) then
      det_J=Aindx_mpq(det_I,i,j,2)
      phase=Aval_mpq(det_I,i,j,2)
      gamma_I_i_j(det_I,jj,ii,i_state) = gamma_I_i_j(det_I,jj,ii,i_state) +  psi_coef(det_J,i_state)* phase 
     end if
    enddo
   enddo
  enddo
 enddo


 END_PROVIDER

 BEGIN_PROVIDER[double precision, delta_gamma_I_i_j, (N_det, mo_num, mo_num, N_states)]
 implicit none
 BEGIN_DOC
 ! sum_J(1 -> n_det) c_J <phi_I | E_ij | phi_J>
 END_DOC 
 integer :: i, j, i_state
 integer :: det_I, det_J
 double precision :: phase
 
 delta_gamma_I_i_j = 0.d0
do i_state=1, N_states
 do i = 1, mo_num
  do j = 1, mo_num

   do det_I = 1, N_det
    delta_gamma_I_i_j(det_I, j, i, i_state) = gamma_I_i_j(det_I,j,i,i_state) - psi_coef(det_I,i_state)* gamma_i_j(j,i,i_state) 
   enddo
  enddo
 enddo
enddo


 END_PROVIDER

