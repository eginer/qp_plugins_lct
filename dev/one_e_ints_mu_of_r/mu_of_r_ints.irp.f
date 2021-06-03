
 BEGIN_PROVIDER [double precision, mu_of_r_for_ints, (n_points_final_grid,N_states) ]
&BEGIN_PROVIDER [double precision, inv_2_mu_of_r_for_ints, (n_points_final_grid,N_states) ]
&BEGIN_PROVIDER [double precision, inv_4_mu_of_r_for_ints, (n_points_final_grid,N_states) ]
&BEGIN_PROVIDER [double precision, grad_mu_of_r_for_ints, (3,n_points_final_grid,N_states) ]
&BEGIN_PROVIDER [double precision, grad_mu_of_r_transp_for_ints, (n_points_final_grid,N_states,3) ]
&BEGIN_PROVIDER [double precision, grad_sq_mu_of_r_for_ints, (n_points_final_grid,N_states) ]
 implicit none 
 BEGIN_DOC
 !
 ! mu(r) and its gradient for evaluation f integrals for the TC hamiltonian
 END_DOC
 integer :: ipoint,istate,mm
 double precision :: wall0,wall1
 print*,'providing mu_of_r ...'
! PROVIDE mo_two_e_integrals_in_map mo_integrals_map big_array_exchange_integrals 
 call wall_time(wall0)

 if(.not.constant_mu)then
  do istate = 1, N_states
   do ipoint = 1, n_points_final_grid
    if(mu_of_r_tc_ints.EQ."basis")then
     mu_of_r_for_ints(ipoint,istate) =  mu_of_r_hf(ipoint)
    else if(mu_of_r_tc_ints.EQ."rsc")then
     mu_of_r_for_ints(ipoint,istate) =  mu_of_r_rs_c(ipoint,istate)
     grad_mu_of_r_for_ints(:,ipoint,istate) = grad_mu_of_r_rs_c(:,ipoint,istate)
    else if(mu_of_r_tc_ints.EQ."lda")then
     mu_of_r_for_ints(ipoint,istate) =  mu_of_r_lda(ipoint,istate)
     grad_mu_of_r_for_ints(:,ipoint,istate) = grad_mu_of_r_lda(:,ipoint,istate)
    else 
     print*,'you requested the following mu_of_r_tc_ints'
     print*,mu_of_r_tc_ints
     print*,'which does not correspond to any of the options for such keyword'
     stop
    endif

! !!!!!!!!!!!!!!!!!!!!!!!!!!!
!    mu_of_r_for_ints(ipoint,istate) = mu_erf
! !!!!!!!!!!!!!!!!!!!!!!!!!!!
    inv_2_mu_of_r_for_ints(ipoint,istate) = 1.d0/(mu_of_r_for_ints(ipoint,istate))**2
    inv_4_mu_of_r_for_ints(ipoint,istate) = 1.d0/(mu_of_r_for_ints(ipoint,istate))**4

    do mm = 1, 3
     grad_mu_of_r_transp_for_ints(ipoint,istate,mm) = grad_mu_of_r_for_ints(mm,ipoint,istate)
    enddo
!    do mm = 1, 3
!     grad_mu_of_r_transp_for_ints(ipoint,istate,mm) = 0.d0
!    enddo

    grad_sq_mu_of_r_for_ints(ipoint,istate) = 0.d0
    do mm = 1, 3
     grad_sq_mu_of_r_for_ints(ipoint,istate) += grad_mu_of_r_for_ints(mm,ipoint,istate)**2.d0
    enddo
   enddo
  enddo
 else
  do istate = 1, N_states
   do ipoint = 1, n_points_final_grid
    mu_of_r_for_ints(ipoint,istate) =  mu_erf 
    grad_mu_of_r_for_ints(:,ipoint,istate) = 0.d0
    grad_sq_mu_of_r_for_ints(ipoint,istate) = 0.d0
   enddo
  enddo
 endif


 call wall_time(wall1)
 print*,'Time to provide mu_of_r_for_ints = ',wall1-wall0
 END_PROVIDER 


 BEGIN_PROVIDER [double precision, mu_of_r_extra_grid_for_ints, (n_points_extra_final_grid,N_states) ]
 implicit none 
 BEGIN_DOC
 !
 ! mu(r) and its gradient for evaluation f integrals for the TC hamiltonian
 END_DOC
 integer :: ipoint,istate,mm
 double precision :: wall0,wall1
 print*,'providing mu_of_r_extra_grid ...'
! PROVIDE mo_two_e_integrals_in_map mo_integrals_map big_array_exchange_integrals 
 call wall_time(wall0)

 if(.not.constant_mu)then

  do istate = 1, N_states
   do ipoint = 1, n_points_extra_final_grid
    if(mu_of_r_tc_ints.EQ."basis")then
     mu_of_r_extra_grid_for_ints(ipoint,istate) =  mu_of_r_extra_grid_hf(ipoint)
    else if(mu_of_r_tc_ints.EQ."rsc")then
     mu_of_r_extra_grid_for_ints(ipoint,istate) =  mu_of_r_extra_grid_rs_c(ipoint,istate)
    else if(mu_of_r_tc_ints.EQ."lda")then
     mu_of_r_extra_grid_for_ints(ipoint,istate) =  mu_of_r_extra_grid_lda(ipoint,istate)
    else 
     print*,'you requested the following mu_of_r_extra_grid_tc_ints'
     print*,mu_of_r_tc_ints
     print*,'which does not correspond to any of the options for such keyword'
     stop
    endif
!    mu_of_r_extra_grid_for_ints(ipoint,istate) = mu_erf
   enddo
  enddo
 else
  do istate = 1, N_states
   do ipoint = 1, n_points_extra_final_grid
    mu_of_r_extra_grid_for_ints(ipoint,istate) =  mu_erf 
   enddo
  enddo
 endif


 call wall_time(wall1)
 print*,'Time to provide mu_of_r_extra_grid_for_ints = ',wall1-wall0
 END_PROVIDER 


