program on_top
 implicit none
 BEGIN_DOC
 !
 END_DOC

 read_wf = .true.
 touch read_wf
 ! total one-e integrals 
 io_mo_one_e_integrals = "None"
 touch io_mo_one_e_integrals  
 ! Vne integrals on the MO basis 
 io_mo_integrals_n_e = "None"
 touch io_mo_integrals_n_e
 ! kinetic integrals on the MO basis 
 io_mo_integrals_kinetic = "None"
 touch io_mo_integrals_kinetic 
 ! Vne integrals on the AO basis 
 io_ao_integrals_n_e = "None"
 touch io_ao_integrals_n_e 
 ! kinetic integrals on the AO basis 
 io_ao_integrals_kinetic = "None"
 touch io_ao_integrals_kinetic 

 ! regular 1/r12 integrals  on the MO basis
 io_mo_two_e_integrals = "None"
 touch io_mo_two_e_integrals
 ! regular 1/r12 integrals  on the AO basis
 io_ao_two_e_integrals = "None"
 touch io_ao_two_e_integrals
 ! integral of the effective potential 
 io_mo_int_mu_of_r = "None" 
 touch io_mo_int_mu_of_r 
 call test_compare_on_top
end program

subroutine test_compare_on_top
 implicit none
 BEGIN_DOC
 !
 END_DOC

 integer          :: ipoint, istate, nx
 double precision :: r(3), xmax, dx
 double precision :: rho2, rho2_ex
 double precision :: mu_correction_of_on_top, mu_of_r, f_psi

 nx = 500
 xmax = 4.d0
 dx = xmax/dble(nx)
 r(:) = nucl_coord_transp(:,1)
 r(3) = 0.5d0 * (nucl_coord_transp(3,1) + nucl_coord_transp(3,2) ) ! for the test on h2
 !r(3) = nucl_coord_transp(3,1)                                    ! for the test on he
 r(3) += - xmax * 0.5d0
 write(33,'((A400))')'#       istate     r(3)      mu_of_r     rho2    rho2_e'
 
 do ipoint=1, nx

  do istate=1, N_states
   ! on-top exact
   call give_mu_of_r_cas(r,istate,mu_of_r,f_psi,rho2)

   ! on-top extrapolated
   rho2_ex = mu_correction_of_on_top(mu_of_r,rho2)
  
   !write(33,'(100(F16.10,X))')istate,r(3),mu_of_r,rho2,rho2_ex
   write(33,*)istate,r(3),mu_of_r,rho2,rho2_ex
  enddo 
  r(3) += dx  
 enddo

end
