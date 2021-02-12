program density_at_r 
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

 double precision :: r(3), xmax, dx
 double precision :: dm_a(N_states), dm_b(N_states), mu_correction_of_on_top, rho2, rho2_ex, mu_of_r, f_psi
 double precision :: n_at_r
 integer :: istate, ipoint, nx

 istate = 1
 print *, 'r(3)  density_at_r  exact_on_top  extrap_on_top'


 r(1) = nucl_coord(1, 1) ! # of the nucleus in the .xyz file
 r(2) = nucl_coord(1, 2)
 r(3) = nucl_coord(1, 3)
 
 nx = 500
 xmax = 7.d0
 dx = xmax/dble(nx)
 r(3) += -xmax*0.5d0
 
 do ipoint=0, nx
   call dm_dft_alpha_beta_at_r(r, dm_a, dm_b) 
   n_at_r = dm_a(istate) + dm_b(istate)
   call give_mu_of_r_cas(r,istate,mu_of_r,f_psi,rho2) 
   rho2_ex = mu_correction_of_on_top(mu_of_r,rho2)

   print *, r(3), n_at_r, rho2, rho2_ex
   r(3) += dx
 enddo

end program


