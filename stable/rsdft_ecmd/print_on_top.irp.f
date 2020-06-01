program print_on_top
  implicit none
  BEGIN_DOC
  !
  END_DOC
  read_wf = .True.
  touch read_wf

  io_mo_one_e_integrals = "None"
  touch io_mo_one_e_integrals  
  io_mo_two_e_integrals = "None"
  touch io_mo_two_e_integrals
  io_ao_two_e_integrals = "None"
  touch io_ao_two_e_integrals
  io_mo_two_e_integrals_erf = "None" 
  touch io_mo_two_e_integrals_erf
  io_ao_two_e_integrals_erf = "None" 
  touch io_ao_two_e_integrals_erf
 
  io_mo_integrals_n_e = "None"
  touch io_mo_integrals_n_e
  io_mo_integrals_kinetic = "None"
  touch io_mo_integrals_kinetic 
  io_ao_integrals_n_e = "None"
  touch io_ao_integrals_n_e 
  io_ao_integrals_kinetic = "None"
  touch io_ao_integrals_kinetic 
 
  call print_on_top_routine
  
end

subroutine print_on_top_routine
 implicit none
 BEGIN_DOC
 !
 END_DOC

 integer          :: nx, istate, ipoint
 double precision :: xmax, r(3), weight, dx
 double precision :: rho_a, rho_b, rho, g0, dg0drho
 double precision :: on_top, on_top_extrap, on_top_UEG, on_top_xc, on_top_extrap_xc, on_top_UEG_xc
 double precision :: mu_correction_of_on_top, mu_of_r, f_psi
  nx = 500
  istate = 1
  xmax = 4.d0
  dx = xmax/dble(nx)
  r(:) = nucl_coord_transp(:,1)
  r(3) += -xmax*0.5d0
  write(34,*)  '#r(3)   on_top   on_top_extrap   on_top_UEG   on_top_xc  on_top_extrap_xc   ontop_UEG_xc    rho'
   do ipoint = 1, nx
    
    mu_of_r = mu_erf_dft
   ! call give_mu_of_r_cas(r,istate,mu_of_r,f_psi,on_top)

    call dm_dft_alpha_beta_at_r(r,rho_a,rho_b)
    rho   = rho_a + rho_b
    call g0_dg0(rho, rho_a, rho_b, g0, dg0drho)
   
    call give_on_top_in_r_one_state(r,istate,on_top)
    on_top = on_top*2.d0
    on_top_extrap = 2.d0 * mu_correction_of_on_top(mu_of_r,on_top)
    on_top_UEG    = (rho**2)*g0
   
    on_top_xc        = on_top        - rho**2
    on_top_extrap_xc = on_top_extrap - rho**2 
    on_top_UEG_xc    = on_top_UEG    - rho**2

    write(34,*) r(3), on_top, on_top_extrap, on_top_UEG, on_top_xc, on_top_extrap_xc, on_top_UEG_xc, rho
    
    r(3) += dx
   enddo
end subroutine
