
 BEGIN_PROVIDER [double precision, mu_of_r_prov, (n_points_final_grid,N_states) ]
 implicit none 
 BEGIN_DOC
 ! general variable for mu(r) 
 !
 ! corresponds to Eq. (37) of J. Chem. Phys. 149, 194301 (2018) 
 !
 ! !!!!!! WARNING !!!!!! if no_core_density == .True. then all contributions from the core orbitals 
 !
 ! in the two-body density matrix are excluded
 END_DOC
 integer :: ipoint,istate
 double precision :: wall0,wall1
 print*,'providing mu_of_r ...'
 call wall_time(wall0)
 do istate = 1, N_states
  do ipoint = 1, n_points_final_grid
   if(mu_of_r_potential.EQ."hf")then
    mu_of_r_prov(ipoint) =  mu_of_r_hf(ipoint)
   else if(mu_of_r_potential.EQ."cas_ful".or.mu_of_r_potential.EQ."cas_truncated")then
    mu_of_r_prov(ipoint) =  mu_of_r_psi_cas(ipoint)
   else if(mu_of_r_potential.EQ."Read")then
    mu_of_r_prov(ipoint) =  mu_of_r_array(ipoint)
   else 
     print*,'you requested the following mu_of_r_potential'
     print*,mu_of_r_potential
     print*,'which does not correspond to any of the options for such keyword'
     stop
   endif
  enddo
 enddo
 call wall_time(wall1)
 print*,'Time to provide mu_of_r = ',wall1-wall0
 END_PROVIDER 


 BEGIN_PROVIDER [double precision, mu_of_r_hf, (n_points_final_grid) ]
 implicit none 
 BEGIN_DOC
 ! mu(r) computed with a HF wave function (assumes that HF MOs are stored in the EZFIO)
 !
 ! corresponds to Eq. (37) of J. Chem. Phys. 149, 194301 (2018) but for \Psi^B = HF^B
 !
 ! !!!!!! WARNING !!!!!! if no_core_density == .True. then all contributions from the core orbitals 
 !
 ! in the two-body density matrix are excluded
 END_DOC
 integer :: ipoint
 double precision :: wall0,wall1,f_hf,on_top,w_hf,pi
 print*,'providing mu_of_r_hf ...'
 call wall_time(wall0)
 pi = dacos(-1.d0)
 provide f_psi_hf_ab 
 !$OMP PARALLEL DO &
 !$OMP DEFAULT (NONE)  &
 !$OMP PRIVATE (ipoint,f_hf,on_top,w_hf) & 
 !$OMP ShARED (n_points_final_grid,mu_of_r_hf,f_psi_hf_ab,on_top_psi_hf,pi) 
 do ipoint = 1, n_points_final_grid
  f_hf   = f_psi_hf_ab(ipoint)
  on_top = on_top_psi_hf(ipoint)
  if(on_top.le.1.d-12.or.f_hf.le.0.d0.or.f_hf * on_top.lt.0.d0)then
    w_hf   = 1.d+10
  else 
    w_hf  = f_hf /  on_top
  endif
  mu_of_r_hf(ipoint) =  w_hf * pi * 0.5d0
 enddo
 !$OMP END PARALLEL DO
 call wall_time(wall1)
 print*,'Time to provide mu_of_r_hf = ',wall1-wall0
 END_PROVIDER 

 BEGIN_PROVIDER [double precision, mu_of_r_psi_cas, (n_points_final_grid,N_states) ]
 implicit none 
 BEGIN_DOC
 ! mu(r) computed with a wave function developped in an active space
 !
 ! corresponds to Eq. (37) of J. Chem. Phys. 149, 194301 (2018) 
 !
 ! !!!!!! WARNING !!!!!! if no_core_density == .True. then all contributions from the core orbitals 
 !
 ! in the two-body density matrix are excluded
 END_DOC
 integer :: ipoint
 double precision :: wall0,wall1,f_psi,on_top,w_psi,pi
 print*,'providing mu_of_r_psi_cas ...'
 call wall_time(wall0)
 pi = dacos(-1.d0)

 provide total_cas_on_top_density f_psi_cas_ab
 !$OMP PARALLEL DO &
 !$OMP DEFAULT (NONE)  &
 !$OMP PRIVATE (ipoint,f_psi,on_top,w_psi) & 
 !$OMP ShARED (n_points_final_grid,mu_of_r_hf,f_psi_hf_ab,on_top_psi_hf,pi) 
 do istate = 1, N_states
  do ipoint = 1, n_points_final_grid
   f_psi  = f_psi_cas_ab(ipoint,istate) 
   on_top = total_cas_on_top_density(ipoint,istate)
   if(on_top.le.1.d-12.or.f_psi.le.0.d0.or.f_psi * on_top.lt.0.d0)then
     w_psi   = 1.d+10
   else 
     w_psi  = f_psi /  on_top
   endif
   mu_of_r_psi_cas(ipoint) = w_psi * pi * 0.5d0
  enddo
 enddo
 !$OMP END PARALLEL DO
 call wall_time(wall1)
 print*,'Time to provide mu_of_r_psi_cas = ',wall1-wall0
 END_PROVIDER 

