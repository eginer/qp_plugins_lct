
 BEGIN_PROVIDER [double precision, mu_of_r_prov, (n_points_final_grid,N_states) ]
 implicit none 
 BEGIN_DOC
 ! general variable for mu(r) 
 !
 ! corresponds 
 END_DOC
 integer :: ipoint,istate
 double precision :: cpu0,cpu1
 print*,'providing the mu_of_r ...'
 call wall_time(cpu0)
 do istate = 1, N_states
  do ipoint = 1, n_points_final_grid
   if(mu_of_r_potential.EQ."hf")then
    mu_of_r_prov(ipoint) =  mu_of_r_hf(ipoint)
   else if(mu_of_r_potential.EQ."cas_ful".or.mu_of_r_potential.EQ."cas_truncated")then
    mu_of_r_prov(ipoint) =  mu_of_r_psi_cas_ful(ipoint)
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
 call wall_time(cpu1)
 print*,'Time to provide mu_of_r = ',cpu1-cpu0
 END_PROVIDER 

 BEGIN_PROVIDER [double precision, mu_of_r, (n_points_integration_angular,n_points_radial_grid,nucl_num,N_states) ]
 implicit none
 do istate = 1, N_states
  do ipoint = 1, n_points_final_grid
   k = index_final_points(1,ipoint)
   i = index_final_points(2,ipoint)
   j = index_final_points(3,ipoint)
   mu_of_r(k,i,j) = mu_of_r_vector(ipoint)
  enddo
 enddo
 END_PROVIDER 
