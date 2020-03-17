BEGIN_PROVIDER [ double precision, cis_amplitudes, (n_act_orb,n_act_orb,N_states)]
 implicit none
 integer :: i,j,istate,degree
 integer          :: exc(0:2,2,2)
 double precision :: phase
 integer :: h1,p1,h2,p2,s1,s2
   
   
 
 do i = 1, N_det
  call get_excitation_degree(ref_bitmask,psi_det(1,1,i),degree,N_int) 
  if(degree.gt.1)then
   print*,'This provider is defined for CIS, but there are more than singles ...'
   print*,'stoping'
  endif
 enddo

 do i = 1, N_det
  call get_excitation_degree(ref_bitmask,psi_det(1,1,i),degree,N_int) 
  if(degree.eq.0)cycle
  call get_excitation(ref_bitmask,psi_det(1,1,i),exc,degree,phase,N_int)
  call decode_exc(exc,degree,h1,p1,h2,p2,s1,s2)
  do istate = 1, N_states
   cis_amplitudes(h1,p1,istate) = psi_coef(i,istate)
  enddo
 enddo

END_PROVIDER 


