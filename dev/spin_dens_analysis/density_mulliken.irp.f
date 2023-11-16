

BEGIN_PROVIDER [double precision, density_population, (ao_num,ao_num,N_states)]
 implicit none
 BEGIN_DOC
! density population on the ao basis :
! density_population(i,j) = rho_AO(alpha)(i,j) + rho_AO(beta)(i,j) * <AO_i|AO_j>
 END_DOC
 density_population = 0.d0
 double precision :: accu
 integer :: i,j,istate
 do istate = 1, N_states
  do i = 1, ao_num
   do j = 1, ao_num
    density_population(j,i,istate) = (one_e_dm_ao_alpha_nstates(j,i,istate) + one_e_dm_ao_beta_nstates(j,i,istate)) * ao_overlap(j,i)
   enddo
  enddo
 enddo
END_PROVIDER

BEGIN_PROVIDER [double precision, density_gross_orbital_product, (ao_num,N_states)]
 implicit none
 density_gross_orbital_product = 0.d0
 integer :: i,j,istate
 BEGIN_DOC
! gross orbital product for the density population
 END_DOC
 do istate= 1, N_states
  do i = 1, ao_num
   do j = 1, ao_num
    density_gross_orbital_product(i,istate) += density_population(j,i,N_states)
   enddo
  enddo
 enddo

END_PROVIDER

BEGIN_PROVIDER [double precision, mulliken_density_densities, (nucl_num,N_states)]
 implicit none
 integer :: i,j,istate
 BEGIN_DOC
!ATOMIC density POPULATION (ALPHA MINUS BETA)
 END_DOC
 mulliken_density_densities = 0.d0
 do istate = 1, N_states
  do i = 1, ao_num
   mulliken_density_densities(ao_nucl(i),istate) += density_gross_orbital_product(i,istate)
  enddo
 enddo

END_PROVIDER
