

BEGIN_PROVIDER [double precision, density_population, (ao_num,ao_num)]
 implicit none
 integer :: i,j
 BEGIN_DOC
! density population on the ao basis :
! density_population(i,j) = rho_AO(alpha)(i,j) + rho_AO(beta)(i,j) * <AO_i|AO_j>
 END_DOC
 density_population = 0.d0
 do i = 1, ao_num
  do j = 1, ao_num
   density_population(j,i) = (one_e_dm_ao_alpha(j,i) + one_e_dm_ao_beta(j,i)) * ao_overlap(j,i)
  enddo
 enddo
END_PROVIDER

BEGIN_PROVIDER [double precision, density_gross_orbital_product, (ao_num)]
 implicit none
 density_gross_orbital_product = 0.d0
 integer :: i,j
 BEGIN_DOC
! gross orbital product for the density population
 END_DOC
 do i = 1, ao_num
  do j = 1, ao_num
   density_gross_orbital_product(i) += density_population(j,i)
  enddo
 enddo

END_PROVIDER

BEGIN_PROVIDER [double precision, mulliken_density_densities, (nucl_num)]
 implicit none
 integer :: i,j
 BEGIN_DOC
!ATOMIC density POPULATION (ALPHA MINUS BETA)
 END_DOC
 mulliken_density_densities = 0.d0
 do i = 1, ao_num
  mulliken_density_densities(ao_nucl(i)) += density_gross_orbital_product(i)
 enddo

END_PROVIDER
