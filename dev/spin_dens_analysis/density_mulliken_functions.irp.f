
subroutine density_mulliken_density_mat(density_mat, mull_dens_atoms)
 implicit none
 integer :: i
 double precision, intent(in) :: density_mat(ao_num, ao_num)
 double precision, intent(out):: mull_dens_atoms(nucl_num)

 double precision :: density_gross_orb_prod(ao_num), density_pop(ao_num, ao_num)

! print*,'density_mat'
! print*, density_mat 
 mull_dens_atoms = 0.d0
 call get_dens_pop(density_mat,density_pop)
 call give_dens_orbital_prod(density_pop,density_gross_orb_prod)
 do i = 1, ao_num
  mull_dens_atoms(ao_nucl(i)) += density_gross_orb_prod(i)
 enddo


end

subroutine  give_dens_orbital_prod(density_pop,density_gross_orb_prod)
 implicit none
 double precision, intent(in) :: density_pop(ao_num, ao_num)
 double precision, intent(out) :: density_gross_orb_prod(ao_num)
 integer :: i,j
 BEGIN_DOC
! gross orbital product for the density population
 END_DOC
 density_gross_orb_prod = 0.d0
 do i = 1, ao_num
  do j = 1, ao_num
   density_gross_orb_prod(i) += density_pop(j,i)
  enddo
 enddo
end

subroutine get_dens_pop(density_mat,density_pop)
 implicit none
 double precision, intent(in) :: density_mat(ao_num, ao_num)
 double precision, intent(out):: density_pop(ao_num, ao_num)
 integer :: i,j
 BEGIN_DOC
! density population on the ao basis :
! density_population(i,j) = rho_AO(alpha)(i,j) + rho_AO(beta)(i,j) * <AO_i|AO_j>
 END_DOC
 density_pop = 0.d0
 do i = 1, ao_num
  do j = 1, ao_num
!   print*,i,j,density_mat(j,i) , ao_overlap(j,i)
   density_pop(j,i) = density_mat(j,i) * ao_overlap(j,i)
  enddo
 enddo

end
