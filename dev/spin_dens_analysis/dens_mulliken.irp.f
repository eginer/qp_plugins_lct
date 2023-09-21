program density_mulliken
 implicit none
  read_wf = .True.
  touch read_wf
  call routine 

end 

subroutine routine
 implicit none
 double precision :: accu
 integer :: i
 integer :: j
 print*,'Mulliken density analysis '
 accu= 0.d0
 do i = 1, nucl_num
  print*,i,nucl_charge(i)-mulliken_density_densities(i),nucl_coord(i,3)
  accu += nucl_charge(i)-mulliken_density_densities(i)
 enddo
 double precision, allocatable :: mull_dens_atoms(:)
 allocate(mull_dens_atoms(nucl_num))

! call density_mulliken_density_mat(one_e_dm_ao, mull_dens_atoms) 
! print*,'test '
! do i = 1, nucl_num
!  print*,i,mulliken_density_densities(i),mull_dens_atoms(i), dabs(mulliken_density_densities(i)-mull_dens_atoms(i))
! enddo
end
