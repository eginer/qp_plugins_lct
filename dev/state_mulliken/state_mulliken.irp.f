program state_mulliken
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  print *, 'Hello world'
  read_wf = .True.
  touch read_wf
  call routine
! call print_non_ad_coupl
end


subroutine routine
 implicit none
 integer :: i,j
 double precision, allocatable :: mull_dens_atoms(:), density_mat(:,:), pop_mull_groups(:), density_mat_mo(:,:)
 double precision, allocatable :: pop_group_1(:)
 integer, allocatable :: iorder(:)
 allocate(mull_dens_atoms(nucl_num),density_mat(ao_num, ao_num),density_mat_mo(mo_num, mo_num),pop_mull_groups(n_groups))
 allocate(pop_group_1(n_states),iorder(n_states))
 print*,''
 print*,''
 print*,''
 do i = 1, n_states
 print*,''
  density_mat_mo(:,:) = ( one_e_dm_mo_alpha(:,:,i) + one_e_dm_mo_beta(:,:,i) )
  call mo_to_ao_no_overlap(density_mat_mo,mo_num,density_mat,ao_num)
  call density_mulliken_density_mat(density_mat, mull_dens_atoms)
  call get_pop_mull_groups(mull_dens_atoms, pop_mull_groups)

  print*,'state ',i
  print*,'Mulliken '
  do j = 1, nucl_num
   print*,nucl_charge(j),nucl_charge(j) - mull_dens_atoms(j),mull_dens_atoms(j)
  enddo
  print*,'Group mulliken '
  write(*,'(100(F16.10,X))')pop_mull_groups(:)
 print*,''
  pop_group_1(i) = -pop_mull_groups(1)
  iorder(i) = i
 enddo

!print*,'Sorting '
!call dsort(pop_group_1,iorder,n_states)
!do i = 1, n_states
! print*,i,iorder(i),-pop_group_1(i)
!enddo
end


subroutine get_pop_mull_groups(mull_dens_atoms, pop_mull_groups)
 implicit none
 double precision, intent(in) :: mull_dens_atoms(nucl_num)
 double precision, intent(out):: pop_mull_groups(n_groups)
 integer :: i,j,k
 pop_mull_groups = 0.d0
 do i = 1, n_groups
  do j = 1, n_atoms_per_group(i)
   k = list_atoms(j,i)
   pop_mull_groups(i) += nucl_charge(k) - mull_dens_atoms(k)
  enddo
 enddo
end


subroutine print_non_ad_coupl
 implicit none
 integer :: i,j
 print*,'non adiabatic coupling'
 do i = 1, n_states
  do j  = i+1, n_states
   print*,i,j,non_ad_coupling(j,i)
  enddo
 enddo

end
