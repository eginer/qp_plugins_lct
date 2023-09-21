BEGIN_PROVIDER [integer, n_groups]
 implicit none
 n_groups = 2 
END_PROVIDER 


BEGIN_PROVIDER [integer, n_atoms_per_group, (n_groups)]
 implicit none
 n_atoms_per_group(1) = 2
 n_atoms_per_group(2) = 4
! n_atoms_per_group(1) = 1
! n_atoms_per_group(2) = 1
END_PROVIDER 

BEGIN_PROVIDER [integer, list_atoms, (nucl_num, n_groups)]
 implicit none
 integer :: i,j

!  list_atoms(1,1) = 1
!  list_atoms(1,2) = 2

 list_atoms(1,1) = 2
 list_atoms(2,1) = 3

 list_atoms(1,2) = 1
 list_atoms(2,2) = 4
 list_atoms(3,2) = 5
 list_atoms(4,2) = 6

 print*,''
 print*,'The groups are ...'
 print*,''
 do i = 1, n_groups
  print*,'Group ',i
  do j = 1, n_atoms_per_group(i)
   print*,list_atoms(j,i),nucl_charge(list_atoms(j,i))
  enddo
 enddo

END_PROVIDER 


