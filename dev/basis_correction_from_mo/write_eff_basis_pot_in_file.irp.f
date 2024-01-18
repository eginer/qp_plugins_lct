program write_eff_basis_pot_in_file
  implicit none
  BEGIN_DOC
! READ the one-e density matrix (in MO basis )in a file and write it in aux_quantities/
! also convert from MO to AO basis
! Input files are labeled mo_alpha_1rdm.txt and mo_beta_1rdm.txt
  END_DOC

  !read_wf=.True.
  !touch read_wf 
  call build_dm_for_aux_quantities

  end 

  subroutine build_dm_for_aux_quantities
  implicit none
  ! Faire plutôt des allocates
  integer :: i, ii, nst
 
  provide pot_basis_alpha_mo
  provide pot_basis_beta_mo

  nst=1

  ! File for alpha potential in MO basis:
  open(12, file="pot_basis_alpha_mo.txt") 

  ! File for alpha potential in MO basis:
  open(13, file="pot_basis_beta_mo.txt") 

  do i=1, mo_num
   write (12, "(*(f14.10))") pot_basis_alpha_mo(i,:,nst) 
   write (13, "(*(f14.10))") pot_basis_beta_mo(i,:,nst) 
  enddo
   
  close(12)
  close(13)
end
