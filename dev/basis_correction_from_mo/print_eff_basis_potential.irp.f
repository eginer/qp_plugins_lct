program print_eff_basis_potential
  implicit none
  BEGIN_DOC
! Compute <Psi | V^B |Psi> 
  END_DOC
  double precision, allocatable :: pot_basis_alpha_mo_from_file(:,:)
  double precision, allocatable :: pot_basis_beta_mo_from_file(:,:)
  double precision :: potential_tot_alpha, potential_tot_beta
  integer :: i, j, ii, jj, nst 

  allocate(pot_basis_alpha_mo_from_file(mo_num, mo_num))
  allocate(pot_basis_beta_mo_from_file(mo_num, mo_num))

  ! READ alpha potential in MO basis:
  open(12, file="pot_basis_alpha_mo.txt") 
  read(12,*) pot_basis_alpha_mo_from_file

  ! READ beta potential in MO basis:
  open(13, file="pot_basis_beta_mo.txt") 
  read(13,*) pot_basis_beta_mo_from_file

  provide data_one_e_dm_alpha_mo
  provide data_one_e_dm_beta_mo
  nst=1 ! state

  potential_tot_alpha = 0.d0
  potential_tot_beta = 0.d0

  do ii=1, n_act_orb
   i = list_act(ii)
   do jj=1, n_act_orb
    j = list_act(jj)
    potential_tot_alpha+= data_one_e_dm_alpha_mo(ii,jj,nst)*pot_basis_alpha_mo_from_file(i,j)
    potential_tot_beta+= data_one_e_dm_beta_mo(ii,jj,nst)*pot_basis_beta_mo_from_file(i,j)
   enddo
  enddo

  close(12)
  close(13)

 print*,'POTENTIAL ALPHA =', potential_tot_alpha
 print*,'POTENTIAL BETA =', potential_tot_beta
 print*,'POTENTIAL TOT =', potential_tot_alpha + potential_tot_beta

end
