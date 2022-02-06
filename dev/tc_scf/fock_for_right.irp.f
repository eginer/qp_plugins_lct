
 BEGIN_PROVIDER [ double precision, good_hermit_tc_fock_mat, (mo_num, mo_num)]
  implicit none
  integer :: i,j
  good_hermit_tc_fock_mat = Fock_matrix_tc_mo_tot
 do j = 1, mo_num
  do i = 1, j-1
   good_hermit_tc_fock_mat(i,j) = Fock_matrix_tc_mo_tot(j,i) 
  enddo
 enddo
 END_PROVIDER 

 BEGIN_PROVIDER [ double precision, grad_good_hermit_tc_fock_mat]
 implicit none
 integer :: i,j
 grad_good_hermit_tc_fock_mat = 0.d0
 do i = 1, elec_alpha_num
  do j = elec_alpha_num+1, mo_num
   grad_good_hermit_tc_fock_mat += dabs(good_hermit_tc_fock_mat(i,j))
  enddo
 enddo
 END_PROVIDER 

 subroutine save_good_hermit_tc_eigvectors
  implicit none
  integer             :: sign
  character*(64)     :: label
  logical            :: output
  sign = 1
  label = "Canonical"
  output = .True.
  
  call mo_as_eigvectors_of_mo_matrix(good_hermit_tc_fock_mat,mo_num,mo_num,label,sign,output)                                                    
 end

 BEGIN_PROVIDER [ double precision, TC_right_HF_energy]
&BEGIN_PROVIDER [ double precision, TC_right_HF_one_electron_energy]
&BEGIN_PROVIDER [ double precision, TC_right_HF_two_e_hermit_energy]
&BEGIN_PROVIDER [ double precision, TC_right_HF_two_e_n_hermit_energy]
 implicit none
 BEGIN_DOC
 ! Hartree-Fock energy containing the nuclear repulsion, and its one- and two-body components.
 END_DOC
 integer :: i,j
 TC_right_HF_energy = nuclear_repulsion
 TC_right_HF_one_electron_energy = 0.d0
 TC_right_HF_two_e_hermit_energy = 0.d0
 TC_right_HF_two_e_n_hermit_energy = 0.d0
 do j=1,ao_num
   do i=1,ao_num
    TC_right_HF_two_e_hermit_energy += 0.5d0 * ( two_e_tc_hermit_integral_alpha(i,j) * SCF_density_matrix_ao_alpha(i,j) &
                                       +two_e_tc_hermit_integral_beta(i,j)  * SCF_density_matrix_ao_beta(i,j) )
    TC_right_HF_two_e_n_hermit_energy += 0.5d0 * ( two_e_tc_non_hermit_integral_alpha(i,j) * SCF_density_matrix_ao_alpha(i,j) &
                                       +two_e_tc_non_hermit_integral_beta(i,j)  * SCF_density_matrix_ao_beta(i,j) )
    TC_right_HF_one_electron_energy += ao_one_e_integrals(i,j) * (SCF_density_matrix_ao_alpha(i,j) + SCF_density_matrix_ao_beta (i,j) )
   enddo
 enddo
 TC_right_HF_energy += TC_right_HF_one_electron_energy + TC_right_HF_two_e_hermit_energy + TC_right_HF_two_e_n_hermit_energy 
END_PROVIDER
