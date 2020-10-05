program pouet
 implicit none
 read_wf = .True. 
 touch read_wf
 call print_local_energy

end

 BEGIN_PROVIDER [double precision, kin_e_array, (N_states,ntheta_psi_ex,nr1_psi_ex)]
&BEGIN_PROVIDER [double precision, non_hermit_e_array, (N_states,ntheta_psi_ex,nr1_psi_ex)]
&BEGIN_PROVIDER [double precision, pot_ee_array, (ntheta_psi_ex,nr1_psi_ex)]
&BEGIN_PROVIDER [double precision, pot_en_array, (ntheta_psi_ex,nr1_psi_ex)]
&BEGIN_PROVIDER [double precision, loc_e_array, (N_states,ntheta_psi_ex,nr1_psi_ex)]
&BEGIN_PROVIDER [double precision, psi_in_r1_r2_array, (N_states,ntheta_psi_ex,nr1_psi_ex)]
 implicit none
 integer :: i,j
 double precision :: kin_e(n_states),pot_ee,pot_en,non_hermit_e(n_states),e_loc(n_states),psi(n_states),loc_e(n_states),mu_in
 mu_in = mu_erf
 integer :: psi_occ(2, N_det)
 call give_occ_two_e_psi(psi_det,N_det,psi_occ)
 do j = 1, nr1_psi_ex
  do i = 1, ntheta_psi_ex
   call local_energy_htilde(r1_psi_ex(1,j),r2_psi_ex(1,i,j),mu_in,n_states,psi_occ,reigvec_trans(1,1),N_det,loc_e,kin_e,pot_ee,pot_en,non_hermit_e,psi)
   kin_e_array(:,i,j) = kin_e(:)
   non_hermit_e_array(:,i,j) = non_hermit_e(:)
   loc_e_array(:,i,j) = loc_e(:)
   psi_in_r1_r2_array(:,i,j) = psi(:)
   pot_ee_array(i,j) = pot_ee
   pot_en_array(i,j) = pot_en
  enddo
 enddo
END_PROVIDER 


subroutine print_local_energy
 implicit none
 integer :: i
 integer                        :: i_unit_output,getUnitAndOpen
 character*(128)                :: output
 PROVIDE ezfio_filename

BEGIN_TEMPLATE
 output=trim(ezfio_filename)//'.cusp_eloc_$X'
 print*,'output = ',trim(output)
 i_unit_output = getUnitAndOpen(output,'w')
 do i = 1, ntheta_psi_ex
  write(i_unit_output,'(100(F16.10,X))')theta_array_psi_ex(i,$X),r12_psi_ex(i,$X), & 
                                        loc_e_array(1,i,$X),psi_in_r1_r2_array(1,i,$X),& 
                                        kin_e_array(1,i,$X), kin_e_array(1,i,$X)/psi_in_r1_r2_array(1,i,$X),& 
                                        non_hermit_e_array(:,i,$X),non_hermit_e_array(:,i,$X)/psi_in_r1_r2_array(1,i,$X),& 
                                        pot_ee_array(i,$X),pot_en_array(i,$X)
 enddo
SUBST [ X]
 1;; 
 2;;
 3;;
 4;;
 5;;
END_TEMPLATE


end
