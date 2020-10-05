program pouet
 implicit none
 read_wf = .True. 
 touch read_wf
! call print_local_energy
 call routine
end

subroutine routine
 implicit none
 double precision :: r1(3), r2(3),x,dx,xmax,r12
 integer :: i,nx
 double precision :: kin_e(n_states),pot_ee,pot_en,non_hermit_e(n_states),e_loc(n_states),psi(n_states),loc_e(n_states),mu_in
 integer :: psi_occ(2, N_det)
 xmax = 5.d0
 nx = 100
 dx = xmax/dble(nx)
 r1 = 0.d0
 r2 = 0.d0
 call give_occ_two_e_psi(psi_det,N_det,psi_occ)
 x = -2.5d0
! do i = 1, nx
!  r1(1) = x
!  call local_energy_htilde(r1,r2,mu_in,n_states,psi_occ,psi_coef,N_det,loc_e,kin_e,pot_ee,pot_en,non_hermit_e,psi)
!  write(33,'(100(F16.10,X))') x,loc_e,kin_e,pot_ee,pot_en,non_hermit_e,psi
!  x += dx
! enddo
 integer :: j
 j = 3
 r1(:) = r1_psi_ex(:,j)
 mu_in = 1.d+9
 do i = 1, ntheta_psi_ex
  r2(:) = r2_psi_ex(:,i,j)
  r12 = dsqrt( (r1(1) - r2(1))**2.d0 + (r1(2) - r2(2))**2.d0 + (r1(3) - r2(3))**2.d0 )
  call local_energy_htilde(r1,r2,mu_in,n_states,psi_occ,psi_coef,N_det,loc_e,kin_e,pot_ee,pot_en,non_hermit_e,psi)
  !                                     1               2    3     4     5      6     7            8
  write(33,'(X,100(F16.10,X))')theta_array_psi_ex(i,j),r12,loc_e,kin_e/psi,pot_ee,pot_en,non_hermit_e/psi,psi
 enddo
 j = 3
 mu_in = mu_erf
 do i = 1, ntheta_psi_ex
  r2(:) = r2_psi_ex(:,i,j)
  r12 = dsqrt( (r1(1) - r2(1))**2.d0 + (r1(2) - r2(2))**2.d0 + (r1(3) - r2(3))**2.d0 )
  !                                     1               2    3     4     5      6     7            8
  call local_energy_htilde(r1,r2,mu_in,n_states,psi_occ,reigvec_trans(1,1),N_det,loc_e,kin_e,pot_ee,pot_en,non_hermit_e,psi)
  write(34,'(X,100(F16.10,X))')theta_array_psi_ex(i,j),r12,loc_e,kin_e,pot_ee,pot_en,non_hermit_e,psi
 enddo
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
