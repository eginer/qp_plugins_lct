BEGIN_PROVIDER [integer, ntheta_psi_ex]
 implicit none 
 ntheta_psi_ex = 100
END_PROVIDER 

BEGIN_PROVIDER [integer, nr1_psi_ex]
 implicit none
 nr1_psi_ex = 7
END_PROVIDER 

BEGIN_PROVIDER [double precision, r1_psi_ex, (3,nr1_psi_ex)]
 implicit none
 r1_psi_ex = 0.d0
 r1_psi_ex(1,1) = 0.1d0
 r1_psi_ex(1,2) = 0.25d0
 r1_psi_ex(1,3) = 0.50d0
 r1_psi_ex(1,4) = 0.75d0
 r1_psi_ex(1,5) = 1.00d0
 r1_psi_ex(1,6) = 1.50d0
 r1_psi_ex(1,7) = 2.00d0
END_PROVIDER 


 BEGIN_PROVIDER [double precision, r2_psi_ex, (3,ntheta_psi_ex,nr1_psi_ex)]
&BEGIN_PROVIDER [double precision, r12_psi_ex, (ntheta_psi_ex,nr1_psi_ex)]
 implicit none
 integer :: i,j
 double precision :: theta,r


 do j = 1, nr1_psi_ex
  do i = 1, ntheta_psi_ex
   r2_psi_ex(3,i,j) = r1_psi_ex(3,j)
   theta = theta_array_psi_ex(i,j)
   r = r1_psi_ex(1,j) 
   r2_psi_ex(1,i,j) = r * dcos(theta)
   r2_psi_ex(2,i,j) = r * dsin(theta)
   r12_psi_ex(i,j) = dsqrt( (r1_psi_ex(1,j) - r2_psi_ex(1,i,j))**2.d0 +  (r1_psi_ex(2,j) - r2_psi_ex(2,i,j))**2.d0 + (r1_psi_ex(3,j) - r2_psi_ex(3,i,j))**2.d0 )
  enddo
 enddo

 
END_PROVIDER 

BEGIN_PROVIDER [double precision, psi_fci_array, (ntheta_psi_ex,nr1_psi_ex)]
 implicit none
 integer :: i , j
!call diagonalize_CI
 double precision :: psi
 do j = 1, nr1_psi_ex
  do i = 1, ntheta_psi_ex
   call get_two_e_psi_at_r1r2(r1_psi_ex(1,j),r2_psi_ex(1,i,j),psi)
   psi_fci_array(i,j) = dabs(psi)
  enddo
 enddo
END_PROVIDER 

BEGIN_PROVIDER [double precision, psi_hf_array, (ntheta_psi_ex,nr1_psi_ex)]
 implicit none
 integer :: i , j
!call diagonalize_CI
 double precision :: mos_array_r1(mo_num), mos_array_r2(mo_num)
 double precision :: r1(3), r2(3)
 do j = 1, nr1_psi_ex
  do i = 1, ntheta_psi_ex
   call give_all_mos_at_r(r1,mos_array_r1)
   call give_all_mos_at_r(r2,mos_array_r2)
   r1(:) = r1_psi_ex(:,j)
   r2(:) = r2_psi_ex(:,i,j)
   psi_hf_array(i,j) = dabs(mos_array_r1(1)*mos_array_r2(1))
  enddo
 enddo
END_PROVIDER 


 BEGIN_PROVIDER [double precision, psi_right_array, (ntheta_psi_ex,nr1_psi_ex)]
&BEGIN_PROVIDER [double precision, psi_right_jastrow_array, (ntheta_psi_ex,nr1_psi_ex)]
&BEGIN_PROVIDER [double precision, psi_right_jastrow_norm_array, (ntheta_psi_ex,nr1_psi_ex)]
 implicit none
 integer :: i,j
 do i = 1, N_det
  psi_coef(i,1) = reigvec_trans(i,1)/dsqrt(reigvec_trans_norm(1))
 enddo
 touch psi_coef
 double precision :: psi,inv_norm,mu,r12,full_jastrow_mu
 mu = mu_erf 
 print*,'<Phi | e^{2 j(r12,mu)} | Phi> = ',norm_n2_jastrow_cst_mu(1)
 inv_norm = 1.d0/dsqrt(norm_n2_jastrow_cst_mu(1))
 do j = 1, nr1_psi_ex
  do i = 1, ntheta_psi_ex
   r12 = r12_psi_ex(i,j)
   call get_two_e_psi_at_r1r2(r1_psi_ex(1,j),r2_psi_ex(1,i,j),psi)
   psi_right_array(i,j) = dabs(psi)
   psi_right_jastrow_array(i,j) = dabs(psi) * full_jastrow_mu(mu,r12)
   psi_right_jastrow_norm_array(i,j) = dabs(psi) * full_jastrow_mu(mu,r12) * inv_norm
  enddo
 enddo
END_PROVIDER 

subroutine print_psi_exc_psi_trans
 implicit none
 integer :: i
 integer                        :: i_unit_output,getUnitAndOpen
 character*(128)                :: output
 PROVIDE ezfio_filename

BEGIN_TEMPLATE
 output=trim(ezfio_filename)//'.cusp_wf_$X'
 print*,'output = ',trim(output)
 i_unit_output = getUnitAndOpen(output,'w')
 do i = 1, ntheta_psi_ex
  write(i_unit_output,'(100(F16.10,X))')theta_array_psi_ex(i,$X),r12_psi_ex(i,$X),array_psi_ex(i,$X), & 
                                        psi_fci_array(i,$X), &
                                        psi_right_array(i,$X), &
                                        psi_right_jastrow_array(i,$X),psi_right_jastrow_norm_array(i,$X)
 enddo
SUBST [ X]
 1;; 
 2;;
 3;;
 4;;
 5;;
END_TEMPLATE


end

subroutine print_psi_hf_psi_trans
 implicit none
 integer :: i
 integer                        :: i_unit_output,getUnitAndOpen
 character*(128)                :: output
 PROVIDE ezfio_filename

BEGIN_TEMPLATE
 output=trim(ezfio_filename)//'.cusp_wf_$X'
 print*,'output = ',trim(output)
 i_unit_output = getUnitAndOpen(output,'w')
 do i = 1, ntheta_psi_ex
  write(i_unit_output,'(100(F16.10,X))')theta_array_psi_ex(i,$X),psi_hf_array(i,$X)
 enddo
SUBST [ X]
 1;; 
 2;;
 3;;
 4;;
 5;;
END_TEMPLATE


end
