BEGIN_PROVIDER [integer, ntheta_psi_ex]
 implicit none 
 ntheta_psi_ex = 21
END_PROVIDER 

 BEGIN_PROVIDER [double precision, theta_array_psi_ex, (ntheta_psi_ex)]
&BEGIN_PROVIDER [double precision, array_psi_ex, (ntheta_psi_ex)]
&BEGIN_PROVIDER [double precision, r1_psi_ex, (3)]
&BEGIN_PROVIDER [double precision, r2_psi_ex, (3,ntheta_psi_ex)]
&BEGIN_PROVIDER [double precision, r12_psi_ex, (ntheta_psi_ex)]
 implicit none
 integer :: i
 double precision :: theta,r
 r = 0.5d0 
 r1_psi_ex = 0.d0
 r1_psi_ex(1) = r

 theta_array_psi_ex(1)  = -3.14159265                                                                                                   
 theta_array_psi_ex(2)  = -2.82743339   
 theta_array_psi_ex(3)  = -2.51327412   
 theta_array_psi_ex(4)  = -2.19911486   
 theta_array_psi_ex(5)  = -1.88495559   
 theta_array_psi_ex(6)  = -1.57079633   
 theta_array_psi_ex(7)  = -1.25663706   
 theta_array_psi_ex(8 ) = -0.9424778    
 theta_array_psi_ex(9 ) = -0.62831853   
 theta_array_psi_ex(10) = -0.31415927   
 theta_array_psi_ex(11) = 0.            
 theta_array_psi_ex(12) = 0.31415927    
 theta_array_psi_ex(13) = 0.62831853    
 theta_array_psi_ex(14) = 0.9424778     
 theta_array_psi_ex(15) = 1.25663706    
 theta_array_psi_ex(16) = 1.57079633    
 theta_array_psi_ex(17) = 1.88495559    
 theta_array_psi_ex(18) = 2.19911486    
 theta_array_psi_ex(19) = 2.51327412    
 theta_array_psi_ex(20) = 2.82743339    
 theta_array_psi_ex(21) = 3.14159265    


 array_psi_ex(1)  = 0.27979608
 array_psi_ex(2)  = 0.279069
 array_psi_ex(3)  = 0.27688511
 array_psi_ex(4)  = 0.27323624
 array_psi_ex(5)  = 0.26810778
 array_psi_ex(6)  = 0.26147715
 array_psi_ex(7)  = 0.25331132
 array_psi_ex(8 ) = 0.24356313
 array_psi_ex(9 ) = 0.23216563
 array_psi_ex(10) = 0.21902378
 array_psi_ex(11) = 0.20400134
 array_psi_ex(12) = 0.21902378
 array_psi_ex(13) = 0.23216563
 array_psi_ex(14) = 0.24356313
 array_psi_ex(15) = 0.25331132
 array_psi_ex(16) = 0.26147715
 array_psi_ex(17) = 0.26810778
 array_psi_ex(18) = 0.27323624
 array_psi_ex(19) = 0.27688511
 array_psi_ex(20) = 0.279069
 array_psi_ex(21) = 0.27979608

 do i = 1, ntheta_psi_ex
  r2_psi_ex(3,i) = r1_psi_ex(3)
  theta = theta_array_psi_ex(i)
  r2_psi_ex(1,i) = r * dcos(theta)
  r2_psi_ex(2,i) = r * dsin(theta)
  r12_psi_ex(i) = dsqrt( (r1_psi_ex(1) - r2_psi_ex(1,i))**2.d0 +  (r1_psi_ex(2) - r2_psi_ex(2,i))**2.d0 + (r1_psi_ex(3) - r2_psi_ex(3,i))**2.d0 )
 enddo

 
END_PROVIDER 

BEGIN_PROVIDER [double precision, psi_fci_array, (ntheta_psi_ex)]
 implicit none
 integer :: i , j
 do i = 1, N_states
  do j = 1, N_det
   psi_coef(j,i) =  CI_eigenvectors(j,i)
  enddo
 enddo
 touch psi_coef
 double precision :: psi
 do i = 1, ntheta_psi_ex
  call get_two_e_psi_at_r1r2(r1_psi_ex,r2_psi_ex(1,i),psi)
  psi_fci_array(i) = dabs(psi)
 enddo
END_PROVIDER 

 BEGIN_PROVIDER [double precision, psi_right_array, (ntheta_psi_ex)]
&BEGIN_PROVIDER [double precision, psi_right_jastrow_array, (ntheta_psi_ex)]
&BEGIN_PROVIDER [double precision, psi_right_jastrow_norm_array, (ntheta_psi_ex)]
 implicit none
 integer :: i
 do i = 1, N_det
  psi_coef(i,1) = reigvec_trans(i,1)/dsqrt(reigvec_trans_norm(1))
 enddo
 touch psi_coef
 double precision :: psi,inv_norm,mu,r12,full_jastrow_mu
 mu = mu_erf 
 inv_norm = 1.d0/dsqrt(norm_n2_jastrow_cst_mu(1))
 do i = 1, ntheta_psi_ex
  r12 = r12_psi_ex(i)
  call get_two_e_psi_at_r1r2(r1_psi_ex,r2_psi_ex(1,i),psi)
  psi_right_array(i) = dabs(psi)
  psi_right_jastrow_array(i) = dabs(psi) * full_jastrow_mu(mu,r12)
  psi_right_jastrow_norm_array(i) = dabs(psi) * full_jastrow_mu(mu,r12) * inv_norm
 enddo
END_PROVIDER 

subroutine print_psi_exc_psi_trans
 implicit none
 integer :: i
 integer                        :: i_unit_output,getUnitAndOpen
 character*(128)                :: output
 PROVIDE ezfio_filename
 output=trim(ezfio_filename)//'.cusp_wf'
 print*,'output = ',trim(output)
 i_unit_output = getUnitAndOpen(output,'w')
 do i = 1, ntheta_psi_ex
  write(i_unit_output,'(100(F16.10,X))')theta_array_psi_ex(i),r12_psi_ex(i),array_psi_ex(i), &
                                        psi_fci_array(i),psi_right_array(i),                 &
                                        psi_right_jastrow_array(i),psi_right_jastrow_norm_array(i)
 enddo
end
