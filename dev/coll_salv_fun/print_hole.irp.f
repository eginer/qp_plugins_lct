program print_hole
 implicit none
 BEGIN_DOC
 ! This programs writes the effective RS-DFT Hamiltonian into the EZFIO folder. 
 ! The next programs that will run unto the EZFIO folder will, by default, 
 !
 ! have the one- and two-body integrals loaded from the EZFIO data. 
 END_DOC
 read_wf = .true.
 touch read_wf
 !! total one-e integrals 
 io_mo_one_e_integrals = "None"
 touch io_mo_one_e_integrals  
 !! Vne integrals on the MO basis 
 io_mo_integrals_n_e = "None"
 touch io_mo_integrals_n_e
 !! kinetic integrals on the MO basis 
 io_mo_integrals_kinetic = "None"
 touch io_mo_integrals_kinetic 
 !! Vne integrals on the AO basis 
 io_ao_integrals_n_e = "None"
 touch io_ao_integrals_n_e 
 !! kinetic integrals on the AO basis 
 io_ao_integrals_kinetic = "None"
 touch io_ao_integrals_kinetic 

 !! regular 1/r12 integrals  on the MO basis
 io_mo_two_e_integrals = "None"
 touch io_mo_two_e_integrals
 !! regular 1/r12 integrals  on the AO basis
 io_ao_two_e_integrals = "None"
 touch io_ao_two_e_integrals
 !! integral of the effective potential 
 call routine_print
 call routine_print_2
end
!

subroutine routine_print
 implicit none
 integer :: i,ntheta,istate,n_taylor
 include 'constants.include.F'
 double precision :: r1(3),r2(3),mu_of_r,f_psi,rho2,n2_psi,psi,mu,full_jastrow_mu,thetamax
 double precision :: r,dtheta,theta,n2_hf,r12,jastrow,slater_ten_no,a0,psi_ex,f_mu
 double precision, allocatable :: mos_array_r1(:) , mos_array_r2(:), g0_rsdft,g0_jastrow,g0_exact
 double precision :: psi0,psi0_j,norm_n2,int_ovlp_n2_jaswtrow2,dm_a,dm_b,dens
 double precision :: psi0_j_new,g0_new,f_mu_new
 allocate(mos_array_r1(mo_num) , mos_array_r2(mo_num) )

 istate = 1
 n_taylor = 4
 ntheta = 10000
 thetamax = 2.d0 * pi
 dtheta = thetamax / dble(ntheta)
 r = 0.5d0 

 r1 = 0.d0
 r1(1) = r

 r2 = r1
 call get_two_e_psi_at_r1r2(r1,r2,psi0)
 psi_ex = 0.20400134d0
 g0_exact = psi_ex/dabs(psi0)
 call give_mu_of_r_cas(r,istate,mu_of_r,f_psi,n2_psi)
 mu = mu_of_r
 print*,''
 print*,''
 print*,''

 call dm_dft_alpha_beta_at_r(r1,dm_a,dm_b)
 dens = 0.5d0 * (dm_a + dm_b)
 norm_n2 = int_ovlp_n2_jaswtrow2(r1,mu_of_r,istate,n_taylor)
 norm_n2 = dens / norm_n2
 print*,'norm_n2 = ',norm_n2

! mu = -1.d0/(2.d0*dsqrt(pi)*dlog(g0_exact))
 a0 = 2.d0 * dsqrt(pi) * mu
 g0_jastrow = full_jastrow_mu(mu,0.d0)
 psi0_j = dabs(psi0 * g0_jastrow)
 g0_new     = dexp(f_mu_new(mu,0.d0))
 psi0_j_new = dabs(psi0 * g0_new)
 g0_rsdft   = 1.d0/(1.d0+ 1.d0/a0)
 print*,'mu = ',mu
 print*,'g0 jastrow / Psi   = ',g0_jastrow, g0_jastrow*dabs(psi0)
 print*,'g0 rsdft   / Psi   = ',g0_rsdft  , g0_rsdft  *dabs(psi0)
 print*,'g0 exact   / Psi   = ',g0_exact  , psi_ex 

 theta = -thetamax * 0.5d0
 do i = 1, ntheta
  r2(1) = r * dcos(theta)
  r2(2) = r * dsin(theta)
  call get_two_e_psi_at_r1r2(r1,r2,psi)
  r12 = dsqrt( (r1(1) - r2(1))**2.d0 +  (r1(2) - r2(2))**2.d0 + (r1(3) - r2(2))**2.d0 )
  write(34,'(100(F16.10,X))')theta, r12, dabs(psi), dabs(psi)*full_jastrow_mu(mu,r12),dabs(psi)*dexp(f_mu_new(mu,r12))
!  , dsqrt(norm_n2) * dabs(psi)*full_jastrow_mu(mu,r12), (1.d0 - derf(mu*r12))/r12
!  write(35,'(100(F16.10,X))')r12,dabs(psi),dabs(psi)*full_jastrow_mu(mu,r12),dabs(psi)/psi0_j*full_jastrow_mu(mu,r12),full_jastrow_mu(mu,r12),f_mu(mu,r12)
  theta += dtheta
!   , n2_hf ,n2_psi 
 enddo
end

subroutine routine_print_2
 implicit none
 integer :: i,ntheta,istate
 include 'constants.include.F'
 double precision, allocatable :: theta_array(:),psi_ex_array(:)
 integer :: n_pts_ex
 double precision :: r1(3),r2(3),mu_of_r,f_psi,rho2,n2_psi,psi,mu,full_jastrow_mu,thetamax
 double precision :: r,dtheta,theta,n2_hf,r12,jastrow,slater_ten_no,a0,psi0_ex,f_mu
 double precision, allocatable :: mos_array_r1(:) , mos_array_r2(:), g0_rsdft,g0_jastrow,g0_exact
 double precision :: psi0,psi0_j
 double precision :: psi0_j_new,g0_new,f_mu_new

 istate = 1
 ntheta = 21
 allocate(theta_array(ntheta),psi_ex_array(ntheta))
 theta_array(1)  = -3.14159265                                                                                                   
 theta_array(2)  = -2.82743339   
 theta_array(3)  = -2.51327412   
 theta_array(4)  = -2.19911486   
 theta_array(5)  = -1.88495559   
 theta_array(6)  = -1.57079633   
 theta_array(7)  = -1.25663706   
 theta_array(8 ) = -0.9424778    
 theta_array(9 ) = -0.62831853   
 theta_array(10) = -0.31415927   
 theta_array(11) = 0.            
 theta_array(12) = 0.31415927    
 theta_array(13) = 0.62831853    
 theta_array(14) = 0.9424778     
 theta_array(15) = 1.25663706    
 theta_array(16) = 1.57079633    
 theta_array(17) = 1.88495559    
 theta_array(18) = 2.19911486    
 theta_array(19) = 2.51327412    
 theta_array(20) = 2.82743339    
 theta_array(21) = 3.14159265    


 psi_ex_array(1)  = 0.27979608
 psi_ex_array(2)  = 0.279069
 psi_ex_array(3)  = 0.27688511
 psi_ex_array(4)  = 0.27323624
 psi_ex_array(5)  = 0.26810778
 psi_ex_array(6)  = 0.26147715
 psi_ex_array(7)  = 0.25331132
 psi_ex_array(8 ) = 0.24356313
 psi_ex_array(9 ) = 0.23216563
 psi_ex_array(10) = 0.21902378
 psi_ex_array(11) = 0.20400134
 psi_ex_array(12) = 0.21902378
 psi_ex_array(13) = 0.23216563
 psi_ex_array(14) = 0.24356313
 psi_ex_array(15) = 0.25331132
 psi_ex_array(16) = 0.26147715
 psi_ex_array(17) = 0.26810778
 psi_ex_array(18) = 0.27323624
 psi_ex_array(19) = 0.27688511
 psi_ex_array(20) = 0.279069
 psi_ex_array(21) = 0.27979608


 r = 0.5d0  ! distance to the nucleus
 r1 = 0.d0
 r1(1) = r
 r2 = r1
 call get_two_e_psi_at_r1r2(r1,r2,psi0)
 psi0_ex = 0.20400134d0
 g0_exact = psi0_ex/dabs(psi0)
 call give_mu_of_r_cas(r,istate,mu_of_r,f_psi,n2_psi)
 print*,''
 print*,''
 print*,''
 mu = mu_of_r
! mu = -1.d0/(2.d0*dsqrt(pi)*dlog(g0_exact))
 a0 = 2.d0 * dsqrt(pi) * mu
 g0_jastrow = full_jastrow_mu(mu,0.d0)
 psi0_j = dabs(psi0 * g0_jastrow)
 g0_new     = dexp(f_mu_new(mu,0.d0))
 psi0_j_new = dabs(psi0 * g0_new)
 g0_rsdft   = 1.d0/(1.d0+ 1.d0/a0)
 print*,'mu = ',mu
 print*,'g0 jastrow / Psi   = ',g0_jastrow, g0_jastrow*dabs(psi0)
 print*,'g0 new     / Psi   = ',g0_new    , g0_new    *dabs(psi0)
 print*,'g0 rsdft   / Psi   = ',g0_rsdft  , g0_rsdft  *dabs(psi0)
 print*,'g0 exact   / Psi   = ',g0_exact  , psi0_ex 

 do i = 1, ntheta
  theta = theta_array(i)
  r2(1) = r * dcos(theta)
  r2(2) = r * dsin(theta)
  call get_two_e_psi_at_r1r2(r1,r2,psi)
  r12 = dsqrt( (r1(1) - r2(1))**2.d0 +  (r1(2) - r2(2))**2.d0 + (r1(3) - r2(2))**2.d0 )
  write(32,'(100(F16.10,X))')theta,r12, psi_ex_array(i),psi_ex_array(i)/psi0_ex
  write(33,'(100(F16.10,X))')r12, dabs(psi)*full_jastrow_mu(mu,r12)/psi0_j,psi_ex_array(i)/psi0_ex,dabs(psi)*dexp(f_mu_new(mu,r12))/psi0_j_new
 enddo
end
