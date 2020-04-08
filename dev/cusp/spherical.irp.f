 subroutine test_spherical_grid
 implicit none
 include 'utils/constants.include.F'
 integer :: i,j,k,istate,nr12
 double precision :: r12max
 double precision :: r1(3),r2(3)
 double precision :: r12,dr12
 double precision :: two_dm_hf,two_dm_hf_laplacian,total_hf_dm,two_rdm
 double precision :: two_dm_psi,two_dm_psi_laplacian,total_psi_dm
 double precision :: dm_a_2,dm_b_2,dm_a_1,dm_b_1
 double precision :: mos_array2(mo_num),mos_array1(mo_num)
 double precision :: n2_hf,accu_hf,accu_psi,accu_j

 istate = 1
 print*,'\\\\\\\\\\\\\\\\\'
 print*,' '
 r1(1) = 0.0d0
 r1(2) = 0.0d0
 r1(3) = r1_z_cups
 call f_HF_valence_ab(r1,r1,f_HF_val_ab,two_bod_dens) 
 r12max = 0.01d0
 nr12 = 500
 dr12 = r12max/nr12
 r12=0.00001d0
 call give_all_mos_at_r(r1,mos_array1)
 double precision :: corr_hole, on_top_00
 write(33,*)'#r12, n2 FCI        ,   n2 HF           ,  n2_c               , J'
 corr_hole = 0.d0
 call two_body_r1r2(r1,r1,istate,on_top_00)
 double precision :: g0_physical,rho_a,rho_b
 call dm_dft_alpha_beta_at_r(r1,rho_a,rho_b)
! g0_physical = 2.d0 * on_top_00/(2.d0 * mos_array1(1)**2)**2
 g0_physical = 2.d0 * on_top_00/(rho_a+rho_b)**2
 double precision :: dmu,g0_ueg,g0_UEG_mu,mu
 dmu = 0.01d0
 mu = 0.000001d0
 print*,'rho ',rho_b + rho_a
 do i = 1, 10000
  g0_ueg = g0_UEG_mu(mu,rho_a,rho_b)
!  if(dabs(g0_ueg - g0_physical).lt.0.01d0)then
!   print*,g0_ueg , g0_physical, mu
!  endif
  mu += dmu
  write(34,'(100(F16.10,X))')mu,g0_ueg,g0_physical
 enddo
 double precision :: f_HF_val_ab,two_bod_dens,w_hf,mu_hf
 call f_HF_valence_ab(r1,r1,f_HF_val_ab,two_bod_dens) 
 w_hf  = f_HF_val_ab/two_bod_dens
 mu_hf =  w_hf * dsqrt(dacos(-1.d0)) * 0.5d0
 print*,'mu_hf               ',mu_hf
 stop
 do i = 1, nr12
  print*,'r12 = ',r12
  ! loop for a given r12 
  accu_psi = 0.d0
  accu_hf  = 0.d0
  accu_j   = 0.d0
  do k = 1, n_points_integration_angular
   r2(1)= r1(1)+ r12* angular_quadrature_points(k,1)
   r2(2)= r1(2)+ r12* angular_quadrature_points(k,2)
   r2(3)= r1(3)+ r12* angular_quadrature_points(k,3)

   call give_all_mos_at_r(r2,mos_array2)
   n2_hf = mos_array2(1)**2 * mos_array1(1)**2

   call two_body_r1r2(r1,r2,istate,two_rdm)

   accu_psi += two_rdm * weights_angular_points(k)
   accu_hf  += n2_hf   * weights_angular_points(k)
   accu_j   += two_rdm / n2_hf * weights_angular_points(k)
  enddo
  accu_psi = accu_psi / (4.d0 * pi)
  accu_hf  = accu_hf  / (4.d0 * pi)
  accu_j   = accu_j   / (4.d0 * pi)
  write(33,'(10(F16.10,X))')r12,accu_psi,accu_hf,accu_psi - accu_hf, accu_psi/accu_hf,accu_j,accu_psi/on_top_00
  r12 += dr12
  corr_hole += (accu_psi - accu_hf) * r12**2 *dr12 * 4.d0 * pi
 enddo
 print*,'corr_hole = ',corr_hole
 call spherical_averaged_two_dm_HF_at_second_order(r1,r12,two_dm_hf,two_dm_hf_laplacian,total_hf_dm)
 print*,'total_hf_dm = ',total_hf_dm
 call spherical_averaged_two_dm_at_second_order(r1,r12,istate,two_dm_psi,two_dm_psi_laplacian,total_psi_dm)
 print*,'total_psi_dm= ',total_psi_dm
 print*,'two_dm_hf_laplacian ',two_dm_hf_laplacian
 print*,'two_dm_psi_laplacian',two_dm_psi_laplacian
 print*,'two_dm_psi          ',two_dm_psi
 print*,'two_dm_hf           ',two_dm_hf
 print*,'J = n2_fci/n2_hf    ',two_dm_psi/two_dm_hf
 print*,'Delta J / J         ',two_dm_psi_laplacian/two_dm_psi - two_dm_hf_laplacian/two_dm_hf
 double precision :: corr_lapl,mu_new
 corr_lapl =  0.5d0 * (two_dm_psi_laplacian - two_dm_hf_laplacian)/(two_dm_psi - two_dm_hf)
 print*,'corr_lapl           ',corr_lapl

 mu_new = 3.d0*dsqrt(dacos(-1.d0))/2.d0  *   corr_lapl
 print*,'mu_new              ',mu_new



end

