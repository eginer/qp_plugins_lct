
subroutine test_int_f_special
 implicit none
 BEGIN_DOC
! test the function ovlp_f_tilde_phi_ij
 END_DOC
 include 'utils/constants.include.F'
 double precision :: A_center(3), B_center(3), alpha,beta      ! basis functions
 integer          :: power_A(3), power_B(3)                    ! basis functions 
 double precision :: r1(3),mu

 double precision :: ovlp_exp_f_phi_ij
 double precision :: numerical,  analytical, weight,primitive_value_explicit
 double precision :: numerical_c,analytical_c,r_12,r2(3)
 double precision :: gauss_a, gauss_b, gauss_d, coulomb, jastrow
 integer          :: ipoint,ao_i,ao_j,num_A,num_B,i,j,k,l
 double precision :: abs_error_av, relat_error_av,icount
 double precision :: abs_error_fit_av, relat_error_fit_av
 double precision :: abs_error_c_av, relat_error_c_av,slater_fit_gam
 double precision :: exp_dl,f_tilde_fit,full_jastrow,f_mu,full_jastrow_mu,tilde_f_mu
 double precision :: numerical_fit,shank_general,jastrow_fit,full_jastrow_fit
 double precision :: de_abs,de_relat,erf_exp_f_phi_ij,ao_ints(ao_num,ao_num),aos_array(ao_num),ao_erf_ints(ao_num,ao_num)
 double precision :: mo_ints(mo_num,mo_num),mos_array(mo_num),mo_erf_ints(mo_num,mo_num)
 double precision, allocatable :: array(:)
 integer :: n_taylor
 double precision :: a0,exponent_exp
 exponent_exp = -1.d0
 n_taylor = 8
 allocate(array(0:n_taylor))
 mu = 1.5d0
 a0 = 1.d0/(2.d0*dsqrt(pi)*mu)
!call test_fit(mu)
 r1 = 0.1d0
! r1(1) = 0.1d0

 abs_error_av   = 0.d0
 relat_error_av = 0.d0
 abs_error_c_av   = 0.d0
 relat_error_c_av = 0.d0
 icount = 0.d0
! call give_jastrow2_ovlp_ints_ao(mu,r1,n_taylor,ao_ints)
! call give_jastrow2_erf_ints_ao(mu,r1,n_taylor,ao_erf_ints)
 call give_jastrow2_ovlp_ints_mo(mu,r1,n_taylor,mo_ints,exponent_exp)
 call give_jastrow2_erf_ints_mo(mu,mu,r1,n_taylor,mo_erf_ints)
 do i = 1, mo_num
  do j = 1, mo_num
   print*,'i,j,mu',i,j,mu
!    num_A = ao_nucl(i)
!    power_A(1:3)= ao_power(i,1:3)
!    A_center(1:3) = nucl_coord(num_A,1:3)
!   
!    num_B = ao_nucl(j)
!    power_B(1:3)= ao_power(j,1:3)
!    B_center(1:3) = nucl_coord(num_B,1:3)
!    do k = 1, ao_prim_num(i)
!     alpha = ao_expo_ordered_transp(k,i)     
!     do l = 1, ao_prim_num(j)
!      beta = ao_expo_ordered_transp(l,j)     
       ! analytical integral 
!       analytical   = ovlp_exp_f_phi_ij(mu,r1,A_center,B_center,power_A,power_B,alpha,beta,n_taylor,exponent_exp)
!       analytical_c =  erf_exp_f_phi_ij(mu,mu,r1,A_center,B_center,power_A,power_B,alpha,beta,n_taylor)
!        analytical   = ao_ints(j,i)
!        analytical_c = ao_erf_ints(j,i)
        analytical   = mo_ints(j,i)
        analytical_c = mo_erf_ints(j,i)

       ! numerical  integral 
       numerical   = 0.d0
       numerical_c = 0.d0
       numerical_fit = 0.d0
       icount += 1.d0
       do ipoint = 1, n_points_final_grid
        r2(1) = final_grid_points(1,ipoint)
        r2(2) = final_grid_points(2,ipoint)
        r2(3) = final_grid_points(3,ipoint)
!        call give_all_aos_at_r(r2,aos_array)
        call give_all_mos_at_r(r2,mos_array)
        weight = final_weight_at_r_vector(ipoint)

!        gauss_a = primitive_value_explicit(power_A,A_center,alpha,r2)
!        gauss_b = primitive_value_explicit(power_B,B_center,beta ,r2)
!        gauss_a = aos_array(i)
!        gauss_b = aos_array(j)

        gauss_a = mos_array(i)
        gauss_b = mos_array(j)
        r_12     = dsqrt( (r1(1) - r2(1))**2.d0 + (r1(2) - r2(2))**2.d0 + (r1(3) - r2(3))**2.d0 )
        if(dabs(r_12).lt.1.d-6)then
         coulomb = 2.d0 * mu / sqpi - 2.d0 * mu**3 * r_12**2 / (3.d0 *sqpi) 
        else
         coulomb = derf(mu * r_12)/r_12
        endif
        jastrow_fit = f_tilde_fit(mu,r_12)
        full_jastrow_fit = dexp(exponent_exp * jastrow_fit*a0) 
        jastrow = tilde_f_mu(mu,r_12)
        full_jastrow = dexp(exponent_exp * jastrow*a0)
        numerical_c   += weight * full_jastrow**2.d0     * gauss_a * gauss_b * coulomb
        numerical     += weight * full_jastrow**2.d0     * gauss_a * gauss_b 
        numerical_fit += weight * full_jastrow_fit**2.d0 * gauss_a * gauss_b 
       enddo
           
       de_abs = dabs(numerical- numerical_fit)
       abs_error_av += de_abs
       if(dabs(numerical).gt.1.d-10)then
        de_relat = de_abs/dabs(numerical)
       endif
       relat_error_av += de_relat
       de_abs = dabs( numerical - analytical)
       abs_error_fit_av += de_abs
       if(dabs(numerical).gt.1.d-10)then
        de_relat = de_abs/dabs(numerical)
       endif
       relat_error_fit_av += de_relat

       de_abs = dabs(analytical_c- numerical_c)
       abs_error_c_av += de_abs
       if(dabs(numerical_c).gt.1.d-10)then
        de_relat = de_abs/dabs(numerical_c)
       endif
       relat_error_c_av += de_relat

!      enddo
!     enddo
  enddo
 enddo


 abs_error_av   = abs_error_av/icount
 relat_error_av = relat_error_av /icount
 print*,'abs_error_av         = ',abs_error_av
 print*,'relat_error_av       = ',relat_error_av
 abs_error_fit_av   = abs_error_fit_av/icount
 relat_error_fit_av = relat_error_fit_av /icount
 print*,'abs_error_fit_av     = ',abs_error_fit_av
 print*,'relat_error_fit_av   = ',relat_error_fit_av
 abs_error_c_av   = abs_error_c_av/icount
 relat_error_c_av = relat_error_c_av /icount
 print*,'abs_error_c_av   = ',abs_error_c_av
 print*,'relat_error_c_av = ',relat_error_c_av

end


subroutine test_fit(mu)
 implicit none
 include 'utils/constants.include.F'
 double precision, intent(in) :: mu
 integer :: i,nx
 double precision :: x,dx,xmax,gama,delta,accu_fit,accu
 double precision :: jastrow_fit,jastrow
 double precision :: full_jastrow, full_jastrow_fit
 double precision :: f_tilde_fit,tilde_f_mu,shank_general
 double precision, allocatable :: array(:)
 double precision :: a0
 a0 = 1.d0/(2.d0*dsqrt(pi)*mu)
 integer :: n_taylor
 n_taylor = 4
 allocate(array(0:n_taylor))
 n_taylor = 4
 xmax = 5.d0
 nx = 100
 dx = xmax/dble(nx)
 x = 0.d0
 accu = 0.d0
 accu_fit = 0.d0
 do i = 1, nx
  jastrow_fit = f_tilde_fit(mu,x)
  full_jastrow_fit = dexp(jastrow_fit*a0)
  jastrow     = tilde_f_mu(mu,x)
  full_jastrow = dexp(jastrow*a0)
!  call exp_dl_rout(jastrow_fit,n_taylor, array)
!  full_jastrow_fit = shank_general(array,n_taylor,n_taylor)
  accu_fit += full_jastrow_fit * dx * dexp(-mu*x*x)
  accu     += full_jastrow     * dx * dexp(-mu*x*x)
  write(33,'(100(F16.10,X))')x,jastrow,jastrow_fit,full_jastrow,full_jastrow_fit
  x += dx
 enddo
 print*,'accu_fit = ',accu_fit
 print*,'accu     = ',accu
 print*,'er abs/rel ',dabs(accu_fit - accu),dabs(accu_fit - accu)/dabs(accu)
 


end

subroutine routine_test_n2_j
 implicit none
 include 'utils/constants.include.F'
 double precision :: r1(3), r2(3),weight,a0,coulomb
 double precision :: numerical,numerical_c,numerical_fit,n2_psi
 double precision :: jastrow_fit,f_tilde_fit,full_jastrow_fit,jastrow,tilde_f_mu,full_jastrow
 double precision :: int_ovlp_n2_jaswtrow2,int_erf_n2_jaswtrow2,analytical,analytical_c
 double precision :: norm,dm_a,dm_b,mu,r_12,numerical_c_fit
 double precision :: exponent_exp
 integer :: ipoint,n_taylor,istate
 exponent_exp = 1.d0
 istate = 1
 mu = 1.5d0
 a0 = 1.d0/(2.d0*dsqrt(pi)*mu)
 n_taylor = 4
 r1 = 0.1d0
 analytical   = int_ovlp_n2_jaswtrow2(r1,mu,istate,n_taylor,exponent_exp)
! analytical_c =  int_erf_n2_jaswtrow2(r1,mu,istate,n_taylor)
 norm = 0.d0
 numerical_c_fit = 0.d0
 numerical_fit   = 0.d0
 numerical       = 0.d0
 call dm_dft_alpha_beta_at_r(r1,dm_a,dm_b)
  do ipoint = 1, n_points_final_grid
   r2(1) = final_grid_points(1,ipoint)
   r2(2) = final_grid_points(2,ipoint)
   r2(3) = final_grid_points(3,ipoint)
   weight = final_weight_at_r_vector(ipoint)
   r_12     = dsqrt( (r1(1) - r2(1))**2.d0 + (r1(2) - r2(2))**2.d0 + (r1(3) - r2(3))**2.d0 )
   call give_n2_cas(r1,r2,istate,n2_psi)
   if(dabs(r_12).lt.1.d-6)then
    coulomb = 2.d0 * mu / sqpi - 2.d0 * mu**3 * r_12**2 / (3.d0 *sqpi) 
   else
    coulomb = derf(mu * r_12)/r_12
   endif
   jastrow_fit = f_tilde_fit(mu,r_12)
   full_jastrow_fit = dexp(jastrow_fit*a0) 
   jastrow = tilde_f_mu(mu,r_12)
   full_jastrow = dexp(jastrow*a0)
   numerical_c     += weight * full_jastrow**2.d0     *n2_psi * coulomb
   numerical       += weight * full_jastrow**2.d0     *n2_psi  
   numerical_fit   += weight * full_jastrow_fit**2.d0 *n2_psi 
   numerical_c_fit += weight * full_jastrow_fit**2.d0 *n2_psi * coulomb
   norm += weight * n2_psi
  enddo
  print*,'norm n2         = ',norm
  print*,'n(r1)/2         = ',(dm_a + dm_b)*0.5d0
  print*,''               
  print*,'n2*J num        = ',numerical
  print*,'n2*J num/fit    = ',numerical_fit
  print*,'n2*J analytical = ',analytical
  print*,'error fit       = ',dabs(numerical_fit-numerical),dabs(numerical_fit-numerical)/dabs(numerical)
  print*,'error analytic  = ',dabs(analytical-numerical),dabs(analytical-numerical)/dabs(numerical)
  print*,''               
  print*,'n2*J erf num    = ',numerical_c
  print*,'n2*J erf/num fit= ',numerical_c_fit
  print*,'n2*J erf/analyti= ',analytical_c
  print*,'error fit       = ',dabs(numerical_c_fit-numerical_c),dabs(numerical_c_fit-numerical_c)/dabs(numerical_c)
  print*,'error analytic  = ',dabs(analytical_c-numerical_c),dabs(analytical_c-numerical_c)/dabs(numerical_c)


end

subroutine routine_test_n2_inv_j
 implicit none
 include 'utils/constants.include.F'
 double precision :: r1(3), r2(3),weight,a0,coulomb
 double precision :: numerical,numerical_c,numerical_fit,n2_psi
 double precision :: jastrow_fit,f_tilde_fit,full_jastrow_fit,jastrow,tilde_f_mu,full_jastrow
 double precision :: int_ovlp_n2_jaswtrow2,int_erf_n2_jaswtrow2,analytical,analytical_c
 double precision :: norm,dm_a,dm_b,mu,r_12,numerical_c_fit
 double precision :: exponent_exp
 integer :: ipoint,n_taylor,istate
 exponent_exp = -1.d0
 istate = 1
 mu = 1.5d0
 a0 = 1.d0/(2.d0*dsqrt(pi)*mu)
 n_taylor = 4
 r1 = 0.1d0
 analytical   = int_ovlp_n2_jaswtrow2(r1,mu,istate,n_taylor,exponent_exp)
 analytical_c =  int_erf_n2_jaswtrow2(r1,mu,istate,n_taylor)
 norm = 0.d0
 numerical_c_fit = 0.d0
 numerical_fit   = 0.d0
 numerical       = 0.d0
 call dm_dft_alpha_beta_at_r(r1,dm_a,dm_b)
  do ipoint = 1, n_points_final_grid
   r2(1) = final_grid_points(1,ipoint)
   r2(2) = final_grid_points(2,ipoint)
   r2(3) = final_grid_points(3,ipoint)
   weight = final_weight_at_r_vector(ipoint)
   r_12     = dsqrt( (r1(1) - r2(1))**2.d0 + (r1(2) - r2(2))**2.d0 + (r1(3) - r2(3))**2.d0 )
   call give_n2_cas(r1,r2,istate,n2_psi)
   if(dabs(r_12).lt.1.d-6)then
    coulomb = 2.d0 * mu / sqpi - 2.d0 * mu**3 * r_12**2 / (3.d0 *sqpi) 
   else
    coulomb = derf(mu * r_12)/r_12
   endif
   jastrow_fit = f_tilde_fit(mu,r_12)
   full_jastrow_fit = dexp(-jastrow_fit*a0) 
   jastrow = tilde_f_mu(mu,r_12)
   full_jastrow = dexp(-jastrow*a0)
   numerical_c     += weight / full_jastrow**2.d0     *n2_psi * coulomb
   numerical       += weight / full_jastrow**2.d0     *n2_psi  
   numerical_fit   += weight / full_jastrow_fit**2.d0 *n2_psi 
   numerical_c_fit += weight / full_jastrow_fit**2.d0 *n2_psi * coulomb
   norm += weight * n2_psi
  enddo
  print*,'norm n2         = ',norm
  print*,'n(r1)/2         = ',(dm_a + dm_b)*0.5d0
  print*,''               
  print*,'n2*J num        = ',numerical
  print*,'n2*J num/fit    = ',numerical_fit
  print*,'n2*J analytical = ',analytical
  print*,'error fit       = ',dabs(numerical_fit-numerical),dabs(numerical_fit-numerical)/dabs(numerical)
  print*,'error analytic  = ',dabs(analytical-numerical),dabs(analytical-numerical)/dabs(numerical)
  print*,''               
  print*,'n2*J erf num    = ',numerical_c
  print*,'n2*J erf/num fit= ',numerical_c_fit
  print*,'n2*J erf/analyti= ',analytical_c
  print*,'error fit       = ',dabs(numerical_c_fit-numerical_c),dabs(numerical_c_fit-numerical_c)/dabs(numerical_c)
  print*,'error analytic  = ',dabs(analytical_c-numerical_c),dabs(analytical_c-numerical_c)/dabs(numerical_c)


end

subroutine routine_test_n2_j_full
 implicit none
 include 'utils/constants.include.F'
 double precision :: r1(3),weight,a0,coulomb
 double precision :: numerical,numerical_c,n2_psi
 double precision :: jastrow_fit,f_tilde_fit,full_jastrow_fit,jastrow,tilde_f_mu,full_jastrow
 double precision :: int_ovlp_n2_jaswtrow2,int_erf_n2_jaswtrow2,analytical,analytical_c
 double precision :: norm,dm_a,dm_b,muj,muc,r_12,numerical_c_renorm,dens,numerical_renorm,exponent_exp
 integer :: ipoint,n_taylor,istate
 exponent_exp = 1.d0
 istate = 1
 n_taylor = 4
 numerical          = 0.d0
 numerical_renorm   = 0.d0
 numerical_c        = 0.d0
 numerical_c_renorm = 0.d0
 muc = 1500.d0
 print*,'n_points_final_grid = ',n_points_final_grid
 pause
 do ipoint = 1, n_points_final_grid
 print*,'ipoint = ',ipoint
  r1(1) = final_grid_points(1,ipoint)
  r1(2) = final_grid_points(2,ipoint)
  r1(3) = final_grid_points(3,ipoint)
  dm_a = one_e_dm_and_grad_alpha_in_r(4,ipoint,istate)
  dm_b =  one_e_dm_and_grad_beta_in_r(4,ipoint,istate)
  dens = 0.5d0 * (dm_a + dm_b)
  muj = mu_of_r_prov(ipoint,istate)
!  muc = mu_of_r_prov(ipoint,istate)
!  muj = 1500.d0
  analytical   = int_ovlp_n2_jaswtrow2(r1,muj,istate,n_taylor,exponent_exp)
  analytical_c =  int_erf_n2_jaswtrow2(r1,muj,muc,istate,n_taylor)
  weight = final_weight_at_r_vector(ipoint)
  numerical_c_renorm += weight * analytical_c * dens/analytical
  numerical_renorm   += weight * dens
  numerical_c     += weight * analytical_c
  numerical       += weight * analytical
 enddo
 print*,'n2*J num/fit           = ',numerical
 print*,'norm_n2_jastrow        = ',norm_n2_jastrow
 print*,'Na*Nb                  = ',elec_alpha_num * elec_beta_num
 print*,''                     
 print*,'n2*J erf num           = ',numerical_c
 print*,'n2*J erf num/ren       = ',numerical_c_renorm
 print*,'************************************************'
 print*,'************************************************'
 print*,'psi_energy_two_e       = ',psi_energy_two_e
 print*,'psi_wee_mu_of_r        = ',psi_wee_mu_of_r
 print*,'---------> '
 print*,'psi_wee_mu_of_r_sr     = ',psi_wee_mu_of_r_sr
 print*,''
 print*,'coulomb_n2_jastrow     = ',coulomb_n2_jastrow
 print*,'coulomb_n2_jastrow_reno= ',coulomb_n2_jastrow_renorm
 print*,'wee_mu_of_r_n2_jastrow = ',wee_mu_of_r_n2_jastrow
 print*,'wee_mu_of_r_n2_jastrow_r ',wee_mu_of_r_n2_jastrow_renorm
 print*,'---------> '
 print*,'wee_mu_of_r_sr_n2_jastrow',wee_mu_of_r_sr_n2_jastrow
 print*,'wee_mu_of_r_sr_n2Rjastrow',wee_mu_of_r_sr_n2_jastrow_renorm
 print*,''
 print*,''
 print*,''


end
