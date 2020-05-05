
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
 double precision :: de_abs,de_relat,erf_exp_f_phi_ij
 double precision, allocatable :: array(:)
 integer :: n_taylor
 double precision :: a0
 n_taylor = 4
 allocate(array(0:n_taylor))
 mu = 1.5d0
 a0 = 1.d0/(2.d0*dsqrt(pi)*mu)
 call test_fit(mu)
 r1 = 0.1d0
! r1(1) = 0.1d0

 abs_error_av   = 0.d0
 relat_error_av = 0.d0
 abs_error_c_av   = 0.d0
 relat_error_c_av = 0.d0
 icount = 0.d0
 do i = 1, ao_num
  do j = 1, ao_num
   print*,'i,j,mu',i,j,mu
    num_A = ao_nucl(i)
    power_A(1:3)= ao_power(i,1:3)
    A_center(1:3) = nucl_coord(num_A,1:3)
   
    num_B = ao_nucl(j)
    power_B(1:3)= ao_power(j,1:3)
    B_center(1:3) = nucl_coord(num_B,1:3)
    do k = 1, ao_prim_num(i)
     alpha = ao_expo_ordered_transp(k,i)     
     do l = 1, ao_prim_num(j)
      beta = ao_expo_ordered_transp(l,j)     
       ! analytical integral 
       analytical   = ovlp_exp_f_phi_ij(mu,r1,A_center,B_center,power_A,power_B,alpha,beta,n_taylor)
       analytical_c =  erf_exp_f_phi_ij(mu,r1,A_center,B_center,power_A,power_B,alpha,beta,n_taylor)

       ! numerical  integral 
       numerical   = 0.d0
       numerical_c = 0.d0
       numerical_fit = 0.d0
       icount += 1.d0
       do ipoint = 1, n_points_final_grid
        r2(1) = final_grid_points(1,ipoint)
        r2(2) = final_grid_points(2,ipoint)
        r2(3) = final_grid_points(3,ipoint)
        weight = final_weight_at_r_vector(ipoint)

        gauss_a = primitive_value_explicit(power_A,A_center,alpha,r2)
        gauss_b = primitive_value_explicit(power_B,B_center,beta ,r2)

        r_12     = dsqrt( (r1(1) - r2(1))**2.d0 + (r1(2) - r2(2))**2.d0 + (r1(3) - r2(3))**2.d0 )
        if(dabs(r_12).lt.1.d-6)then
         coulomb = 2.d0 * mu / sqpi - 2.d0 * mu**3 * r_12**2 / (3.d0 *sqpi) 
        else
         coulomb = derf(mu * r_12)/r_12
        endif
        jastrow_fit = f_tilde_fit(mu,r_12)
        full_jastrow_fit = dexp(jastrow_fit*a0) 
        jastrow = tilde_f_mu(mu,r_12)
        full_jastrow = dexp(jastrow*a0)
        numerical_c += weight * full_jastrow * gauss_a * gauss_b * coulomb
        numerical +=     weight * full_jastrow     * gauss_a * gauss_b 
        numerical_fit += weight * full_jastrow_fit * gauss_a * gauss_b 
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
      enddo
     enddo
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
