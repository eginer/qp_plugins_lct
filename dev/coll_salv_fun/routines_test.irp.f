subroutine test_int_f_tilde
 implicit none
 BEGIN_DOC
! test the function ovlp_f_tilde_phi_ij
 END_DOC
 include 'utils/constants.include.F'
 double precision :: A_center(3), B_center(3), alpha,beta      ! basis functions
 integer          :: power_A(3), power_B(3)                    ! basis functions 
 double precision :: r1(3),mu

 double precision :: ovlp_f_tilde_phi_ij
 double precision :: numerical,  analytical, weight,primitive_value_explicit
 double precision :: numerical_c,analytical_c,r_12,r2(3)
 double precision :: gauss_a, gauss_b, gauss_d, coulomb, jastrow
 integer          :: ipoint,ao_i,ao_j,num_A,num_B,i,j,k,l
 double precision :: abs_error_av, relat_error_av,icount
 double precision :: abs_error_c_av, relat_error_c_av,slater_fit_gam
 double precision :: exp_dl,f_tilde_fit,full_jastrow
 integer :: n_taylor
 n_taylor = 10
 mu = 0.5d0
 r1 = 0.d0
! r1(1) = 0.1d0

 abs_error_av   = 0.d0
 relat_error_av = 0.d0
 abs_error_c_av   = 0.d0
 relat_error_c_av = 0.d0
 icount = 0.d0
 do i = 1, ao_num
  do j = 1, ao_num
   print*,'i,j',i,j
    num_A = ao_nucl(i)
    power_A(1:3)= ao_power(i,1:3)
    A_center(1:3) = nucl_coord(num_A,1:3)
   
    num_B = ao_nucl(j)
    power_B(1:3)= ao_power(j,1:3)
    B_center(1:3) = nucl_coord(num_B,1:3)
    do l=1,ao_prim_num(i)
     alpha = ao_expo_ordered_transp(l,i)     
     do k=1,ao_prim_num(j)
      beta = ao_expo_ordered_transp(k,j)     
   
       ! analytical integral 
       
       analytical = ovlp_f_tilde_phi_ij(mu,r1,A_center,B_center,power_A,power_B,alpha,beta,n_taylor)

       ! numerical  integral 
       numerical   = 0.d0
       numerical_c = 0.d0
       icount += 1.d0
       do ipoint = 1, n_points_final_grid
        r2(1) = final_grid_points(1,ipoint)
        r2(2) = final_grid_points(2,ipoint)
        r2(3) = final_grid_points(3,ipoint)
        weight = final_weight_at_r_vector(ipoint)
        r_12     = dsqrt( (r1(1) - r2(1))**2 + (r1(2) - r2(2))**2 + (r1(3) - r2(3))**2 )
        if(dabs(r_12).lt.1.d-6)then
         coulomb = 2.d0 * mu / sqpi - 2.d0 * mu**3 * r_12**2 / (3.d0 *sqpi) 
        else
         coulomb = derf(mu * r_12)/r_12
        endif
        ! fit of f_tilde with h(r_12)
        jastrow = f_tilde_fit(mu,r_12)
        full_jastrow = exp_dl(-jastrow, n_taylor)
!        full_jastrow = dexp(-jastrow)
        gauss_a = primitive_value_explicit(power_A,A_center,alpha,r2)
        gauss_b = primitive_value_explicit(power_B,B_center,beta ,r2)
        numerical += weight * full_jastrow * gauss_a * gauss_b 
!        numerical_c += weight * full_jastrow * gauss_a * gauss_b * coulomb
       enddo
       abs_error_av += dabs(analytical - numerical)
       if(dabs(numerical).gt.1.d-10)then
        relat_error_av += dabs(analytical - numerical)/dabs(numerical)
        if(dabs(analytical - numerical)/dabs(numerical) .gt. 1.d-5 )then
         print*,'i,j',i,j
         print*,'l,k',l,k
         print*,'power_A = ',power_A
         print*,'power_B = ',power_B
         print*,'alpha, beta', alpha, beta
         print*,'numerical, analytical ',numerical,analytical
         print*,'error      = ',dabs(analytical - numerical),dabs(analytical - numerical)/dabs(numerical)
        endif
       endif 

!       abs_error_c_av += dabs(analytical_c - numerical_c)
!       if(dabs(numerical_c).gt.1.d-10)then
!        relat_error_c_av += dabs(analytical_c - numerical_c)/dabs(numerical_c)
!        if(dabs(analytical_c - numerical_c)/dabs(numerical_c) .gt. 1.d-6 )then
!         print*,'i,j',i,j
!         print*,'l,k',l,k
!         print*,'power_A = ',power_A
!         print*,'power_B = ',power_B
!         print*,'alpha, beta', alpha, beta
!         print*,'numerical_c, analytical_c ',numerical_c,analytical_c
!         print*,'error      = ',dabs(analytical_c - numerical_c),dabs(analytical_c - numerical_c)/dabs(numerical_c)
!        endif
!       endif 

     enddo  ! k
    enddo ! l
  enddo ! j 
 enddo ! i
 abs_error_av   = abs_error_av/icount
 relat_error_av = relat_error_av /icount
 print*,'abs_error_av     = ',abs_error_av
 print*,'relat_error_av   = ',relat_error_av
! abs_error_c_av   = abs_error_c_av/icount
! relat_error_c_av = relat_error_c_av /icount
! print*,'abs_error_c_av   = ',abs_error_c_av
! print*,'relat_error_c_av = ',relat_error_c_av

end
