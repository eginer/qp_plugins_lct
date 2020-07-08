subroutine big_thing
 implicit none 
 double precision :: weight1,r1(3),weight2,r2(3),int_ao
 double precision :: ao_two_e_integral_schwartz_accel_gauss
 double precision :: int_r2(3),int_gauss_num,alpha_r12,coef,r12,int_r2_bis(3)
 integer :: ipoint,i,j,n_pt_in,jpoint,m
 double precision :: aos_array_r1(ao_num),aos_array_r2(ao_num),aos_grad_array_r1(3,ao_num),aos_grad_array_r2(3,ao_num)
 double precision :: accu_relat(3),accu_abs(3),err_relat,err_abs
 double precision :: accu_tmp(3),d_dr12(3),mu_in,derf_mu_x,accu_tmp_bis(3)
 include 'utils/constants.include.F'
 integer :: iao,jao,kao,lao
 double precision :: test, ao_two_e_integral_schwartz_accel_erf,num_int
 mu_in = mu_erf
 accu_abs = 0.d0
 accu_relat = 0.d0
 num_int = 0.d0
 do jao = 1, ao_num ! r1
  do lao = 1, ao_num ! r2
   do iao = 1, ao_num ! r1
    do kao = 1, ao_num ! r2
    num_int += 1.d0
     print*,'<ij|kl> = ',jao,lao,iao,kao
     call ao_two_e_d_dr12_int(iao,jao,kao,lao,mu_in,d_dr12)
     int_gauss_num = 0.d0
     accu_tmp = 0.d0
     accu_tmp_bis = 0.d0
     do ipoint = 1, n_points_final_grid
      r1(1) = final_grid_points(1,ipoint)
      r1(2) = final_grid_points(2,ipoint)
      r1(3) = final_grid_points(3,ipoint)
      call give_all_aos_and_grad_at_r(r1,aos_array_r1,aos_grad_array_r1)
      weight1 = final_weight_at_r_vector(ipoint)
      int_r2 = 0.d0
      int_r2_bis = 0.d0
      do jpoint = 1, n_points_final_grid
       r2(1) = final_grid_points(1,jpoint)
       r2(2) = final_grid_points(2,jpoint)
       r2(3) = final_grid_points(3,jpoint)
       weight2 = final_weight_at_r_vector(jpoint)
       call give_all_aos_and_grad_at_r(r2,aos_array_r2,aos_grad_array_r2)
       r12 = (r1(1) - r2(1))**2.d0 + (r1(2) - r2(2))**2.d0 + (r1(3) - r2(3))**2.d0 
       r12 = dsqrt(r12)
       do m = 1, 3

        int_r2_bis(m) += weight2 * 0.5d0 * derf_mu_x(mu_in,r12) & 
                             * (r1(m) - r2(m)) *  aos_array_r1(jao) * aos_array_r2(lao)  & 
                             * (aos_grad_array_r1(m,iao) * aos_array_r2(kao) - aos_grad_array_r2(m,kao) * aos_array_r1(iao))

!       !!!!!!!!! x1 dx1 i(r1)
        int_r2(m) += weight2 * 0.5d0 * derf_mu_x(mu_in,r12)      & 
                             * r1(m) * aos_grad_array_r1(m,iao)  & 
                             *  aos_array_r2(kao)                & 
                             *  aos_array_r1(jao) * aos_array_r2(lao)  
       !!!!!!!!! x2 dx2 k(r2)
        int_r2(m) += weight2 * 0.5d0 * derf_mu_x(mu_in,r12)      & 
                             * r2(m) * aos_grad_array_r2(m,kao)  & 
                             *  aos_array_r1(iao)                & 
                             *  aos_array_r1(jao) * aos_array_r2(lao)  
!       !!!!!!!!! x1 i(r1) dx2 k(r1)
        int_r2(m) -= weight2 * 0.5d0 * derf_mu_x(mu_in,r12)      & 
                             * r1(m) * aos_array_r1(iao)  & 
                             *  aos_grad_array_r2(m,kao)                & 
                             *  aos_array_r1(jao) * aos_array_r2(lao)  
!       !!!!!!!!! x2 k(r2) dx1 i(r1)
        int_r2(m) -= weight2 * 0.5d0 * derf_mu_x(mu_in,r12)      & 
                             * r2(m) * aos_array_r2(kao)  & 
                             *  aos_grad_array_r1(m,iao)                & 
                             *  aos_array_r1(jao) * aos_array_r2(lao)  
       enddo
      enddo
      do m = 1, 3
       accu_tmp(m) += weight1 * int_r2(m)
       accu_tmp_bis(m) += weight1 * int_r2_bis(m)
      enddo
     enddo
     do m = 1, 3
      int_gauss_num = accu_tmp(m)
      int_ao = d_dr12(m)
      err_abs = dabs(int_gauss_num - int_ao)
      if(dabs(int_gauss_num).gt.1.d-10)then
       err_relat = err_abs/dabs(int_gauss_num)
      else
       err_relat = 0.d0
      endif
      print*,'m = ',m
      print*,'int_gauss_num = ',int_gauss_num
      print*,'accu_tmp_bis(m)=',accu_tmp_bis(m)
      print*,'int_ao        = ',int_ao
      print*,'abs error     = ',err_abs
      print*,'err_relat     = ',err_relat
      if(err_relat .gt. 1.d-2)then
       print*,'AHAHAHAAH'
       stop
      endif
      accu_abs(m) += err_abs
      accu_relat(m) += err_relat
     enddo
    enddo
   enddo
  enddo
 enddo
 print*,''
 print*,''
 print*,''
 print*,'Summary'
 print*,''
 print*,''
 print*,''
 print*,'accu_abs   = ',accu_abs/num_int
 print*,'accu_relat = ',accu_relat/num_int
end

double precision function derf_mu_x(mu,x)
 implicit none
 include 'utils/constants.include.F'
 double precision, intent(in) :: mu,x
  if(dabs(x).gt.1.d-6)then
   derf_mu_x = derf(mu * x)/x
  else
   derf_mu_x =  inv_sq_pi * 2.d0 * mu * (1.d0 - mu*mu*x*x/3.d0)
  endif

end

subroutine test_hermit
 implicit none 
 include 'utils/constants.include.F'
 integer :: iao,jao,kao,lao
 double precision :: mu_in,d_dr12(3),num_int_ik,num_int_jl
 mu_in = mu_erf
 do jao = 1, ao_num ! r1
  do lao = 1, ao_num ! r2
   do iao = 1, ao_num ! r1
    do kao = 1, ao_num ! r2
     if(iao == jao)cycle
     if(kao == lao)cycle
     !                        <jl|ik>
     call ao_two_e_d_dr12_int(iao,jao,kao,lao,mu_in,d_dr12)
     num_int_ik = d_dr12(1)
     num_int_ik+= d_dr12(2)
     num_int_ik+= d_dr12(3)
     if(dabs(num_int_ik).lt.1.d-10)cycle
     print*,'******************'
     print*,'j,l,i,k',jao,lao,iao,kao
     print*,'<jl|ik>'
     print*,num_int_ik
     !                        <lj|ki>
     call ao_two_e_d_dr12_int(kao,lao,iao,jao,mu_in,d_dr12)
     num_int_ik = d_dr12(1)
     num_int_ik+= d_dr12(2)
     num_int_ik+= d_dr12(3)
     print*,'<jl|ki>'
     print*,num_int_ik
    enddo
   enddo
  enddo
 enddo
 
end


subroutine new_test_big_thing
 implicit none
 double precision :: weight,r1(3),mu_large,d_dr12_large(3),d_dr12_tmp(3),d_dr12_large_tmp(3)
 integer :: ipoint,i,j,n_pt_in,jpoint,m
 double precision :: int_r2(3,ao_num,ao_num), int_r2_large(3, ao_num,ao_num),xyz_ints_large(3, ao_num,ao_num)
 double precision :: mu_in,aos_array_r1(ao_num),aos_grad_array_r1(3,ao_num),phi_j_erf_mu_r_phi
 double precision :: xyz_ints(3, ao_num,ao_num), coulomb(ao_num,ao_num), coulomb_large(ao_num,ao_num)
 double precision :: d_dr12(3),int_ao,int_gauss_num,err_relat,err_abs,accu_relat(3),accu_abs(3)
 integer :: iao,jao,kao,lao
 mu_in = mu_erf
 mu_large = 1.d+9
 do kao = 1, ao_num ! r1
  do lao = 1, ao_num ! r2
    int_r2 = 0.d0
    int_r2_large = 0.d0
    do ipoint = 1, n_points_final_grid
     r1(1) = final_grid_points(1,ipoint)
     r1(2) = final_grid_points(2,ipoint)
     r1(3) = final_grid_points(3,ipoint)
     call give_all_aos_and_grad_at_r(r1,aos_array_r1,aos_grad_array_r1)
     weight = final_weight_at_r_vector(ipoint)
     do iao = 1, ao_num ! r1
      do jao = 1, ao_num ! r2
       call phi_j_erf_mu_r_xyz_phi(iao,jao,mu_in, r1, xyz_ints(1, jao,iao))
       coulomb(jao,iao) = phi_j_erf_mu_r_phi(iao,jao,mu_in,r1)
       call phi_j_erf_mu_r_xyz_phi(iao,jao,mu_large, r1, xyz_ints_large(1, jao,iao))
       coulomb_large(jao,iao) = phi_j_erf_mu_r_phi(iao,jao,mu_large,r1)
      enddo
     enddo
     do jao = 1, ao_num
      do iao = 1, ao_num
       do m = 1, 3
        int_r2(m,iao,jao) += weight * aos_array_r1(kao) * r1(m) * aos_grad_array_r1(m,iao) * coulomb(jao,lao)    
        int_r2(m,iao,jao) += weight * aos_array_r1(lao) * r1(m) * aos_grad_array_r1(m,jao) * coulomb(iao,kao)    
        int_r2(m,iao,jao) -= weight * aos_array_r1(lao)         * aos_grad_array_r1(m,jao) * xyz_ints(m,iao,kao) 
        int_r2(m,iao,jao) -= weight * aos_array_r1(kao)         * aos_grad_array_r1(m,iao) * xyz_ints(m,jao,lao) 

        int_r2_large(m,iao,jao) += weight * aos_array_r1(kao) * r1(m) * aos_grad_array_r1(m,iao) * coulomb_large(jao,lao)    
        int_r2_large(m,iao,jao) += weight * aos_array_r1(lao) * r1(m) * aos_grad_array_r1(m,jao) * coulomb_large(iao,kao)    
        int_r2_large(m,iao,jao) -= weight * aos_array_r1(lao)         * aos_grad_array_r1(m,jao) * xyz_ints_large(m,iao,kao) 
        int_r2_large(m,iao,jao) -= weight * aos_array_r1(kao)         * aos_grad_array_r1(m,iao) * xyz_ints_large(m,jao,lao) 
       enddo
      enddo
     enddo
    enddo
    int_r2 *= 0.5d0
    int_r2_large *= 0.5d0

    do jao = 1, ao_num
     do iao = 1, ao_num
      call ao_two_e_eff_dr12_pot(iao,kao,jao,lao,mu_in,d_dr12,d_dr12_large)
      print*,'iao,kao,jao,lao',iao,kao,jao,lao
      do m = 1, 3
       int_ao = d_dr12(m)  - d_dr12_large(m)! <kl|ij>
       int_gauss_num = int_r2(m,iao,jao) - int_r2_large(m,iao,jao) ! <kl|ij>
       err_abs = dabs(int_gauss_num - int_ao)
       if(dabs(int_gauss_num).gt.1.d-10)then
        err_relat = err_abs/dabs(int_gauss_num)
       else
        err_relat = 0.d0
       endif
       print*,'m = ',m
       print*,'int_gauss_num = ',int_gauss_num
       print*,'int_ao        = ',int_ao
       print*,'abs error     = ',err_abs
       print*,'err_relat     = ',err_relat
       if(err_relat .gt. 1.d-2 .and. dabs(int_gauss_num).gt.1.d-7)then
        print*,'AHAHAHAAH'
!        stop
       endif
       accu_abs(m) += err_abs
       accu_relat(m) += err_relat
      enddo
     enddo
    enddo
  enddo
 enddo

 print*,''
 print*,''
 print*,''
 print*,'Summary'
 print*,''
 print*,''
 print*,''
 print*,'accu_abs   = ',accu_abs/(dble(ao_num)**4.d0)
 print*,'accu_relat = ',accu_relat/(dble(ao_num)**4.d0)

end
