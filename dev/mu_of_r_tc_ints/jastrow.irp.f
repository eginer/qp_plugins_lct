double precision function jastrow_mu(r1,r2,mu)
 double precision, intent(in) :: r1(3), r2(3), mu
 include 'constants.include.F'
 double precision :: r12
 r12 = 0.d0
 do m = 1, 3
  r12 += (r1(m) - r2(m)) * (r1(m) - r2(m))
 enddo
 r12 = dsqrt(r12)
 jastrow_mu = 0.5d0 * r12 * (1.d0 - derf(mu*r12)) - 0.5d0 * inv_sq_pi / mu * dexp(-mu*mu*r12*r12)

end


subroutine grad_r1_jastrow_mu(r1,r2,mu,grad_mu,jastrow,grad_jastrow)
 double precision, intent(in) :: mu,grad_mu(3),r1(3),r2(3)
 double precision, intent(out) :: jastrow,grad_jastrow(3)
 double precision :: r12_vec(3),gauss,erfc_mu,r12
 integer :: m
 include 'constants.include.F'
 r12 = 0.d0
 do m = 1, 3
  r12_vec(m) = (r1(m) - r2(m))
  r12 += r12_vec(m) * r12_vec(m)
 enddo
 r12 = dsqrt(r12)
 gauss = dexp(-mu*mu*r12*r12)
 erfc_mu = 1.d0 - derf(mu*r12) 
 jastrow_mu = 0.5d0 * r12 * erfc_mu - 0.5d0 * inv_sq_pi * gauss 
 grad_jastrow = 0.d0
 if(r12.lt.1.d-6)return
 do m = 1, 3
  grad_jastrow(m) = 0.5d0 * r12_vec(m)/r12 * erfc_mu + 0.5d0 * inv_sq_pi /(mu*mu) * gauss * grad_mu(m)
 enddo
end

subroutine grad_r2_jastrow_mu(r1,r2,mu,jastrow,grad_jastrow)
 double precision, intent(in) :: mu,r1(3),r2(3)
 double precision, intent(out) :: jastrow,grad_jastrow(3)
 double precision :: r12_vec(3),gauss,erfc_mu,r12
 integer :: m
 include 'constants.include.F'
 r12 = 0.d0
 do m = 1, 3
  r12_vec(m) = (r1(m) - r2(m))
  r12 += r12_vec(m) * r12_vec(m)
 enddo
 r12 = dsqrt(r12)
 gauss = dexp(-mu*mu*r12*r12)
 erfc_mu = 1.d0 - derf(mu*r12) 
 jastrow_mu = 0.5d0 * r12 * erfc_mu - 0.5d0 * inv_sq_pi * gauss 
 grad_jastrow = 0.d0
 if(r12.lt.1.d-6)return
 do m = 1, 3
  grad_jastrow(m) = -0.5d0 * r12_vec(m)/r12 * erfc_mu 
 enddo
end

subroutine test_grad_jastrow
 implicit none
 include 'constants.include.F'
 integer :: ipoint,i,j,k,l,m,jpoint
 double precision :: r1(3),r2(3),weight1,weight2,weight_tot,ao_prod_r2,ao_prod_r1,mu_min,r12
 double precision :: mu_lda_erf,rho_a_hf,rho_b_hf, mu, mu_plus, mu_minus, grad_mu(3),dx
 double precision :: jastrow,grad_jastrow(3),jastrow_plus,jastrow_minus, grad_jastrow_num(3)
 double precision :: accu(3),jastrow_mu,r1_scal,r2_scal
 do l = 7, 7
  dx = 10.d0**(-l)
  mu_min = mu_erf
  accu = 0.d0
 do jpoint = 1, n_points_final_grid ! r2
!  do jpoint = 91, 91
   weight2 = final_weight_at_r_vector(jpoint)
   r2(:) = final_grid_points(:,jpoint)
   r2_scal = dsqrt(r2(1)**2+r2(2)**2+r2(3)**2)
  do ipoint = 1, n_points_final_grid ! r1
!   do ipoint = 43,43
    weight1 = final_weight_at_r_vector(ipoint)
    r1_scal = dsqrt(r1(1)**2+r1(2)**2+r1(3)**2)
    do m = 1, 3
     r1(:) = final_grid_points(:,ipoint) 
     r1(m) += dx 
     call dm_dft_alpha_beta_at_r(r1,rho_a_hf,rho_b_hf)
     mu_plus = mu_lda_erf(rho_a_hf,rho_b_hf,mu_min)
     r1(:) = final_grid_points(:,ipoint) 
     r1(m) -= dx 
     call dm_dft_alpha_beta_at_r(r1,rho_a_hf,rho_b_hf)
     mu_minus = mu_lda_erf(rho_a_hf,rho_b_hf,mu_min)
     grad_mu(m) = (mu_plus - mu_minus)/(2.d0 * dx)
    enddo
    r1(:) = final_grid_points(:,ipoint) 
    call dm_dft_alpha_beta_at_r(r1,rho_a_hf,rho_b_hf)
    mu = mu_lda_erf(rho_a_hf,rho_b_hf,mu_min)
    call grad_r1_jastrow_mu(r1,r2,mu,grad_mu,jastrow,grad_jastrow)
    do m = 1, 3
     r1(:) = final_grid_points(:,ipoint) 
     r1(m) += dx 
     call dm_dft_alpha_beta_at_r(r1,rho_a_hf,rho_b_hf)
     mu_plus = mu_lda_erf(rho_a_hf,rho_b_hf,mu_min)
     jastrow_plus = jastrow_mu(r1,r2,mu_plus)
 
     r1(:) = final_grid_points(:,ipoint) 
     r1(m) -= dx 
     call dm_dft_alpha_beta_at_r(r1,rho_a_hf,rho_b_hf)
     mu_minus = mu_lda_erf(rho_a_hf,rho_b_hf,mu_min)
     jastrow_minus = jastrow_mu(r1,r2,mu_minus)
     grad_jastrow_num(m) = (jastrow_plus - jastrow_minus) / (2.d0 * dx)
    enddo
    do m = 1, 3
     accu(m) += dabs(grad_jastrow_num(m) - grad_jastrow(m)) * weight2 * weight1 * dexp(-r1_scal*r1_scal) * dexp(-r2_scal*r2_scal)
     if(dabs(grad_jastrow(m)).gt.1.d-7)then
      if(dabs(grad_jastrow_num(m) - grad_jastrow(m))/dabs(grad_jastrow(m)).gt.1.d-1)then
       print*,'' 
       print*,dx,m
       print*,ipoint,jpoint
       print*,grad_jastrow_num(m),grad_jastrow(m),dabs(grad_jastrow_num(m) - grad_jastrow(m))
      endif
     endif
    enddo
!  
   enddo
  enddo
  print*,'accu = '
  print*,accu 
 enddo


end


subroutine test_grad_jastrow_sq
 implicit none
 include 'constants.include.F'
 integer :: ipoint,i,j,k,l,m,jpoint
 double precision :: r1(3),r2(3),weight1,weight2,weight_tot,ao_prod_r2,ao_prod_r1,mu_min,r12
 double precision :: mu_lda_erf,rho_a_hf,rho_b_hf, mu, mu_plus, mu_minus, grad_mu(3),dx
 double precision :: jastrow,grad_jastrow(3),jastrow_plus,jastrow_minus, grad_jastrow_sq_anal
 double precision :: accu,jastrow_mu,r1_scal,r2_scal, grad_r1_jastrow_sq
 double precision :: nabla_r12_bis,nabla_sq_term_general,erf_mu_sq_general
 double precision :: accu1,accu2
 do l = 7, 7
  dx = 10.d0**(-l)
  mu_min = mu_erf
  accu = 0.d0
  accu1 = 0.d0
  accu2 = 0.d0
  do jpoint = 1, n_points_final_grid ! r2
   weight2 = final_weight_at_r_vector(jpoint)
   r2(:) = final_grid_points(:,jpoint)
   r2_scal = dsqrt(r2(1)**2+r2(2)**2+r2(3)**2)
   do ipoint = 1, n_points_final_grid ! r1
     weight1 = final_weight_at_r_vector(ipoint)
     r1_scal = dsqrt(r1(1)**2+r1(2)**2+r1(3)**2)
     r12 = 0.d0
     do m = 1, 3 ! compute grad mu
      r1(:) = final_grid_points(:,ipoint) 
      r12 += (r1(m) - r2(m))*(r1(m) - r2(m))
     enddo
     r12 = dsqrt(r12)
     do m = 1, 3 ! compute grad mu
      r1(:) = final_grid_points(:,ipoint) 
      r1(m) += dx 
      call dm_dft_alpha_beta_at_r(r1,rho_a_hf,rho_b_hf)
      mu_plus = mu_lda_erf(rho_a_hf,rho_b_hf,mu_min)
      r1(:) = final_grid_points(:,ipoint) 
      r1(m) -= dx 
      call dm_dft_alpha_beta_at_r(r1,rho_a_hf,rho_b_hf)
      mu_minus = mu_lda_erf(rho_a_hf,rho_b_hf,mu_min)
      grad_mu(m) = (mu_plus - mu_minus)/(2.d0 * dx)
     enddo
!     grad_mu = 0.d0
     if(r12.lt.1.d-8)cycle
     r1(:) = final_grid_points(:,ipoint) 
     r1_scal = dsqrt(r1(1)**2+r1(2)**2+r1(3)**2)
     call dm_dft_alpha_beta_at_r(r1,rho_a_hf,rho_b_hf)
     mu = mu_lda_erf(rho_a_hf,rho_b_hf,mu_min)
!     mu = mu_erf 
     call grad_r1_jastrow_mu(r1,r2,mu,grad_mu,jastrow,grad_jastrow)
     grad_r1_jastrow_sq = 0.d0
     do m = 1, 3
      grad_r1_jastrow_sq -= 0.5d0 * grad_jastrow(m)*grad_jastrow(m)
     enddo
     call grad_r2_jastrow_mu(r1,r2,mu,jastrow,grad_jastrow)
     do m = 1, 3
      grad_r1_jastrow_sq -= 0.5d0 * grad_jastrow(m)*grad_jastrow(m)
     enddo
     grad_jastrow_sq_anal = nabla_sq_term_general(r1,r2,mu,grad_mu) & 
                          + erf_mu_sq_general(r1,r2,mu)             &  
                          + nabla_r12_bis(r1,r2,mu,grad_mu)         
     accu += dabs(grad_r1_jastrow_sq - grad_jastrow_sq_anal) * weight1 * weight2 * dexp(-r1_scal*r1_scal) * dexp(-r2_scal*r2_scal)
     accu1 += grad_r1_jastrow_sq   * weight1 * weight2 * dexp(-r1_scal*r1_scal) * dexp(-r2_scal*r2_scal)
     accu2 += grad_jastrow_sq_anal * weight1 * weight2 * dexp(-r1_scal*r1_scal) * dexp(-r2_scal*r2_scal)
   enddo
  enddo
 enddo
 print*,'difference = ',accu
 print*,'accu num   = ',accu1
 print*,'accu analyt= ',accu2
end

double precision function nabla_sq_term_general(r1,r2,mu,grad_mu)
 implicit none
 BEGIN_DOC
! -1/(8 pi * mu(r1)^4) (\nabla_1 \mu(r1))^2 exp(-2 (\mu(r1)r12)^2)
 END_DOC
 double precision, intent(in) :: r1(3), r2(3), mu, grad_mu(3)
 double precision:: cst
 include 'constants.include.F'
 double precision :: r12,grad_mu_sq
 integer :: m
 cst = -0.125d0 * inv_pi
 r12 = 0.d0
 grad_mu_sq = 0.d0
 do m = 1, m
  r12  += (r1(m) - r2(m))*(r1(m) - r2(m))
  grad_mu_sq += grad_mu(m)*grad_mu(m)
 enddo
 nabla_sq_term_general = cst * 1/mu**4 * dexp(-2.d0 * mu * mu * r12) * grad_mu_sq

end

double precision function erf_mu_sq_general(r1,r2,mu)
 implicit none
 BEGIN_DOC
! -1/4 (erf(mu(r1) r12) - 1)^2
 END_DOC
 double precision, intent(in) :: r1(3), r2(3), mu
 double precision :: cst
 include 'constants.include.F'
 double precision :: r12,contrib
 cst = -0.25d0
 r12  = (r1(1) - r2(1))*(r1(1) - r2(1))
 r12 += (r1(2) - r2(2))*(r1(2) - r2(2))
 r12 += (r1(3) - r2(3))*(r1(3) - r2(3))
 r12 = dsqrt(r12)
 contrib = ( 1.d0 - derf(mu * r12)) 
 erf_mu_sq_general = cst * contrib * contrib 

end


double precision function nabla_r12_bis(r1,r2,mu,grad_mu)
 implicit none 
 BEGIN_DOC
! - 1/(4 sqrt(pi) (\mu(r1))^2) \nabla_1 \mu(r1) . r_12  e^{(-\mu(r1)r12)^2} (1 - erf(\mu(r1)r12))/r12
 END_DOC
 double precision, intent(in) :: r1(3), r2(3), mu, grad_mu(3)
 double precision :: cst
 include 'constants.include.F'
 double precision :: r12,r12_sq,r12_vec(3)
 cst = -0.25 * inv_sq_pi

 nabla_r12_bis = 0.d0

 r12_vec(1) = (r1(1) - r2(1))
 r12_vec(2) = (r1(2) - r2(2))
 r12_vec(3) = (r1(3) - r2(3))
 r12_sq = r12_vec(1)*r12_vec(1) + r12_vec(2)*r12_vec(2) + r12_vec(3)*r12_vec(3)
 r12 = dsqrt(r12_sq)
 if(r12.lt.1.d-10)return
 nabla_r12_bis = cst / (mu*mu) * dexp(- mu * mu * r12_sq) & 
             * ( grad_mu(1) * r12_vec(1) & 
               + grad_mu(2) * r12_vec(2) &
               + grad_mu(3) * r12_vec(3) ) & 
             * (1.d0 - derf(mu * r12))/r12
end
