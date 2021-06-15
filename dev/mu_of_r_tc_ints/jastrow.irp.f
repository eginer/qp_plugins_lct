double precision function jastrow_mu(r1,r2,mu)
 implicit none
 double precision, intent(in) :: r1(3), r2(3), mu
 include 'constants.include.F'
 double precision :: r12
 integer :: m
 r12 = 0.d0
 do m = 1, 3
  r12 += (r1(m) - r2(m)) * (r1(m) - r2(m))
 enddo
 r12 = dsqrt(r12)
 jastrow_mu = 0.5d0 * r12 * (1.d0 - derf(mu*r12)) - 0.5d0 * inv_sq_pi / mu * dexp(-mu*mu*r12*r12)

end


subroutine grad_r1_jastrow_mu(r1,r2,mu,grad_mu,jastrow,grad_jastrow)
 implicit none
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
 jastrow= 0.5d0 * r12 * erfc_mu - 0.5d0 * inv_sq_pi * gauss 
 grad_jastrow = 0.d0
 if(r12.lt.1.d-6)return
 do m = 1, 3
  grad_jastrow(m) = 0.5d0 * r12_vec(m)/r12 * erfc_mu + 0.5d0 * inv_sq_pi /(mu*mu) * gauss * grad_mu(m)
 enddo
end

subroutine lapl_r1_jastrow_bis(r1,r2,mu_min,dx,lapl)
 implicit none
 double precision, intent(in) :: r1(3),r2(3),mu_min,dx
 double precision, intent(out):: lapl(3)
 double precision :: r(3), grad_mu(3), jastrow, grad_jastrow_p(3), grad_jastrow_m(3)
 double precision :: mu,rho_a_hf,rho_b_hf,mu_lda_damped, grad_mu_r1(3)
 integer :: m
 lapl = 0.d0
 call get_grad_damped_mu_lda(r1,dx,mu_min,grad_mu_r1)
 do m = 1, 3
  r = r1
  r(m) += dx
  if(constant_mu)then
   grad_mu = 0.d0
   mu = mu_erf 
  else
   call get_grad_damped_mu_lda(r,dx,mu_min,grad_mu)
!   grad_mu = grad_mu_r1
   call dm_dft_alpha_beta_at_r(r,rho_a_hf,rho_b_hf)
   mu = mu_lda_damped(rho_a_hf,rho_b_hf,mu_min)
  endif
  call grad_r1_jastrow_mu(r,r2,mu,grad_mu,jastrow,grad_jastrow_p)

  r = r1
  r(m) -= dx
  if(constant_mu)then
   grad_mu = 0.d0
   mu = mu_erf 
  else
   call get_grad_damped_mu_lda(r,dx,mu_min,grad_mu)
!   grad_mu = grad_mu_r1
   call dm_dft_alpha_beta_at_r(r,rho_a_hf,rho_b_hf)
   mu = mu_lda_damped(rho_a_hf,rho_b_hf,mu_min)
  endif
  call grad_r1_jastrow_mu(r,r2,mu,grad_mu,jastrow,grad_jastrow_m)
  lapl(m) = (grad_jastrow_p(m) - grad_jastrow_m(m)) / (2.d0 * dx)
 enddo
end

subroutine lapl_r2_jastrow_bis(r1,r2,mu_min,dx,lapl)
 implicit none
 double precision, intent(in) :: r1(3),r2(3),mu_min,dx
 double precision, intent(out):: lapl(3)
 double precision :: r(3), grad_mu(3), jastrow, grad_jastrow_p(3), grad_jastrow_m(3)
 double precision :: mu,rho_a_hf,rho_b_hf,mu_lda_damped
 integer :: m
 lapl = 0.d0
 if(constant_mu)then
  mu = mu_erf
  grad_mu = 0.d0
 else
  call get_grad_damped_mu_lda(r1,dx,mu_min,grad_mu)
  call dm_dft_alpha_beta_at_r(r1,rho_a_hf,rho_b_hf)
  mu = mu_lda_damped(rho_a_hf,rho_b_hf,mu_min)
 endif
 do m = 1, 3
  r = r2
  r(m) += dx
  call grad_r2_jastrow_mu(r1,r,mu,jastrow,grad_jastrow_p)

  r = r2
  r(m) -= dx
  call grad_r2_jastrow_mu(r1,r,mu,jastrow,grad_jastrow_m)
  lapl(m) = (grad_jastrow_p(m) - grad_jastrow_m(m)) / (2.d0 * dx)
 enddo
end

subroutine lapl_r1_jastrow(r1,r2,mu_min,dx,lapl)
 implicit none
 double precision, intent(in) :: r1(3),r2(3),mu_min,dx
 double precision, intent(out):: lapl(3)
 double precision :: r(3), jastrow, jastrow_plus, jastrow_minus,jastrow_mu,rho_a_hf,rho_b_hf,mu_lda_damped ,mu
 integer :: m
 do m = 1, 3 
  r = r1 
  call dm_dft_alpha_beta_at_r(r,rho_a_hf,rho_b_hf)
  if(constant_mu)then
   mu = mu_erf
  else
   mu = mu_lda_damped(rho_a_hf,rho_b_hf,mu_min)
  endif
  jastrow = jastrow_mu(r,r2,mu)

  r = r1
  r(m) += 2.d0 * dx
  call dm_dft_alpha_beta_at_r(r,rho_a_hf,rho_b_hf)
  if(constant_mu)then
   mu = mu_erf
  else
   mu = mu_lda_damped(rho_a_hf,rho_b_hf,mu_min)
  endif
  jastrow_plus = jastrow_mu(r,r2,mu)

  r = r1
  r(m) -= 2.d0 * dx
  call dm_dft_alpha_beta_at_r(r,rho_a_hf,rho_b_hf)
  if(constant_mu)then
   mu = mu_erf
  else
   mu = mu_lda_damped(rho_a_hf,rho_b_hf,mu_min)
  endif
  jastrow_minus= jastrow_mu(r,r2,mu)
  
  lapl(m) = (jastrow_plus  + jastrow_minus  - 2.d0 * jastrow)/(4.d0 * dx * dx)
 enddo
end


subroutine lapl_r2_jastrow(r1,r2,mu_min,dx,lapl)
 implicit none
 double precision, intent(in) :: r1(3),r2(3),mu_min,dx
 double precision, intent(out):: lapl(3)
 double precision :: r(3), jastrow, jastrow_plus, jastrow_minus,jastrow_mu,rho_a_hf,rho_b_hf,mu_lda_damped ,mu
 integer :: m
 do m = 1, 3 
  r = r2 
  call dm_dft_alpha_beta_at_r(r1,rho_a_hf,rho_b_hf)
  if(constant_mu)then
   mu = mu_erf
  else
   mu = mu_lda_damped(rho_a_hf,rho_b_hf,mu_min)
  endif
  jastrow = jastrow_mu(r1,r,mu)

  r = r2
  r(m) += 2.d0 * dx
  jastrow_plus = jastrow_mu(r1,r,mu)

  r = r2
  r(m) -= 2.d0 * dx
  jastrow_minus= jastrow_mu(r1,r,mu)
  
  lapl(m) = (jastrow_plus  + jastrow_minus  - 2.d0 * jastrow)/(4.d0 * dx * dx)
 enddo
end


subroutine grad_r2_jastrow_mu(r1,r2,mu,jastrow,grad_jastrow)
 implicit none
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
 jastrow = 0.5d0 * r12 * erfc_mu - 0.5d0 * inv_sq_pi * gauss 
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
 double precision :: mu_lda_damped,rho_a_hf,rho_b_hf, mu, mu_plus, mu_minus, grad_mu(3),dx
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
     mu_plus = mu_lda_damped(rho_a_hf,rho_b_hf,mu_min)
     r1(:) = final_grid_points(:,ipoint) 
     r1(m) -= dx 
     call dm_dft_alpha_beta_at_r(r1,rho_a_hf,rho_b_hf)
     mu_minus = mu_lda_damped(rho_a_hf,rho_b_hf,mu_min)
     grad_mu(m) = (mu_plus - mu_minus)/(2.d0 * dx)
    enddo
    r1(:) = final_grid_points(:,ipoint) 
    call dm_dft_alpha_beta_at_r(r1,rho_a_hf,rho_b_hf)
    mu = mu_lda_damped(rho_a_hf,rho_b_hf,mu_min)
    call grad_r1_jastrow_mu(r1,r2,mu,grad_mu,jastrow,grad_jastrow)
    do m = 1, 3
     r1(:) = final_grid_points(:,ipoint) 
     r1(m) += dx 
     call dm_dft_alpha_beta_at_r(r1,rho_a_hf,rho_b_hf)
     mu_plus = mu_lda_damped(rho_a_hf,rho_b_hf,mu_min)
     jastrow_plus = jastrow_mu(r1,r2,mu_plus)
 
     r1(:) = final_grid_points(:,ipoint) 
     r1(m) -= dx 
     call dm_dft_alpha_beta_at_r(r1,rho_a_hf,rho_b_hf)
     mu_minus = mu_lda_damped(rho_a_hf,rho_b_hf,mu_min)
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


subroutine test_nabla_1_sq_ao
 implicit none
 integer :: i,j,k,l
 double precision :: pure_num, comp_num(3),delta,tot, analy
 provide scalar_mu_r_pot_chemist_ao 
 do i = 1, ao_num
  do j = i,i
   do k = i,i
    do l = j,j
     call get_num_ints_j_sq(i,k,j,l,pure_num, comp_num, delta)
     tot = comp_num(1) + comp_num(2) + comp_num(3)
     analy = scalar_mu_r_pot_chemist_ao(i,k,j,l)
     print*,'i,k,j,l',i,k,j,l
     print*,pure_num, tot, delta
     print*,analy,dabs(pure_num - analy), dabs(tot - analy)
    enddo
   enddo
  enddo
 enddo
end

subroutine get_num_ints_j_sq(i_ao,k_ao,j_ao,l_ao,pure_num, comp_num, delta)
 implicit none
 BEGIN_DOC
! you enter with (i,k|j,l) in the AO basis for chemist notation 
!
! you get out with pure_num = \int dr1 dr2 AO_i(r1) AO_k(r1) (-1/2 (\nabla_1 j(r1,r2))^2 + (-1/2 (\nabla_2 j(r1,r2))^2  AO_j(r2) AO_l(r2)
! 
! and the various components of (-1/2 (\nabla_1 j(r1,r2))^2 + (-1/2 (\nabla_2 j(r1,r2))^2 
!
! and delta = difference between both
 END_DOC
 include 'constants.include.F'
 integer, intent(in) :: i_ao,k_ao,j_ao,l_ao
 double precision, intent(out):: pure_num, comp_num(3),delta
 integer :: ipoint,i,j,k,l,m,jpoint
 double precision :: r1(3),r2(3),weight1,weight2,weight_tot,ao_prod_r2,ao_prod_r1,mu_min,r12
 double precision :: mu_lda_damped,rho_a_hf,rho_b_hf, mu, mu_plus, mu_minus, grad_mu(3),dx
 double precision :: jastrow,grad_jastrow(3),jastrow_plus,jastrow_minus, grad_jastrow_sq_anal
 double precision :: accu,jastrow_mu, grad_r1_jastrow_sq
 double precision :: nabla_r12_bis,nabla_sq_term_general,erf_mu_sq_general
 double precision :: contrib(3)
 do l = 7, 7
  dx = 10.d0**(-l)
  mu_min = mu_erf
  pure_num = 0.d0
  comp_num = 0.d0
  delta    = 0.d0
  do jpoint = 1, n_points_final_grid ! r2
   weight2 = final_weight_at_r_vector(jpoint)
   r2(:) = final_grid_points(:,jpoint)
   ao_prod_r2 = aos_in_r_array_transp(jpoint,j_ao) * aos_in_r_array_transp(jpoint,l_ao) * weight2
   do ipoint = 1, n_points_final_grid ! r1
     weight1 = final_weight_at_r_vector(ipoint)
     ao_prod_r1 = aos_in_r_array_transp(ipoint,i_ao) * aos_in_r_array_transp(ipoint,k_ao) * weight1
     r1(:) = final_grid_points(:,ipoint) 
     r12 = 0.d0
     do m = 1, 3 ! compute grad mu
      r12 += (r1(m) - r2(m))*(r1(m) - r2(m))
     enddo
     r12 = dsqrt(r12)
     if(r12.lt.1.d-8)cycle
     grad_mu(:) = grad_mu_of_r_for_ints(:,ipoint,1)
     mu = mu_of_r_for_ints(ipoint,1)
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
     pure_num += grad_r1_jastrow_sq                         * ao_prod_r1 * ao_prod_r2
     contrib(1) = erf_mu_sq_general(r1,r2,mu)               
     contrib(2) = nabla_sq_term_general(r1,r2,mu,grad_mu)   
     contrib(3) = nabla_r12_bis(r1,r2,mu,grad_mu)           
     comp_num(1) += contrib(1) * ao_prod_r1 * ao_prod_r2 
     comp_num(2) += contrib(2) * ao_prod_r1 * ao_prod_r2 
     comp_num(3) += contrib(3) * ao_prod_r1 * ao_prod_r2  

     grad_jastrow_sq_anal = contrib(1) + contrib(2) + contrib(3) 
     delta    += dabs(grad_r1_jastrow_sq - grad_jastrow_sq_anal) * ao_prod_r1 * ao_prod_r2
   enddo
  enddo
 enddo
end

subroutine test_non_hermit_ao
 implicit none
 integer :: i,j,k,l
 double precision :: pure_num, comp_num(2),delta,tot, analy
 double precision :: accu
 provide scalar_mu_r_pot_chemist_ao 
 accu = 0.d0
 do i = 1, ao_num
  do j = 1, ao_num
   do k = 1, ao_num
    do l = 1, ao_num
     call get_num_ints_j_non_hermit(i,k,j,l,pure_num, comp_num, delta)
     tot = comp_num(1) + comp_num(2) 
     analy = deriv_mu_r_pot_chemist_ao(i,k,j,l)
     print*,'i,k,j,l',i,k,j,l
     print*,pure_num, tot, delta
     print*,analy,dabs(pure_num - analy), dabs(tot - analy)
     accu += dabs(tot - analy)
    enddo
   enddo
  enddo
 enddo
 print*,'accu = ',accu/dble(ao_num**4)
 
end

subroutine get_num_ints_j_non_hermit(i_ao,k_ao,j_ao,l_ao,pure_num, comp_num, delta)
 implicit none
 BEGIN_DOC
! you enter with (i,k|j,l) in the AO basis for chemist notation 
!
! you get out with pure_num = \int dr1 dr2 AO_l(r2) AO_k(r1) (\nabla_1 u(r1,r2) + \naba_2 u(r1,r2))  AO_j(r2) AO_i(r1) 
! 
! and the various components of (\nabla_1 u(r1,r2) + \naba_2 u(r1,r2))
!
! and delta = difference between both
 END_DOC
 include 'constants.include.F'
 integer, intent(in) :: i_ao,k_ao,j_ao,l_ao
 double precision, intent(out):: pure_num, comp_num(2),delta
 integer :: ipoint,i,j,k,l,m,jpoint
 double precision :: r1(3),r2(3),weight1,weight2,weight_tot,mu_min,r12
 double precision :: ao_prod_r2,ao_prod_r1,ao_grad_r1(3), ao_grad_r2(3)
 double precision :: mu_lda_damped,rho_a_hf,rho_b_hf, mu, mu_plus, mu_minus, grad_mu(3),dx
 double precision :: jastrow,grad_jastrow_r1(3),grad_jastrow_r2(3),jastrow_plus,jastrow_minus
 double precision :: accu,jastrow_mu, grad_r1_jastrow_sq
 double precision :: nabla_r12_bis,nabla_sq_term_general,erf_mu_sq_general, non_hermit_anal,non_hermit_num
 double precision :: contrib(2),func_non_hermit_at_r1_bis,vec(3),func_non_hermit_at_r1
 do l = 7, 7
  dx = 10.d0**(-l)
  mu_min = mu_erf
  pure_num = 0.d0
  comp_num = 0.d0
  delta    = 0.d0
  pure_num = 0.d0
  do jpoint = 1, n_points_final_grid ! r2
   weight2 = final_weight_at_r_vector(jpoint)
   r2(:) = final_grid_points(:,jpoint)
   ao_prod_r2 = aos_in_r_array_transp(jpoint,j_ao) * aos_in_r_array_transp(jpoint,l_ao) * weight2
   do m = 1, 3
    ao_grad_r2(m) = aos_in_r_array_transp(jpoint,l_ao) * weight2 * aos_grad_in_r_array_transp_3(m,jpoint,j_ao)
   enddo
   do ipoint = 1, n_points_final_grid ! r1
     weight1 = final_weight_at_r_vector(ipoint)
     ao_prod_r1 = aos_in_r_array_transp(ipoint,i_ao) * aos_in_r_array_transp(ipoint,k_ao) * weight1
     do m = 1, 3
      ao_grad_r1(m) = aos_in_r_array_transp(ipoint,k_ao) * weight1 * aos_grad_in_r_array_transp_3(m,ipoint,i_ao)
     enddo
     r1(:) = final_grid_points(:,ipoint) 
     r12 = 0.d0
     do m = 1, 3 ! compute r12 
      r12 += (r1(m) - r2(m))*(r1(m) - r2(m))
     enddo
     r12 = dsqrt(r12)
     if(r12.lt.1.d-8)cycle
     grad_mu(:) = grad_mu_of_r_for_ints(:,ipoint,1)
     mu = mu_of_r_for_ints(ipoint,1)
!     mu = mu_erf 

     call grad_r1_jastrow_mu(r1,r2,mu,grad_mu,jastrow,grad_jastrow_r1)
     non_hermit_num = 0.d0
     do m = 1, 3
      non_hermit_num -= grad_jastrow_r1(m) * ao_prod_r2 * ao_grad_r1(m)
     enddo
     call grad_r2_jastrow_mu(r1,r2,mu,jastrow,grad_jastrow_r2)
     do m = 1, 3
      non_hermit_num -= grad_jastrow_r2(m) * ao_prod_r1 * ao_grad_r2(m)
     enddo
!     contrib(1) = func_non_hermit_at_r1_bis(ipoint,jpoint,i_ao,k_ao,j_ao,l_ao,mu) * weight1 * weight2
     contrib(1) = func_non_hermit_at_r1(ipoint,jpoint,i_ao,k_ao,j_ao,l_ao) * weight1 * weight2
     call gamma_at_r1(r1,r2,mu,grad_mu,vec)
     contrib(2) = 0.d0
     do m = 1, 3
      contrib(2) -= vec(m) * ao_prod_r2 * ao_grad_r1(m)
     enddo
     non_hermit_anal = contrib(1) + contrib(2)
     pure_num += non_hermit_num
     comp_num(1) += contrib(1) + contrib(2) 
     delta    += dabs(non_hermit_anal - non_hermit_num) 
   enddo
  enddo
 enddo
end

subroutine test_lapl_j_ao
 implicit none
 integer :: i,j,k,l
 double precision :: pure_num, comp_num(4),delta,tot, analy
 double precision :: accu,get_ao_two_e_integral_eff_pot,get_ao_two_e_integral_erf
 provide scalar_mu_r_pot_chemist_ao 
 accu = 0.d0
! do i = 1, ao_num
!  do j = 1, ao_num
!   do k = 1, ao_num
!    do l = 1, ao_num
 do i = 1, ao_num
  do j = 4, 4
   do k = j,j
    do l = i,i 
     call get_num_ints_lapl_j(i,k,j,l,pure_num, comp_num, delta)
     tot = comp_num(1) + comp_num(2)  + comp_num(3) + comp_num(4)
      analy = scalar_mu_r_pot_chemist_ao(i,k,j,l)
!     analy = get_ao_two_e_integral_erf(i,j,k,l,ao_integrals_erf_map)   & 
!           + get_ao_two_e_integral_eff_pot(i,j,k,l,ao_integrals_eff_pot_map) 
     if(dabs(pure_num).lt.1.d-10)cycle
!     if(dabs(pure_num - tot)/dabs(pure_num).lt.1.d-4)cycle
     print*,'i,k,j,l',i,k,j,l
     print*,pure_num, tot, dabs(pure_num - tot)
     print*,analy,dabs(tot - analy), dabs(tot - analy)/dabs(tot)
     accu += dabs(tot - pure_num)
!     accu += dabs(tot - analy)
    enddo
   enddo
  enddo
 enddo
 print*,'accu = ',accu/dble(ao_num**4)
 
end

subroutine get_num_ints_lapl_j(i_ao,k_ao,j_ao,l_ao,pure_num, comp_num, delta)
 implicit none
 BEGIN_DOC
! you enter with (i,k|j,l) in the AO basis for chemist notation 
!
! you get out with pure_num = \int dr1 dr2 AO_i(r1) AO_k(r1) (-1/2 (\nabla_1)^2 j(r1,r2) + (-1/2 (\nabla_2)^2 j(r1,r2)  AO_j(r2) AO_l(r2)
! 
! and the various components of (-1/2 (\nabla_1 j(r1,r2))^2 + (-1/2 (\nabla_2 j(r1,r2))^2 
!
! and delta = difference between both
 END_DOC
 include 'constants.include.F'
 integer, intent(in) :: i_ao,k_ao,j_ao,l_ao
 double precision, intent(out):: pure_num, comp_num(4),delta
 integer :: ipoint,i,j,k,l,m,jpoint
 double precision :: r1(3),r2(3),weight1,weight2,weight_tot,ao_prod_r2,ao_prod_r1,mu_min,r12
 double precision :: mu_lda_damped,rho_a_hf,rho_b_hf, mu, mu_plus, mu_minus, grad_mu(3),dx
 double precision :: jastrow,grad_jastrow(3),jastrow_plus,jastrow_minus, lapl_jastrow_sq_anal
 double precision :: accu,jastrow_mu, lapl_r1_jast, lapl_r2_jast,lapl_tot_jast
 double precision :: nabla_r12_bis,nabla_sq_term_general,erf_mu_sq_general,lapl_j(3)
 double precision :: contrib(4),gauss_mu_r12,derf_mu_x,gauss_grad_mu_r1,grad_gamma_r1,grad_gamma_r1_bis
 do l = 5, 5
  dx = 10.d0**(-l)
  mu_min = mu_erf
  pure_num = 0.d0
  comp_num = 0.d0
  delta    = 0.d0
  do jpoint = 1, n_points_final_grid ! r2
   weight2 = final_weight_at_r_vector(jpoint)
   r2(:) = final_grid_points(:,jpoint)
   ao_prod_r2 = aos_in_r_array_transp(jpoint,j_ao) * aos_in_r_array_transp(jpoint,l_ao) * weight2
   do ipoint = 1, n_points_final_grid ! r1
     weight1 = final_weight_at_r_vector(ipoint)
     ao_prod_r1 = aos_in_r_array_transp(ipoint,i_ao) * aos_in_r_array_transp(ipoint,k_ao) * weight1
     r1(:) = final_grid_points(:,ipoint) 
     r12 = 0.d0
     do m = 1, 3 ! compute grad mu
      r12 += (r1(m) - r2(m))*(r1(m) - r2(m))
     enddo
     r12 = dsqrt(r12)
     if(r12.lt.1.d-8)cycle
     grad_mu(:) = grad_mu_of_r_for_ints(:,ipoint,1)
     mu = mu_of_r_for_ints(ipoint,1)
!     call lapl_r1_jastrow(r1,r2,mu_min,dx,lapl_j)
     call lapl_r1_jastrow_bis(r1,r2,mu_min,dx,lapl_j)
     lapl_r1_jast = 0.d0
     do m = 1, 3
      lapl_r1_jast -= 0.5d0 * lapl_j(m)
     enddo
!     call lapl_r2_jastrow(r1,r2,mu_min,dx,lapl_j)
     call lapl_r2_jastrow_bis(r1,r2,mu_min,dx,lapl_j)
     lapl_r2_jast = 0.d0
     do m = 1, 3
      lapl_r2_jast -= 0.5d0 * lapl_j(m)
     enddo

     lapl_tot_jast = lapl_r1_jast + lapl_r2_jast +1.d0/r12 
     pure_num += lapl_tot_jast * ao_prod_r1 * ao_prod_r2
     contrib = 0.d0
     contrib(1) = derf_mu_x(mu,r12)
     contrib(2) = gauss_mu_r12(r12,mu)
     contrib(3) = gauss_grad_mu_r1(r1,r2,mu,grad_mu)
     contrib(4) = -0.5d0 * grad_gamma_r1_bis(r1,r2, dx,mu_min)
     comp_num(1) += contrib(1) * ao_prod_r1 * ao_prod_r2 
     comp_num(2) += contrib(2) * ao_prod_r1 * ao_prod_r2 
     comp_num(3) += contrib(3) * ao_prod_r1 * ao_prod_r2  
     comp_num(4) += contrib(4) * ao_prod_r1 * ao_prod_r2  

     lapl_jastrow_sq_anal = contrib(1) + contrib(2) + contrib(3) + contrib(4)
!     if(dabs(lapl_tot_jast - lapl_jastrow_sq_anal)/dabs(lapl_tot_jast).gt.1.d-3)then
!      print*,'r1,r2'
!      print*,r1
!      print*,r2
!      print*,ipoint,jpoint
!      print*,lapl_tot_jast,lapl_jastrow_sq_anal,dabs(lapl_tot_jast - lapl_jastrow_sq_anal)/dabs(lapl_tot_jast)
!     endif
     delta    += dabs(lapl_tot_jast - lapl_jastrow_sq_anal) * ao_prod_r1 * ao_prod_r2
   enddo
  enddo
 enddo
end

subroutine test_grad_gamma_r1
 implicit none
 integer :: i,j,k,l
 double precision :: pure_num, analy
 double precision :: accu, int_pp
 provide scalar_mu_r_pot_chemist_ao 
 accu = 0.d0
 do i = 1, ao_num
  do j = 1, ao_num
!   do k = 1, ao_num
   do k = j,j
!    do l = 1, ao_num
    do l = i,i
! do i = 1, 1
!  do j = 1, 1
!   do k = 1, 1
!    do l = 1, 2
     call get_num_ints_grad_gamma_r1(i,k,j,l,pure_num,int_pp)
     analy = scalar_mu_r_pot_chemist_ao(i,k,j,l)
     if(dabs(pure_num).lt.1.d-10)cycle
!     if(dabs(pure_num - analy )/dabs(pure_num).lt.1.d-4)cycle
     if(dabs(pure_num - int_pp)/dabs(pure_num).lt.1.d-4)cycle
     print*,'i,k,j,l',i,k,j,l
!     print*,pure_num, analy, dabs(pure_num - analy)
!     accu += dabs(pure_num - analy)
     print*,pure_num, int_pp, dabs(pure_num - int_pp)
     print*,analy, dabs(analy - pure_num)
     accu += dabs(pure_num - int_pp)
    enddo
   enddo
  enddo
 enddo
 print*,'accu = ',accu/dble(ao_num**4)
 
end

subroutine get_num_ints_grad_gamma_r1(i_ao,k_ao,j_ao,l_ao,pure_num,int_pp)
 implicit none
 BEGIN_DOC
! you enter with (i,k|j,l) in the AO basis for chemist notation 
!
! you get out with pure_num = \int dr1 dr2 AO_i(r1) AO_k(r1) -1/2 * grad_r1 . gamma(r1) AO_j(r2) AO_l(r2)
! 
 END_DOC
 include 'constants.include.F'
 integer, intent(in) :: i_ao,k_ao,j_ao,l_ao
 double precision, intent(out):: pure_num,int_pp
 integer :: ipoint,i,j,k,l,m,jpoint
 double precision :: r1(3),r2(3),weight1,weight2
 double precision :: weight_tot,ao_prod_r2,ao_prod_r1,mu_min,r12
 double precision :: grad_gamma_r1_bis,dx,cst_lapl_gamma,lapl_gamma
 cst_lapl_gamma = 0.25d0 * inv_sq_pi
 do l = 5, 5
  dx = 10.d0**(-l)
  mu_min = mu_erf
  pure_num = 0.d0
  int_pp = 0.d0
  do jpoint = 1, n_points_final_grid ! r2
   weight2 = final_weight_at_r_vector(jpoint)
   r2(:) = final_grid_points(:,jpoint)
   ao_prod_r2 = aos_in_r_array_transp(jpoint,j_ao) * aos_in_r_array_transp(jpoint,l_ao) * weight2
   do ipoint = 1, n_points_final_grid ! r1
     weight1 = final_weight_at_r_vector(ipoint)
     ao_prod_r1 = aos_in_r_array_transp(ipoint,i_ao) * aos_in_r_array_transp(ipoint,k_ao) * weight1
     r1(:) = final_grid_points(:,ipoint) 
     r12 = 0.d0
     do m = 1, 3
      r12 += (r1(m) - r2(m))*(r1(m) - r2(m))
     enddo
     r12 = dsqrt(r12)
     if(r12.lt.1.d-8)cycle
     pure_num += -0.5d0 * grad_gamma_r1_bis(r1,r2, dx,mu_min) * ao_prod_r1 * ao_prod_r2
     int_pp += lapl_gamma(ipoint,jpoint,i_ao,k_ao,cst_lapl_gamma) * ao_prod_r2 * weight1 
   enddo
  enddo
 enddo
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
 do m = 1, 3
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

subroutine gamma_at_r1(r1,r2,mu,grad_mu,vec)
 implicit none
 double precision, intent(in) :: r1(3), r2(3), mu, grad_mu(3)
 double precision, intent(out):: vec(3)
 double precision :: contrib
 double precision :: cst,r12
 integer :: m
 r12 = 0.d0
 do m = 1, 3
  r12 += (r1(m) - r2(m))*(r1(m) - r2(m))
 enddo
 include 'constants.include.F'
 cst = 0.5d0 * inv_sq_pi /(mu*mu) * dexp(-mu*mu*r12)
 vec = grad_mu * cst 
end

double precision function gamma_scal_at_r1(r1,r2,mu)
 implicit none
 double precision, intent(in) :: r1(3), r2(3), mu
 double precision :: r12
 integer :: m
 r12 = 0.d0
 do m = 1, 3
  r12 += (r1(m) - r2(m))*(r1(m) - r2(m))
 enddo
 include 'constants.include.F'
 gamma_scal_at_r1 = 0.5d0 * inv_sq_pi /(mu*mu) * dexp(-mu*mu*r12)
end

double precision function gauss_grad_mu_r1(r1,r2,mu,grad_mu)
 implicit none
 double precision, intent(in) :: r1(3), r2(3), mu, grad_mu(3)
 include 'constants.include.F'
 double precision :: cst,r12,r12_vec(3)
 integer :: m
 r12 = 0.d0
 do m = 1, 3
  r12_vec(m) = r1(m) - r2(m)
  r12 += r12_vec(m) * r12_vec(m) 
 enddo
 cst = 0.5d0 * inv_sq_pi * dexp(-mu*mu*r12)
 gauss_grad_mu_r1 = 0.d0
 do m = 1, 3
  gauss_grad_mu_r1 += r12_vec(m) * grad_mu(m)
 enddo
 gauss_grad_mu_r1 *= cst 
end


double precision function grad_gamma_r1(r1,r2,mu, grad_mu, dx,mu_min)
 implicit none
 double precision, intent(in) :: r1(3), r2(3), mu, grad_mu(3), dx, mu_min
 double precision :: r(3), gamma_scal_at_r1, mu_lda_damped, rho_a_hf,rho_b_hf
 double precision :: deriv_f(3), lapl_mu(3), f_p, f_m, mu_p, mu_m,f
 integer :: m
 f = gamma_scal_at_r1(r1,r2,mu)
 do m = 1, 3
  r = r1
  r(m) += dx
  call dm_dft_alpha_beta_at_r(r,rho_a_hf,rho_b_hf)
  mu_p = mu_lda_damped(rho_a_hf,rho_b_hf,mu_min)
  f_p = gamma_scal_at_r1(r,r2,mu_p)

  r = r1
  r(m) -= dx
  call dm_dft_alpha_beta_at_r(r,rho_a_hf,rho_b_hf)
  mu_m = mu_lda_damped(rho_a_hf,rho_b_hf,mu_min)
  f_m = gamma_scal_at_r1(r,r2,mu_m)

  deriv_f(m) = (f_p - f_m)/(2.d0 *dx)
 enddo
 call get_lapl_mu_lda(r1,dx,mu_min,lapl_mu)
 grad_gamma_r1 = 0.d0
 do m = 1, 3
  grad_gamma_r1 += deriv_f(m) * grad_mu(m) + f * lapl_mu(m)
 enddo
end

double precision function gamma_r1_comp(r1,r2,mu,mu_min,dx,m)
 implicit none
 double precision, intent(in) :: r1(3), r2(3), mu, mu_min,dx
 integer, intent(in) :: m
 double precision :: grad_mu_r1(3),gamma_scal_at_r1
 call get_grad_damped_mu_lda(r1,dx,mu_min,grad_mu_r1)
 gamma_r1_comp = gamma_scal_at_r1(r1,r2,mu) * grad_mu_r1(m)

end

double precision function grad_gamma_r1_bis(r1,r2, dx,mu_min)
 implicit none
 double precision, intent(in) :: r1(3), r2(3), dx, mu_min
 integer :: m
 double precision :: r(3),gamma_r1_comp, gamma_m, gamma_p,mu,rho_a_hf,rho_b_hf,mu_lda_damped
 grad_gamma_r1_bis = 0.d0
 do m = 1, 3
  r = r1
  r(m) += dx
  call dm_dft_alpha_beta_at_r(r,rho_a_hf,rho_b_hf)
  mu = mu_lda_damped(rho_a_hf,rho_b_hf,mu_min)
  gamma_p = gamma_r1_comp(r,r2,mu,mu_min,dx,m)

  r = r1
  r(m) -= dx
  call dm_dft_alpha_beta_at_r(r,rho_a_hf,rho_b_hf)
  mu = mu_lda_damped(rho_a_hf,rho_b_hf,mu_min)
  gamma_m = gamma_r1_comp(r,r2,mu,mu_min,dx,m)

  grad_gamma_r1_bis += (gamma_p - gamma_m)/(2.d0 * dx)
 enddo
end

double precision function gauss_mu_r12(r12,mu)
 implicit none
 double precision, intent(in) :: r12, mu
 include 'constants.include.F'
 gauss_mu_r12 = mu * inv_sq_pi * dexp(-mu*mu*r12*r12)
end

