program print_wee
 implicit none
 read_wf = .true.
 touch read_wf
!call routine
 call routine_ab
end

subroutine routine
 implicit none
 integer :: i,j,nr
 double precision :: r12, dr12,r12max
 double precision :: spherical_av_n2,n2,mu
 double precision :: spherical_av_f_paper,f_paper
 double precision :: r(3),local_potential,two_bod
 nr = 1000 
 r12max = 2.d0
 dr12 = r12max/dble(nr)
 r = 0.d0
 r(1) = 0.d0
 r(2) = 0.d0
 r(3) = 0.5d0
!!!!!!!!! definition of mu
 r12 = 0.001d0
 call give_eff_inter_alpha_alpha_hf_at_r1_r12(r,r12,f_paper,n2)
 f_paper = f_paper/n2
 call give_mu_r12(f_paper,r12,mu)
 print*,'muuuuuuuuuuuuuuuuuuuu = ',mu


 r12 = 0.d0
 do i = 1, nr
  r12 += dr12 
  n2 = spherical_av_n2(r,r12)
  f_paper = spherical_av_f_paper(r,r12)
  write(33,'(10(F16.10,X))')r12,f_paper/n2,erf(mu * r12)/r12,f_paper,n2
 enddo
end



double precision  function spherical_av_n2(r1,r12)
 implicit none
 include 'utils/constants.include.F'
 double precision, intent(in) :: r1(3),r12
 integer :: k
 double precision :: n2_hf_alpha_alpha,r2(3)
 spherical_av_n2 = 0.d0
 do k = 1,n_points_integration_angular
  r2(1)= r1(1)+r12*angular_quadrature_points(k,1)
  r2(2)= r1(2)+r12*angular_quadrature_points(k,2)
  r2(3)= r1(3)+r12*angular_quadrature_points(k,3)
  call give_n2_alpha_alpha_hf_at_r1_r2(r1,r2,n2_hf_alpha_alpha) 
  spherical_av_n2 += n2_hf_alpha_alpha*weights_angular_points(k)
 enddo
 spherical_av_n2 = spherical_av_n2 * 0.25d0 / pi 
end

double precision  function spherical_av_f_paper(r1,r12)
 implicit none
 include 'utils/constants.include.F'
 double precision, intent(in) :: r1(3),r12
 integer :: k
 double precision :: f_hf_alpha_alpha,r2(3)
 spherical_av_f_paper = 0.d0
 do k = 1,n_points_integration_angular
  r2(1)= r1(1)+r12*angular_quadrature_points(k,1)
  r2(2)= r1(2)+r12*angular_quadrature_points(k,2)
  r2(3)= r1(3)+r12*angular_quadrature_points(k,3)
  call give_f_alpha_alpha_hf_at_r1_r2(r1,r2,f_hf_alpha_alpha) 
  spherical_av_f_paper += f_hf_alpha_alpha * weights_angular_points(k)
 enddo
 spherical_av_f_paper = spherical_av_f_paper * 0.25d0 / pi 
end



subroutine routine_ab
 implicit none
 integer :: i,j,nr
 double precision :: r12, dr12,r12max
 double precision :: spherical_av_n2_ab,n2,mu
 double precision :: spherical_av_f_paper_ab,f_paper
 double precision :: r(3),local_potential,two_bod
 nr = 100 
 r12max = 2.d0
 dr12 = r12max/dble(nr)
 r = 0.d0
 r(1) = 0.5d0
 r(2) = 1.5d0
 r(3) = 1.5d0
!!!!!!!!! definition of mu
 r12 = 0.01d0
!r12 = 0.d0

!double precision :: accu_ana_1,f2_0,f_deriv2,f_deriv4,accu_ana_2,n2_0,n2_deriv2,n2_deriv4
!call give_f_paper_alpha_beta_hf_at_r1_r12(r,r12,accu_ana_1,f2_0,f_deriv2,f_deriv4)
!call give_n2_alpha_beta_hf_at_r1_r12(r,r12,accu_ana_2,n2_0,n2_deriv2,n2_deriv4)
!print*,'w paper 1', accu_ana_1/accu_ana_2


 call give_eff_inter_alpha_beta_hf_at_r1_r12(r,r12,f_paper,n2)
 f_paper = f_paper/n2
 print*,'w paper',f_paper
 call give_mu_r12(f_paper,r12,mu)
 print*,'muuuuuuuuuuuuuuuuuuuu = ',mu


 r12 = 0.d0
 do i = 1, nr
  r12 += dr12 
  n2 = spherical_av_n2_ab(r,r12)
  f_paper = spherical_av_f_paper_ab(r,r12)
  write(33,'(10(F16.10,X))')r12,f_paper/n2,erf(mu * r12)/r12,f_paper,n2
 enddo
end



double precision  function spherical_av_n2_ab(r1,r12)
 implicit none
 include 'utils/constants.include.F'
 double precision, intent(in) :: r1(3),r12
 integer :: k
 double precision :: n2_hf_alpha_beta,r2(3)
 spherical_av_n2_ab = 0.d0
 do k = 1,n_points_integration_angular
  r2(1)= r1(1)+r12*angular_quadrature_points(k,1)
  r2(2)= r1(2)+r12*angular_quadrature_points(k,2)
  r2(3)= r1(3)+r12*angular_quadrature_points(k,3)
  call give_n2_alpha_beta_hf_at_r1_r2(r1,r2,n2_hf_alpha_beta) 
  spherical_av_n2_ab += n2_hf_alpha_beta*weights_angular_points(k)
 enddo
 spherical_av_n2_ab = spherical_av_n2_ab * 0.25d0 / pi 
end

double precision  function spherical_av_f_paper_ab(r1,r12)
 implicit none
 include 'utils/constants.include.F'
 double precision, intent(in) :: r1(3),r12
 integer :: k
 double precision :: f_hf_alpha_beta,r2(3)
 spherical_av_f_paper_ab = 0.d0
 do k = 1,n_points_integration_angular
  r2(1)= r1(1)+r12*angular_quadrature_points(k,1)
  r2(2)= r1(2)+r12*angular_quadrature_points(k,2)
  r2(3)= r1(3)+r12*angular_quadrature_points(k,3)
  call give_f_alpha_beta_hf_at_r1_r2(r1,r2,f_hf_alpha_beta) 
  spherical_av_f_paper_ab += f_hf_alpha_beta * weights_angular_points(k)
 enddo
 spherical_av_f_paper_ab = spherical_av_f_paper_ab * 0.25d0 / pi 
end
