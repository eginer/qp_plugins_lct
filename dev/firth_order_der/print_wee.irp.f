program print_wee
 implicit none
 read_wf = .true.
 touch read_wf
 call write_wee_aa
!call routine_ab
end

subroutine write_wee_aa
 implicit none
 include 'utils/constants.include.F'
 integer :: i, npoints
 double precision :: accu_num,accu_ana_1,accu_ana_2
 double precision :: r1(3),r2(3)
 double precision :: r12,accu_n2_aa,accu_f_aa
 double precision :: f_hf_alpha_alpha,n2_hf_alpha_alpha 
 double precision :: r12max,delta,n2_deriv2,n2_deriv4,f_deriv2,f_deriv4,wee_paper,f_paper,n2_paper
 double precision :: mu_num,mu_ana,spherical_av_wee_paper,spherical_av_f,spherical_av_n2
 print*,'\\\\\\\\\\\\\\\\\'
 print*,' '
 r1(1) = x_center_wee
 r1(2) = y_center_wee
 r1(3) = z_center_wee
 r12=0.d0

 ! computation of mu with the ratio of spherical averaged quantities 
 call give_f_paper_alpha_alpha_hf_at_r1_r12(r1,delta_r12,accu_ana_1,f_deriv2,f_deriv4)
 call give_n2_alpha_alpha_hf_at_r1_r12(r1,delta_r12,accu_ana_2,n2_deriv2,n2_deriv4)
 wee_paper = (accu_ana_1)/(accu_ana_2)
 call give_mu_r12(wee_paper,delta_r12,mu_num)
 print*,'mu_num = ',mu_num
 ! computation of mu with the spherical average of the ratio 
 wee_paper = spherical_av_wee_paper(r1,delta_r12)
 call give_mu_r12(wee_paper,delta_r12,mu_ana)
 print*,'mu_ana = ',mu_ana

 r12max = 10.d0
 npoints = 1000
 delta = r12max / dble(npoints)
 do i = 1,npoints
  r12 += delta
  accu_num = 0.d0 
  accu_n2_aa = 0.d0 
  accu_f_aa = 0.d0 
  wee_paper = spherical_av_wee_paper(r1,r12)
  f_paper   = spherical_av_f(r1,r12)
  n2_paper   = spherical_av_n2(r1,r12)
  write(33,'(100(F16.10,X))')r12,wee_paper,erf(mu_num * r12)/r12,erf(mu_ana * r12)/r12,f_paper/n2_paper,f_paper,n2_paper
 
 enddo
 
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



