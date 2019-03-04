program print_wee
 implicit none
 read_wf = .true.
 touch read_wf
 call write_mu_r
!call write_wee_aa
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
 double precision :: spherical_av_f_manu

 character*(128) :: output
 integer :: i_unit_output,getUnitAndOpen
 provide ezfio_filename
 character*(128) :: filename
 write (filename, "(F3.1,A1,F3.1,A1,F3.1)")x_center_wee,'_',y_center_wee,'_',z_center_wee

 output=trim(ezfio_filename)//'.'//trim(filename)
 output=trim(output)
 print*,'output = ',trim(output)
 i_unit_output = getUnitAndOpen(output,'w')
 
 
 f_hf_alpha_alpha = spherical_av_f(r1,delta_r12)
 print*,'\\\\\\\\\\\\\\\\\'
 print*,' '
 r1(1) = x_center_wee
 r1(2) = y_center_wee
 r1(3) = z_center_wee

 ! computation of mu with the ratio of spherical averaged quantities 
 call give_f_paper_alpha_alpha_hf_at_r1_r12(r1,delta_r12,accu_ana_1,f_deriv2,f_deriv4)
 call give_n2_alpha_alpha_hf_at_r1_r12(r1,delta_r12,accu_ana_2,n2_deriv2,n2_deriv4)
 wee_paper = (accu_ana_1)/(accu_ana_2)
 print*,' wee_paper ',wee_paper
 print*,'F_analytic, n2_analytic'
 print*,accu_ana_1,accu_ana_2
 f_hf_alpha_alpha = spherical_av_f(r1,delta_r12)
 n2_hf_alpha_alpha= spherical_av_n2(r1,delta_r12)
 print*,'F_numeric , n2_numeric '
 print*,f_hf_alpha_alpha,n2_hf_alpha_alpha
!call give_mu_r12(wee_paper,delta_r12,mu_num)
!print*,'mu_num = ',mu_num
 ! computation of mu with the spherical average of the ratio 
 wee_paper = spherical_av_wee_paper(r1,delta_r12)
 print*,' wee_paper ',wee_paper
!call give_mu_r12(wee_paper,delta_r12,mu_ana)
 print*,'mu_ana = ',mu_ana
!stop
 r12max = 5.d0
 npoints = 100
 delta = r12max / dble(npoints)
 do i = 1,npoints
  r12 += delta
  accu_num = 0.d0 
  accu_n2_aa = 0.d0 
  accu_f_aa = 0.d0 
  wee_paper = spherical_av_wee_paper(r1,r12)
  f_paper   = spherical_av_f(r1,r12)
  n2_paper   = spherical_av_n2(r1,r12)
 !write(33,'(100(F16.10,X))')r12,wee_paper,erf(mu_num * r12)/r12,erf(mu_ana * r12)/r12,f_paper/n2_paper,f_paper,n2_paper
  write(i_unit_output,'(100(F16.10,X))')r12,wee_paper,f_paper/n2_paper,f_paper,n2_paper
 
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


subroutine write_mu_r
 implicit none
 integer :: n_points,i_point
 double precision :: dr,rmin,rmax,r(3),local_potential,two_bod
 double precision, allocatable :: mu_average_z(:),weight_average_z(:),z_tab(:),density_z(:,:),mu_average_z_dens(:)
 character*(128) :: output
 integer :: i_unit_output,getUnitAndOpen
 double precision :: dm_a,dm_b,ec,mu
 provide ezfio_filename
 character*(128) :: filename
 write (filename, "(A3)")"wee"
 output=trim(ezfio_filename)//'.'//trim(filename)
 print*,'output = ',trim(output)
 i_unit_output = getUnitAndOpen(output,'w')
 n_points = 10000
 rmin = nucl_coord(1,3)-3.d0
 
 rmax = nucl_coord(1,3)+6.d0
 dr   = (rmax - rmin)/dble(n_points) 
 r = 0.d0


 allocate(mu_average_z_dens(n_points+1),mu_average_z(n_points+1),weight_average_z(n_points+1),z_tab(n_points+1),density_z(2,n_points+1))
 z_tab = 0.d0
 r = 0.d0
 r(3) = rmin
 do i_point = 1, n_points
  z_tab(i_point) = r(3)
  r(3) += dr
 enddo
 mu_average_z_dens = 0.d0
 mu_average_z = 0.d0
 weight_average_z = 0.d0
 density_z = 0.d0
 do i_point = 1, n_points_final_grid
  r(:) = final_grid_points(:,i_point)
 !print*,'r(3)',r(3)
  if(r(3).lt.rmin)then
   i_tab = 1
  else if(r(3).gt.rmax)then
   i_tab = n_points_final_grid
  else
   integer :: i_tab
   i_tab = int((r(3) - rmin)/dr)+1
   i_tab = max(i_tab,1)
   i_tab = min(i_tab,n_points_final_grid)
  endif
   mu_average_z_dens(i_tab) += mu_of_r_vector(i_point) * final_weight_at_r_vector(i_point) * ( one_e_dm_alpha_at_r(i_point,1) + one_e_dm_beta_at_r(i_point,1) )
   mu_average_z(i_tab)      += mu_of_r_vector(i_point) * final_weight_at_r_vector(i_point) 
   weight_average_z(i_tab) += final_weight_at_r_vector(i_point) 
   density_z(1,i_tab) +=  one_e_dm_alpha_at_r(i_point,1)  * final_weight_at_r_vector(i_point)
   density_z(2,i_tab) +=  one_e_dm_beta_at_r(i_point,1)   * final_weight_at_r_vector(i_point)
  !print*,i_tab,z_tab(i_tab)
  !pause
 enddo
 
 r = 0.d0
 r(3) = rmin
 double precision :: accu
 accu = 0.d0
 do i_point = 1, n_points
! call f_HF_valence_ab(r,r,local_potential,two_bod)
  call f_HF_ab(r,r,local_potential,two_bod)
  call dm_dft_alpha_beta_at_r(r,dm_a,dm_b)
  if(two_bod.le.1.d-12.or.local_potential.le.0.d0.or.local_potential * two_bod.lt.0.d0)then
    local_potential = 1.d-10
  else
    local_potential = local_potential /  two_bod
  endif
  mu =  local_potential * dsqrt(dacos(-1.d0)) * 0.5d0

  call ESRC_MD_LDAERF (mu,dm_a,dm_b,.True.,ec)
  if(two_bod.ne.0.d0)then
   write(i_unit_output,'(100(F16.10,X))')r(3),mu,ec,local_potential,two_bod,dm_a+dm_b
  endif
  r(3) += dr
 !if(weight_average_z(i_point).gt.1.d-10)then
 ! dm_a = density_z(1,i_tab)/weight_average_z(i_point)
 ! dm_b = density_z(2,i_tab)/weight_average_z(i_point)
 ! mu = mu_average_z(i_point)/weight_average_z(i_point)
 ! if(.not. isnan(mu))then
 !  call ESRC_MD_LDAERF (mu,dm_a,dm_b,.True.,ec)
 ! endif
 !endif
 !if(weight_average_z(i_point).gt.1.d-10)then
 ! write(i_unit_output,'(100(F16.10,X))')z_tab(i_point),mu_average_z(i_point)/weight_average_z(i_point),ec*weight_average_z(i_point),mu_average_z_dens(i_point),density_z(1,i_point)
 !endif
 !accu += mu_average_z_dens(i_point) 
 enddo
!print*,'accu  = ',accu/dble(elec_alpha_num + elec_beta_num)
 print*,'mu av = ',mu_average
end

