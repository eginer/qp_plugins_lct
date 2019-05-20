program projected_operators
  implicit none
  BEGIN_DOC
! TODO
  END_DOC
  read_wf = .True.
  touch read_wf
 no_core_density = "no_core_dm"
 touch no_core_density

! provide mo_two_e_integrals_jj
  provide mo_class
! call routine_v
! call routine_rho 
! call routine_final
! call routine_valence
! call routine_core
! call routine_core_valence
! call test_f_hf_ao
! call test_f_hf_ao_per_atom
! call test_ovlp
  provide f_hf_ab_ao_per_atom

end

subroutine routine_core_valence
 implicit none
 integer :: ipoint,k,l,i,j
 double precision :: accu_core_val,accu_ful, weight,r(3),integral_psi_core_val,integral_psi,r2(3),two_bod
 accu_core_val = 0.d0
 accu_ful = 0.d0
 do ipoint  = 1, n_points_final_grid
  weight=final_weight_at_r_vector(ipoint)
  r(1) = final_grid_points(1,ipoint)
  r(2) = final_grid_points(2,ipoint)
  r(3) = final_grid_points(3,ipoint)
  call integral_f_HF_core_valence_ab(r,integral_psi_core_val)
  accu_core_val += integral_psi_core_val * weight
 enddo

 print*,'accu_core_val        = ',accu_core_val

 double precision :: accu_2
 integer :: i_i,j_j
 accu_2 = 0.d0
 !'alpha core + beta val'
 do i = 1, n_core_orb
  i_i = list_core(i)
  do j = 1, n_valence_orb_for_hf(2)
   j_j = list_valence_orb_for_hf(j,2)
   accu_2 += mo_two_e_integrals_jj(j_j,i_i)
  enddo
 enddo

 !'beta  core + alpha val'
 do i = 1, n_core_orb
  i_i = list_core(i)
  do j = 1, n_valence_orb_for_hf(1)
   j_j = list_valence_orb_for_hf(j,1)
   accu_2 += mo_two_e_integrals_jj(j_j,i_i)
  enddo
 enddo

 print*,'accu_2               = ',accu_2

 print*,'***************'
 accu_core_val = 0.d0
 r = 0.d0
 r(3) = 0.137d0
 call integral_f_HF_core_valence_ab(r,integral_psi)
 do ipoint  = 1, n_points_final_grid
  weight=final_weight_at_r_vector(ipoint)
  r2(1) = final_grid_points(1,ipoint)
  r2(2) = final_grid_points(2,ipoint)
  r2(3) = final_grid_points(3,ipoint)
  call f_HF_core_valence_ab(r,r2,integral_psi_core_val,two_bod)
  accu_core_val += integral_psi_core_val * weight
 enddo
 !print*,'integral_f_hf        = ',integral_f_hf
 print*,'accu_core_val        = ',accu_core_val
 print*,'integral_psi         = ',integral_psi
end


subroutine routine_core
 implicit none
 integer :: ipoint,k,l,i,j
 double precision :: accu_core,accu_ful, weight,r(3),integral_psi_core,integral_psi,r2(3)
 double precision :: accu_2,two_bod
 accu_core = 0.d0
 do ipoint  = 1, n_points_final_grid
  weight=final_weight_at_r_vector(ipoint)
  r(1) = final_grid_points(1,ipoint)
  r(2) = final_grid_points(2,ipoint)
  r(3) = final_grid_points(3,ipoint)
  call integral_f_HF_core_ab(r,integral_psi_core)
  accu_core += integral_psi_core * weight
 enddo
 print*,'accu_core            = ',accu_core
 accu_2 = 0.d0
 ! alpha beta 
 do i = 1, n_core_orb
  do j = 1, n_core_orb
   accu_2 += mo_two_e_integrals_jj(j,i)
  enddo
 enddo
 print*,'accu_2               = ',accu_2

 print*,'***************'
 r = 0.d0
 call integral_f_HF_core_ab(r,integral_psi)
 accu_core = 0.d0
 do ipoint  = 1, n_points_final_grid
  weight=final_weight_at_r_vector(ipoint)
  r2(1) = final_grid_points(1,ipoint)
  r2(2) = final_grid_points(2,ipoint)
  r2(3) = final_grid_points(3,ipoint)
  call f_HF_core_ab(r,r2,integral_psi_core,two_bod)
  accu_core += integral_psi_core * weight
 enddo
 print*,'accu_core            = ',accu_core
 print*,'integral_psi         = ',integral_psi
 print*,'***************'
end


subroutine routine_valence
 implicit none
 integer :: ipoint,k,l,i,j,i_i,j_j
 double precision :: accu_val,accu_ful, weight,r(3),integral_psi_val,integral_psi,r2(3),two_bod
 accu_val = 0.d0
 do ipoint  = 1, n_points_final_grid
  weight=final_weight_at_r_vector(ipoint)
  r(1) = final_grid_points(1,ipoint)
  r(2) = final_grid_points(2,ipoint)
  r(3) = final_grid_points(3,ipoint)
  call integral_f_HF_valence_ab(r,integral_psi_val)
  accu_val += integral_psi_val * weight
 enddo
 print*,'**************************'
 print*,'accu_val             = ',accu_val
 double precision :: accu_2
 accu_2 = 0.d0
 do i = 1, n_valence_orb_for_hf(1)
  i_i = list_valence_orb_for_hf(i,1)
  do j = 1, n_valence_orb_for_hf(2)
   j_j = list_valence_orb_for_hf(j,2)
   accu_2 += mo_two_e_integrals_jj(j_j,i_i)
  enddo
 enddo
 print*,'accu_2               = ',accu_2
 print*,'**************************'
 accu_val = 0.d0
 r = 0.d0
 call integral_f_HF_valence_ab(r,integral_psi)
 do ipoint  = 1, n_points_final_grid
  weight=final_weight_at_r_vector(ipoint)
  r2(1) = final_grid_points(1,ipoint)
  r2(2) = final_grid_points(2,ipoint)
  r2(3) = final_grid_points(3,ipoint)
  call f_HF_valence_ab(r,r2,integral_psi_val,two_bod)
  accu_val += integral_psi_val * weight
 enddo
 print*,'integral_psi         = ',integral_psi
 print*,'accu_val             = ',accu_val
 print*,'**************************'
end

subroutine routine_rho
 implicit none
 print*,'integral_r1r2_f_HF_aa = ',integral_r1r2_f_HF_aa
 print*,'psi_energy_two_e      = ',psi_energy_two_e

end
subroutine routine_v
 implicit none
 integer :: ipoint,k,l
 double precision :: accu, weight,r(3),integral_of_f_PSI_ab_over_2,f_HF_aa_integrated
 accu = 0.d0
 do ipoint  = 1, n_points_final_grid
  weight=final_weight_at_r_vector(ipoint)
  r(1) = final_grid_points(1,ipoint)
  r(2) = final_grid_points(2,ipoint)
  r(3) = final_grid_points(3,ipoint)
  accu += f_HF_aa_integrated(r) * weight
! accu += integral_of_f_PSI_ab_over_2(r) * weight
 enddo
 print*,'accu                   = ',accu
 print*,'psi_energy_two_e      = ',psi_energy_two_e
!print*,'integral_f_hf          = ',integral_f_hf
 accu = 0.d0
 integer :: i,j
 do i = 1, mo_num
  do j = 1, mo_num
   accu +=  two_bod_alpha_beta_mo(j,j,i,i,1)
  enddo
 enddo
 print*,'accu = ',accu
end
 

subroutine routine_final
 implicit none
 integer :: i,j,k,nx
 double precision :: dx,xmax,r12
 double precision :: r1(3),r2(3),numerator,denominator,f_HF_aa,HF_two_body_dm_aa,potential
 xmax = 2.d0
 nx = 10000
 dx = xmax/dble(nx)
 r1 = 0.d0
 r2 = r1
 
 double precision :: laplacian_f,value_f_HF_aa
 double precision :: laplacian_hf,HF_two_bod
 call f_HF_aa_spherical_averaged(r1,0.d0,laplacian_f,value_f_HF_aa)
 call HF_two_body_dm_aa_spherical_laplacian(r1,0.d0,HF_two_bod,laplacian_HF)
 double precision :: mu
 print*,'laplacian_f = ',laplacian_f
 print*,'laplacian_HF= ',laplacian_HF
 mu = laplacian_f/laplacian_HF
 print*,'mu = ',mu
 do i = 1, nx
  r2(1) += dx
  r2(2) += dx
  r12 = dsqrt((r1(1) - r2(1))**2 + (r1(2) - r2(2))**2  + (r1(3) - r2(3))**2 )
  numerator = f_HF_aa(r1,r2) 
  denominator = HF_two_body_dm_aa(r1,r2)
  potential = numerator/denominator
  write(33,'(100(F16.10,X))')r12,potential,erf(mu * r12)/r12,1.d0/r12,numerator,denominator,numerator/denominator
 enddo

end

subroutine test_f
 implicit none
 integer :: ipoint
 double precision :: accuex,accuc,accu,weight
 accuex = 0.d0
 accuc  = 0.d0
 accu   = 0.d0
 do ipoint  = 1, n_points_final_grid
  weight=final_weight_at_r_vector(ipoint)
  accuc += core_inact_act_f_psi_ab(ipoint) * weight
  accuex += f_psi_ab(ipoint) * weight
  accu += dabs(accuc - accuex)
 enddo
 print*,'accuex = ',accuex
 print*,'accuc  = ',accuc 
 print*,'accu   = ',accu  
end


subroutine test_f_hf_ao
 implicit none
 integer :: ipoint,i,inucl
 double precision :: accuex,accuerror,accuao,weight, r(3)
 double precision :: integral_psi,two_bod,accu_beta
 accuex      = 0.d0
 accuao      = 0.d0
 accuerror   = 0.d0
 accu_beta   = 0.d0
 do ipoint  = 1, n_points_final_grid
  weight=final_weight_at_r_vector(ipoint)
  r(:) = final_grid_points(:,ipoint) 
  call f_HF_valence_ab(r,r,integral_psi,two_bod)
  accuex += integral_psi * weight
  accuao += f_hf_ab_ao(ipoint) * weight 
  accuerror += dabs(integral_psi - f_hf_ab_ao(ipoint)) * weight
 enddo
 print*,'accuex      = ',accuex
 print*,'accuao      = ',accuao
 print*,'accuerror   = ',accuerror
end


subroutine test_f_hf_ao_per_atom
 implicit none
 integer :: ipoint,i,inucl
 double precision :: accuex,accuerror,accuao,weight, r(3)
 double precision :: integral_psi,two_bod,accu_beta
 accuex      = 0.d0
 accuao      = 0.d0
 accuerror   = 0.d0
 accu_beta   = 0.d0
 do inucl = 1, nucl_num
  do ipoint  = 1, n_pts_per_atom(inucl)
   weight = final_weight_at_r_vector_per_atom(ipoint,inucl)
   r(:) = final_grid_points_per_atom(:,ipoint,inucl) 
   call f_HF_valence_ab(r,r,integral_psi,two_bod)
   accuex += integral_psi * weight
   accuao += f_hf_ab_ao_per_atom(ipoint,inucl) * weight 
   accuerror += dabs(integral_psi - f_hf_ab_ao_per_atom(ipoint,inucl)) * weight
  enddo
 enddo
 print*,'accuex      = ',accuex
 print*,'accuao      = ',accuao
 print*,'accuerror   = ',accuerror
end


subroutine test_ovlp
 implicit none
 provide n_good_pairs_density_per_atom



end
