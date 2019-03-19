 BEGIN_PROVIDER [ integer, n_z_points_grid]
&BEGIN_PROVIDER [ double precision, index_z_tab, (n_points_final_grid)]
&BEGIN_PROVIDER [ double precision, z_points_grid, (n_points_final_grid)]
 implicit none
 integer :: i,j,k,l
 double precision, allocatable :: z_points(:)
 integer, allocatable :: iorder(:)
 allocate(z_points(n_points_final_grid), iorder(n_points_final_grid))
 do i = 1, n_points_final_grid
  z_points(i) = final_grid_points(3,i)
  iorder(i) = i
 enddo
 call dsort(z_points,iorder, n_points_final_grid)
 n_z_points_grid = 1
 index_z_tab(iorder(1)) = n_z_points_grid
 z_points_grid(1) = z_points(1)
 do i = 2, n_points_final_grid
  if (z_points(i).ne.z_points(i-1))then
   n_z_points_grid += 1
   z_points_grid(n_z_points_grid) = z_points(i)
  endif
  index_z_tab(iorder(i)) = n_z_points_grid
 enddo

END_PROVIDER 

 BEGIN_PROVIDER [double precision, e_c_lda_val_hf_z_tab, (n_z_points_grid)]
&BEGIN_PROVIDER [double precision, e_c_lda_val_hf_z_sum, (n_z_points_grid)]
&BEGIN_PROVIDER [double precision, e_c_lda_val_hf_sum ]
 implicit none
 integer :: ipoint_z,i
 double precision :: mu,ec,dm_a,dm_b
 e_c_lda_val_hf_z_tab = 0.d0
 do i = 1, n_points_final_grid
  ipoint_z = index_z_tab(i) ! to which z point does it coincide
  mu = mu_of_r_hf_valencecoal_vector(i)
  dm_a = one_e_dm_no_core_and_grad_alpha_in_r(4,i,1)
  dm_b = one_e_dm_no_core_and_grad_beta_in_r(4,i,1)
  call ESRC_MD_LDAERF (mu,dm_a,dm_b,.True.,ec)
  e_c_lda_val_hf_z_tab(ipoint_z) += ec * final_weight_at_r_vector(i)
 enddo
 e_c_lda_val_hf_sum = 0.d0
 e_c_lda_val_hf_z_sum = 0.d0
 do i = 1, n_z_points_grid
  e_c_lda_val_hf_sum += e_c_lda_val_hf_z_tab(i) 
  e_c_lda_val_hf_z_sum(i) = e_c_lda_val_hf_sum
 enddo
 
END_PROVIDER

 BEGIN_PROVIDER [double precision, e_c_lda_ful_hf_z_tab, (n_z_points_grid)]
&BEGIN_PROVIDER [double precision, e_c_lda_ful_hf_z_sum, (n_z_points_grid)]
&BEGIN_PROVIDER [double precision, e_c_lda_ful_hf_sum ]
 implicit none
 integer :: ipoint_z,i
 double precision :: mu,ec,dm_a,dm_b
 e_c_lda_ful_hf_z_tab = 0.d0
 do i = 1, n_points_final_grid
  ipoint_z = index_z_tab(i) ! to which z point does it coincide
  mu = mu_of_r_hf_coal_vector(i)
  dm_a = one_e_dm_and_grad_alpha_in_r(4,i,1)
  dm_b = one_e_dm_and_grad_beta_in_r(4,i,1)
  call ESRC_MD_LDAERF (mu,dm_a,dm_b,.True.,ec)
  e_c_lda_ful_hf_z_tab(ipoint_z) += ec * final_weight_at_r_vector(i)
 enddo
 e_c_lda_ful_hf_sum = 0.d0
 do i = 1, n_z_points_grid
  e_c_lda_ful_hf_sum += e_c_lda_ful_hf_z_tab(i) 
  e_c_lda_ful_hf_z_tab(i) = e_c_lda_ful_hf_sum
 enddo
 
END_PROVIDER





 BEGIN_PROVIDER [double precision, z_max_grid_cyl]
&BEGIN_PROVIDER [double precision, z_min_grid_cyl]
&BEGIN_PROVIDER [integer, n_z_grid_cyl]
&BEGIN_PROVIDER [double precision, dz_grid_cyl]
 implicit none
 z_min_grid_cyl = nucl_coord(1,3) - 4.5d0
 z_max_grid_cyl = nucl_coord(1,3) + 7.5d0
 n_z_grid_cyl = 1000
 dz_grid_cyl = ( z_max_grid_cyl - z_min_grid_cyl ) /dble(n_z_grid_cyl)
END_PROVIDER 

 BEGIN_PROVIDER [double precision, z_points_grid_cylindr, (n_z_grid_cyl)]
&BEGIN_PROVIDER [integer, index_z0 ]
 implicit none
 integer :: i,j 
 logical :: test
 test = .True.
 z_points_grid_cylindr(1) = z_min_grid_cyl
 do i = 2, n_z_grid_cyl
  z_points_grid_cylindr(i) = z_points_grid_cylindr(i-1) + dz_grid_cyl
  if(z_points_grid_cylindr(i).gt.0.d0.and.test)then
   index_z0 = i
   test = .False.
  endif
 enddo
END_PROVIDER 

 BEGIN_PROVIDER [integer, n_radial_points_grid_cyl]
&BEGIN_PROVIDER [double precision, rmax_grid_cyl ]
&BEGIN_PROVIDER [double precision, dr_grid_cyl]
 implicit none
 n_radial_points_grid_cyl = 900
 rmax_grid_cyl = 9.D0
 dr_grid_cyl = rmax_grid_cyl / dble(n_radial_points_grid_cyl)
END_PROVIDER 

 BEGIN_PROVIDER [integer, n_angular_points_grid_cyl]
&BEGIN_PROVIDER [double precision, d_theta_grid_cyl]
 implicit none
 double precision :: pi
 pi = dacos(-1.D0)
 n_angular_points_grid_cyl = 2 
 d_theta_grid_cyl = 2.d0 * pi / dble(n_angular_points_grid_cyl) 
END_PROVIDER 

 BEGIN_PROVIDER [double precision, grid_points_grid_cyl, (3,n_angular_points_grid_cyl, n_radial_points_grid_cyl , n_z_grid_cyl )]
 implicit none
 integer :: i,j,k
 do i = 1, n_z_grid_cyl
  do j = 1, n_radial_points_grid_cyl
   do k = 1, n_angular_points_grid_cyl
    grid_points_grid_cyl(1,k,j,i) = dble(j-1) * dr_grid_cyl * dcos((k-1) * d_theta_grid_cyl)
    grid_points_grid_cyl(2,k,j,i) = dble(j-1) * dr_grid_cyl * dsin((k-1) * d_theta_grid_cyl)
    grid_points_grid_cyl(3,k,j,i) = z_points_grid_cylindr(i)
   enddo
  enddo
 enddo

END_PROVIDER 

BEGIN_PROVIDER [double precision, mu_ful_hf_grid_cyl, ( n_angular_points_grid_cyl, n_radial_points_grid_cyl , n_z_grid_cyl )]
 implicit none 
 integer :: i,j,k
 double precision :: local_potential,two_bod,r(3)
 mu_ful_hf_grid_cyl = 0.d0
 do i = 1, n_z_grid_cyl
  do j = 1, n_radial_points_grid_cyl
   do k = 1, n_angular_points_grid_cyl
    r(:) = grid_points_grid_cyl(:,k,j,i)
    call f_HF_ab(r,r,local_potential,two_bod)
    if(two_bod.le.1.d-12.or.local_potential.le.0.d0.or.local_potential * two_bod.lt.0.d0)then
      local_potential = 0.d+00
    else 
      local_potential = local_potential /  two_bod
    endif
    mu_ful_hf_grid_cyl(k,j,i) =  local_potential * dsqrt(dacos(-1.d0)) * 0.5d0
   enddo
  enddo
 enddo

END_PROVIDER


BEGIN_PROVIDER [double precision, mu_val_hf_grid_cyl, ( n_angular_points_grid_cyl, n_radial_points_grid_cyl , n_z_grid_cyl )]
 implicit none 
 integer :: i,j,k
 double precision :: local_potential,two_bod,r(3)
 mu_val_hf_grid_cyl = 0.d0
 do i = 1, n_z_grid_cyl
  do j = 1, n_radial_points_grid_cyl
   do k = 1, n_angular_points_grid_cyl
    r(:) = grid_points_grid_cyl(:,k,j,i)
    call f_HF_valence_ab(r,r,local_potential,two_bod)
    if(two_bod.le.1.d-12.or.local_potential.le.0.d0.or.local_potential * two_bod.lt.0.d0)then
      local_potential = 0.d+00
    else 
      local_potential = local_potential /  two_bod
    endif
    mu_val_hf_grid_cyl(k,j,i) =  local_potential * dsqrt(dacos(-1.d0)) * 0.5d0
   enddo
  enddo
 enddo

END_PROVIDER


 BEGIN_PROVIDER [double precision, one_e_dm_alpha_grid_cyl , ( n_angular_points_grid_cyl, n_radial_points_grid_cyl , n_z_grid_cyl )]
&BEGIN_PROVIDER [double precision, one_e_dm_beta_grid_cyl , ( n_angular_points_grid_cyl, n_radial_points_grid_cyl , n_z_grid_cyl )]
&BEGIN_PROVIDER [double precision, one_e_dm_alpha_no_core_grid_cyl , ( n_angular_points_grid_cyl, n_radial_points_grid_cyl , n_z_grid_cyl )]
&BEGIN_PROVIDER [double precision, one_e_dm_beta_no_core_grid_cyl , ( n_angular_points_grid_cyl, n_radial_points_grid_cyl , n_z_grid_cyl )]
 implicit none 
 integer :: i,j,k
 double precision :: dm_a, dm_b, dm_a_no_core, dm_b_no_core,r(3)
 one_e_dm_alpha_grid_cyl = 0.d0
 one_e_dm_beta_grid_cyl  = 0.d0
 one_e_dm_alpha_no_core_grid_cyl = 0.d0
 one_e_dm_beta_no_core_grid_cyl  = 0.d0
 do i = 1, n_z_grid_cyl
  do j = 1, n_radial_points_grid_cyl
   do k = 1, n_angular_points_grid_cyl
    r(:) = grid_points_grid_cyl(:,k,j,i)
    call dm_dft_alpha_beta_at_r(r,dm_a,dm_b)
    one_e_dm_alpha_grid_cyl(k,j,i) = dm_a
    one_e_dm_beta_grid_cyl(k,j,i)  = dm_b
    call dm_dft_alpha_beta_no_core_at_r(r,dm_a_no_core,dm_b_no_core)
    one_e_dm_alpha_no_core_grid_cyl(k,j,i) = dm_a_no_core
    one_e_dm_beta_no_core_grid_cyl(k,j,i)  = dm_b_no_core
   enddo
  enddo
 enddo

END_PROVIDER

 BEGIN_PROVIDER [double precision, e_c_lda_hf_ful_grid_cyl, ( n_angular_points_grid_cyl, n_radial_points_grid_cyl , n_z_grid_cyl )]
&BEGIN_PROVIDER [double precision, e_c_lda_hf_val_grid_cyl, ( n_angular_points_grid_cyl, n_radial_points_grid_cyl , n_z_grid_cyl )]
&BEGIN_PROVIDER [double precision, e_c_lda_hf_ful_grid_cyl_z_tab, ( n_z_grid_cyl )]
&BEGIN_PROVIDER [double precision, e_c_lda_hf_val_grid_cyl_z_tab, ( n_z_grid_cyl )]
&BEGIN_PROVIDER [double precision, dm_grid_cyl_z_tab, ( n_z_grid_cyl )]
&BEGIN_PROVIDER [double precision, dm_no_core_grid_cyl_z_tab, ( n_z_grid_cyl )]
&BEGIN_PROVIDER [double precision, mu_hf_ful_grid_cyl_z_tab, ( n_z_grid_cyl )]
&BEGIN_PROVIDER [double precision, mu_hf_val_grid_cyl_z_tab, ( n_z_grid_cyl )]
&BEGIN_PROVIDER [double precision, e_c_lda_hf_val_grid_cyl_z_sum ]
&BEGIN_PROVIDER [double precision, e_c_lda_hf_ful_grid_cyl_z_sum ]
&BEGIN_PROVIDER [double precision, elec_alpha_num_grid_cyl]
&BEGIN_PROVIDER [double precision, elec_beta_num_grid_cyl]
&BEGIN_PROVIDER [double precision, elec_alpha_num_no_core_grid_cyl]
&BEGIN_PROVIDER [double precision, elec_beta_num_no_core_grid_cyl]
 implicit none
 integer :: i,j,k
 double precision :: weight_tmp
 double precision :: mu,ec,dm_a,dm_b
 double precision :: mu_val,ec_no_core,dm_a_no_core,dm_b_no_core
 e_c_lda_hf_ful_grid_cyl_z_tab = 0.d0
 e_c_lda_hf_ful_grid_cyl_z_sum = 0.d0
 e_c_lda_hf_val_grid_cyl_z_sum = 0.d0
 e_c_lda_hf_val_grid_cyl_z_tab = 0.d0
 elec_alpha_num_grid_cyl = 0.d0
 elec_beta_num_grid_cyl = 0.d0
 elec_alpha_num_no_core_grid_cyl = 0.d0
 elec_beta_num_no_core_grid_cyl = 0.d0
 mu_hf_ful_grid_cyl_z_tab = 0.d0
 mu_hf_val_grid_cyl_z_tab = 0.d0
 dm_grid_cyl_z_tab = 0.d0
 dm_no_core_grid_cyl_z_tab = 0.d0
 do i = 1, n_z_grid_cyl
  do j = 1, n_radial_points_grid_cyl
   do k = 1, n_angular_points_grid_cyl
    mu = mu_ful_hf_grid_cyl(k,j,i)
    dm_a = one_e_dm_alpha_grid_cyl(k,j,i)
    dm_b = one_e_dm_beta_grid_cyl(k,j,i)
    call ESRC_MD_LDAERF (mu,dm_a,dm_b,.True.,ec)
    e_c_lda_hf_ful_grid_cyl(k,j,i) = ec
    weight_tmp = dble(j-1) * dr_grid_cyl * dr_grid_cyl * d_theta_grid_cyl

    e_c_lda_hf_ful_grid_cyl_z_tab(i) += ec * dble(j-1) * dr_grid_cyl * dr_grid_cyl * d_theta_grid_cyl
    mu_hf_ful_grid_cyl_z_tab(i) += mu * (dm_a + dm_b) * dble(j-1) * dr_grid_cyl * dr_grid_cyl * d_theta_grid_cyl
    dm_grid_cyl_z_tab(i) += (dm_a + dm_b) * dble(j-1) * dr_grid_cyl * dr_grid_cyl * d_theta_grid_cyl 
    e_c_lda_hf_ful_grid_cyl_z_sum += ec * weight_tmp * dz_grid_cyl
    elec_beta_num_grid_cyl += dm_b * weight_tmp * dz_grid_cyl
    elec_alpha_num_grid_cyl += dm_a * weight_tmp * dz_grid_cyl

    mu_val = mu_val_hf_grid_cyl(k,j,i)
    dm_a_no_core = one_e_dm_alpha_no_core_grid_cyl(k,j,i)
    dm_b_no_core = one_e_dm_beta_no_core_grid_cyl(k,j,i)
    call ESRC_MD_LDAERF (mu_val,dm_a_no_core,dm_b_no_core,.True.,ec_no_core)
    e_c_lda_hf_val_grid_cyl(k,j,i) = ec_no_core
    weight_tmp = dble(j-1) * dr_grid_cyl * dr_grid_cyl * d_theta_grid_cyl

    dm_no_core_grid_cyl_z_tab(i) += (dm_a_no_core + dm_b_no_core) * dble(j-1) * dr_grid_cyl * dr_grid_cyl * d_theta_grid_cyl 
    mu_hf_val_grid_cyl_z_tab(i) += mu_val * (dm_a_no_core + dm_b_no_core) * dble(j-1) * dr_grid_cyl * dr_grid_cyl * d_theta_grid_cyl
    e_c_lda_hf_val_grid_cyl_z_tab(i) += ec_no_core * dble(j-1) * dr_grid_cyl * dr_grid_cyl * d_theta_grid_cyl
    e_c_lda_hf_val_grid_cyl_z_sum += ec_no_core * weight_tmp * dz_grid_cyl
    elec_beta_num_no_core_grid_cyl += dm_b_no_core * weight_tmp * dz_grid_cyl
    elec_alpha_num_no_core_grid_cyl += dm_a_no_core * weight_tmp * dz_grid_cyl

   enddo
  enddo
 enddo
 mu_hf_ful_grid_cyl_z_tab = mu_hf_ful_grid_cyl_z_tab / (elec_beta_num_grid_cyl + elec_alpha_num_grid_cyl)
 mu_hf_val_grid_cyl_z_tab = mu_hf_val_grid_cyl_z_tab / (elec_beta_num_no_core_grid_cyl + elec_alpha_num_no_core_grid_cyl)

 END_PROVIDER 

 BEGIN_PROVIDER [integer, imax_ec]
 implicit none 
 imax_ec = minloc(e_c_lda_hf_val_grid_cyl_z_tab, DIM=1)
 print*,'imax_ec = ',imax_ec

 END_PROVIDER 

 BEGIN_PROVIDER [double precision, e_c_lda_ful_sym, ( n_z_grid_cyl )  ]
 implicit none
 integer :: i,j
 integer :: iref
 iref = 549
 e_c_lda_ful_sym = 0.d0
 do i = 1, n_z_grid_cyl
  j = i + iref - imax_ec
  if(j.le.n_z_grid_cyl.and.j.gt.0)then
   e_c_lda_ful_sym(j) = e_c_lda_hf_ful_grid_cyl_z_tab(i)
  endif
 enddo

 END_PROVIDER 

 BEGIN_PROVIDER [double precision, e_c_lda_val_sym, ( n_z_grid_cyl )  ]
 implicit none
 integer :: i,j
 integer :: iref
 iref = 549
 e_c_lda_val_sym = 0.d0
 do i = 1, n_z_grid_cyl
  j = i + iref - imax_ec
  if(j.le.n_z_grid_cyl.and.j.gt.0)then
   e_c_lda_val_sym(j) = e_c_lda_hf_val_grid_cyl_z_tab(i)
  endif
 enddo

 END_PROVIDER 

 BEGIN_PROVIDER [double precision, mu_hf_ful_sym, ( n_z_grid_cyl )  ]
 implicit none
 integer :: i,j
 integer :: iref
 iref = 549
 mu_hf_ful_sym = 0.d0
 do i = 1, n_z_grid_cyl
  j = i + iref - imax_ec
  if(j.le.n_z_grid_cyl.and.j.gt.0)then
   mu_hf_ful_sym(j) = mu_hf_ful_grid_cyl_z_tab(i)
  endif
 enddo

 END_PROVIDER 

 BEGIN_PROVIDER [double precision, mu_hf_val_sym, ( n_z_grid_cyl )  ]
 implicit none
 integer :: i,j
 integer :: iref
 iref = 549
 mu_hf_val_sym = 0.d0
 do i = 1, n_z_grid_cyl
  j = i + iref - imax_ec
  if(j.le.n_z_grid_cyl.and.j.gt.0)then
   mu_hf_val_sym(j) = mu_hf_val_grid_cyl_z_tab(i)
  endif
 enddo

 END_PROVIDER 

 BEGIN_PROVIDER [double precision, dm_ful_sym, ( n_z_grid_cyl )  ]
 implicit none
 integer :: i,j
 integer :: iref
 iref = 549
 dm_ful_sym = 0.d0
 do i = 1, n_z_grid_cyl
  j = i + iref - imax_ec
  if(j.le.n_z_grid_cyl.and.j.gt.0)then
   dm_ful_sym(j) = dm_grid_cyl_z_tab(i)
  endif
 enddo

 END_PROVIDER 

 BEGIN_PROVIDER [double precision, dm_val_sym, ( n_z_grid_cyl )  ]
 implicit none
 integer :: i,j
 integer :: iref
 iref = 549
 dm_val_sym = 0.d0
 do i = 1, n_z_grid_cyl
  j = i + iref - imax_ec
  if(j.le.n_z_grid_cyl.and.j.gt.0)then
   dm_val_sym(j) = dm_no_core_grid_cyl_z_tab(i)
  endif
 enddo

 END_PROVIDER 
