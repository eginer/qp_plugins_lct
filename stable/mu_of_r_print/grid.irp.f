 BEGIN_PROVIDER [integer, i_nucl_first]
&BEGIN_PROVIDER [integer, i_nucl_end]
 implicit none
 if(nucl_num.gt.1)then
  if(nucl_coord(1,3) .lt. nucl_coord(2,3))then
   i_nucl_first = 1
   i_nucl_end   = 2
  else
   i_nucl_first = 2
   i_nucl_end   = 1
  endif
 else 
   i_nucl_first = 1
   i_nucl_end   = 1
 endif

END_PROVIDER 


BEGIN_PROVIDER [ double precision, z_min_print  ]
  implicit none
  BEGIN_DOC
! beginning of the grid for mu(r) prints
  END_DOC

  z_min_print = nucl_coord(i_nucl_first,3) - spread_grid 

END_PROVIDER

BEGIN_PROVIDER [ double precision, z_max_print  ]
  implicit none
  BEGIN_DOC
! beginning of the grid for mu(r) prints
  END_DOC

  z_max_print = nucl_coord(i_nucl_end,3) + spread_grid 

END_PROVIDER

BEGIN_PROVIDER [double precision, dz_grid_mur]
 implicit none
 dz_grid_mur = (z_max_print - z_min_print)/dble(n_points_print_mur)

END_PROVIDER 

BEGIN_PROVIDER [double precision, grid_points_mur, (3, n_points_print_mur)]
 implicit none
 integer :: i
 double precision :: accu
 grid_points_mur = 0.d0
 accu = z_min_print
 do i = 1, n_points_print_mur
  accu += dz_grid_mur
  grid_points_mur(3,i)  = accu
 enddo
END_PROVIDER 

 BEGIN_PROVIDER[double precision, mos_on_grid_mur_array, (mo_num,n_points_print_mur)]
&BEGIN_PROVIDER[double precision, mos_on_grid_mur_array_transp,(n_points_print_mur,mo_num)]
 implicit none
 BEGIN_DOC
 ! mos_on_grid_mur_array(i,j)        = value of the ith mo on the jth grid point
 !
 ! mos_on_grid_mur_array_transp(i,j) = value of the jth mo on the ith grid point
 END_DOC
 integer :: i,j
 double precision :: mos_array(mo_num), r(3)
 do i = 1, n_points_print_mur
  r(1) = grid_points_mur(1,i)
  r(2) = grid_points_mur(2,i)
  r(3) = grid_points_mur(3,i)
  call give_all_mos_at_r(r,mos_array)
  do j = 1, mo_num
   mos_on_grid_mur_array(j,i) = mos_array(j)
   mos_on_grid_mur_array_transp(i,j) = mos_array(j)
  enddo
 enddo
 END_PROVIDER


 BEGIN_PROVIDER[double precision, core_inact_act_mos_in_r_array_grid_mur, (n_core_inact_act_orb,n_points_print_mur)]
&BEGIN_PROVIDER[double precision, core_inact_act_mos_in_r_array_grid_mur_transp,(n_points_print_mur,n_core_inact_act_orb)]
 implicit none
 integer :: i,j,k
 do i = 1, n_core_inact_act_orb
  j = list_core_inact_act(i)
  do k = 1, n_points_print_mur
   core_inact_act_mos_in_r_array_grid_mur_transp(k,i) = mos_on_grid_mur_array_transp(k,j)
  enddo
 enddo

 do k = 1, n_points_print_mur
  do i = 1, n_core_inact_act_orb
   core_inact_act_mos_in_r_array_grid_mur(i,k) = core_inact_act_mos_in_r_array_grid_mur_transp(k,i)
  enddo
 enddo

END_PROVIDER 

