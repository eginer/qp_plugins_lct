 BEGIN_PROVIDER [double precision, on_top_of_r_vector,(n_points_final_grid,N_states) ]
&BEGIN_PROVIDER [double precision, mu_of_r_cusp_condition_vector,(n_points_final_grid,N_states) ]
 implicit none
 integer :: i_point,istate
 double precision :: two_dm_in_r_selected_points,dpi,r(3),two_dm,two_dm_laplacian,total_dm
 double precision :: two_dm_HF,two_dm_laplacian_HF,total_dm_HF
 double precision :: corr_hole_2,alpha,mu0
 istate = 1
 double precision :: wall_0,wall_1
 double precision :: alpha_bis,delta,beta,mu
 print*,'providing the on_top_of_r_vector'
 i_point = 1
 dpi = 1.5d0 * dsqrt(dacos(-1.d0))
 call wall_time(wall_0)
 do i_point = 1, n_points_final_grid
  r(1) = final_grid_points(1,i_point)
  r(2) = final_grid_points(2,i_point)
  r(3) = final_grid_points(3,i_point)
  call spherical_averaged_two_dm_at_second_order(r,0.d0,istate,two_dm,two_dm_laplacian,total_dm)
  call spherical_averaged_two_dm_HF_at_second_order(r,0.d0,two_dm_HF,two_dm_laplacian_HF,total_dm_HF)
 !corr_hole_2 = (two_dm_laplacian - two_dm_laplacian_HF * two_dm/total_dm_HF)/two_dm_HF
  two_dm = max(two_dm,1.d-15)
  two_dm_HF = max(two_dm_HF,1.d-15)
  on_top_of_r_vector(i_point,istate) = two_dm
!! approximated polynom 
  mu0 =  dpi * (two_dm_laplacian / two_dm - two_dm_laplacian_HF / two_dm_HF)
  if(mu0.lt.0.d0)then
   print*,r
   print*,two_dm_laplacian / two_dm, two_dm_laplacian_HF / two_dm_HF, mu0
  endif
  mu0 = max(mu0,1.d-15)
  alpha = 1.d0 + 2.d0/(dsqrt(dacos(-1.d0)) * mu0)
  mu = mu0*alpha
!! new version with exact polynom 
!  alpha_bis = (two_dm_laplacian / two_dm - two_dm_laplacian_HF / two_dm_HF)
!  alpha_bis = max(alpha_bis,1.d-15)
!  beta = 2.d0/(3.d0*dacos(-1.d0))
!  delta = 2.d0/dsqrt(dacos(-1.d0))
!  mu = 1.d0/(2.d0*beta)*(alpha_bis + dsqrt(alpha_bis*alpha_bis + 4.d0 * alpha_bis * beta * delta))

  mu_of_r_cusp_condition_vector(i_point,istate) = mu
 enddo
 call wall_time(wall_1)
 print*,'provided the on_top_of_r_vector'
 print*,'Time to provide :',wall_1 - wall_0
 END_PROVIDER 
 

