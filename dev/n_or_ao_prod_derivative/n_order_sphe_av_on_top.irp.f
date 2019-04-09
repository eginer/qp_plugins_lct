
 BEGIN_PROVIDER[double precision, mos_nabla_n_in_r_array, (mo_num,mo_num,n_points_final_grid,0:order_derivative_ontop)]
&BEGIN_PROVIDER[double precision, mos_nabla_n_in_r_array_transp, (order_derivative_ontop,n_points_final_grid,mo_num,mo_num)]
 implicit none
 BEGIN_DOC
 ! blablablabla
 END_DOC
 integer :: i,j,k,m,n,n_exp
 double precision :: r(3)
 double precision :: delta_n_at_r(ao_num,ao_num)
 mos_nabla_n_in_r_array =0.d0
 mos_nabla_n_in_r_array_transp =0.d0
 do n_exp = 0,order_derivative_ontop 
  do i = 1, n_points_final_grid
   r(1) = final_grid_points(1,i)
   r(2) = final_grid_points(2,i)
   r(3) = final_grid_points(3,i) 
   call give_delta_n_at_r(r,n_exp,delta_n_at_r)
   do j = 1, mo_num
    do k = 1, mo_num
     do m = 1, ao_num
      do n = 1, ao_num
       mos_nabla_n_in_r_array(k,j,i,n_exp) += mo_coef(m,j)*mo_coef(n,k)*delta_n_at_r(n,m)
!!     mos_nabla_n_in_r_array_transp(i,k,j) += mo_coef(m,j)*mo_coef(n,k)*aos_nabla_4_in_r_array(n,m,i) 
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo
 END_PROVIDER



 BEGIN_PROVIDER[double precision, n2_sphe_av_exp,(n_points_final_grid,0:order_derivative_ontop,n_states)]
 implicit none
 BEGIN_DOC
 ! blablablabla
 END_DOC
 integer :: i,j,k,l,n,istate,r
 
 n2_sphe_av_exp = 0.d0 

 do istate = 1, n_states
  do n = 0, order_derivative_ontop      
   do r = 1, n_points_final_grid
    do l = 1, mo_num       
     do k = 1, mo_num        
      do j = 1, mo_num 
       do i = 1, mo_num 
        n2_sphe_av_exp(r,n,istate) += (2.d0/(2.d0*dble(n)+1.d0))*( two_bod_alpha_beta_mo_physicist(i,j,k,l,istate)*mos_in_r_array(i,r)* mos_in_r_array(k,r) * mos_nabla_n_in_r_array(j,l,r,n))
       enddo
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo  

 END_PROVIDER



 BEGIN_PROVIDER[double precision, coeff_Hx_exp,(0:order_derivative_ontop,n_states)]
 implicit none
 BEGIN_DOC
 ! blablablabla
 END_DOC
 include 'utils/constants.include.F'
 integer :: i,j,k,l,n,istate,r
 double precision :: fact
 
 coeff_Hx_exp = 0.d0 

 do istate = 1, n_states
  do n = 0, order_derivative_ontop      
   do r = 1, n_points_final_grid
    coeff_Hx_exp(n,istate) += n2_sphe_av_exp(r,n,istate) * final_weight_at_r_vector(r) 
   enddo
   coeff_Hx_exp(n,istate)= 2.d0*sqrt(pi) * gamma((2.d0*dble(n)+3.d0)/2.d0) / ( fact(2*n) * (2.d0*dble(n)+2)) * coeff_Hx_exp(n,istate) 
  enddo
 enddo  

 END_PROVIDER


 BEGIN_PROVIDER[double precision, e_hx_sr_by_order_exp,(n_states,0:order_derivative_ontop)]
 implicit none
 BEGIN_DOC
 ! blablablabla
 END_DOC
 integer :: n,istate
 
 e_hx_sr_by_order_exp= 0.d0 

 do istate = 1, n_states
  do n = 0, order_derivative_ontop      
   e_hx_sr_by_order_exp(istate,n) += coeff_Hx_exp(n,istate)/(mu_erf**(2.d0*dble(n)+2.d0))
  enddo
 enddo  

 END_PROVIDER

 BEGIN_PROVIDER[double precision, e_hx_sr_tot_exp,(n_states)]
 implicit none
 BEGIN_DOC
 ! blablablabla
 END_DOC
 integer :: n,istate
 
 e_hx_sr_tot_exp = 0.d0 

 do istate = 1, n_states
  do n = 0, order_derivative_ontop      
   e_hx_sr_tot_exp(istate) += coeff_Hx_exp(n,istate)/(mu_erf**(2.d0*dble(n)+2.d0))
  enddo
 enddo  

 END_PROVIDER



