

 BEGIN_PROVIDER[double precision, coeff_H_exp,(0:order_derivative_ontop,n_states)]
 implicit none
 BEGIN_DOC
 ! blablablabla
 END_DOC
 include 'utils/constants.include.F'
 integer :: i,j,k,l,n,istate,r
 double precision :: fact
 double precision :: dm_a,dm_b,r1(3)

 coeff_H_exp = 0.d0 

 do istate = 1, n_states
  do n = 0, order_derivative_ontop      
   do r = 1, n_points_final_grid
    r1(1) = final_grid_points(1,r)
    r1(2) = final_grid_points(2,r)
    r1(3) = final_grid_points(3,r)
    call dm_dft_alpha_beta_at_r(r1,dm_a,dm_b)

   !coeff_H_exp(n,istate) += (one_e_dm_alpha_at_r(r,istate)+one_e_dm_beta_at_r(r,istate))* n_sphe_av_exp(r,n,istate) * final_weight_at_r_vector(r) 
    coeff_H_exp(n,istate) += (dm_a+dm_b)* n_sphe_av_exp(r,n,istate) * final_weight_at_r_vector(r) 
   enddo
   coeff_H_exp(n,istate)= 2.d0*sqrt(pi) * gamma((2.d0*dble(n)+3.d0)/2.d0) / ( fact(2*n) * (2.d0*dble(n)+2)) * coeff_H_exp(n,istate) 
  enddo
 enddo  

 END_PROVIDER


 BEGIN_PROVIDER[double precision, e_h_sr_by_order_exp,(n_states,0:order_derivative_ontop)]
 implicit none
 BEGIN_DOC
 ! blablablabla
 END_DOC
 integer :: n,istate
 
 e_h_sr_by_order_exp= 0.d0 

 do istate = 1, n_states
  do n = 0, order_derivative_ontop      
   e_h_sr_by_order_exp(istate,n) += coeff_H_exp(n,istate)/(mu_erf**(2.d0*dble(n)+2.d0))
  enddo
 enddo  

 END_PROVIDER

 BEGIN_PROVIDER[double precision, e_h_sr_tot_exp,(n_states)]
 implicit none
 BEGIN_DOC
 ! blablablabla
 END_DOC
 integer :: n,istate
 
 e_h_sr_tot_exp = 0.d0 

 do istate = 1, n_states
  do n = 0, order_derivative_ontop      
   e_h_sr_tot_exp(istate) += coeff_H_exp(n,istate)/(mu_erf**(2.d0*dble(n)+2.d0))
  enddo
 enddo  

 END_PROVIDER



 BEGIN_PROVIDER[double precision, e_h_pade_diag_coef,(n_states)]
 implicit none
 BEGIN_DOC
 ! blablablabla
 END_DOC
 integer :: n,m,istate,INFO,i,j
 n=(order_derivative_ontop+1)/2  
 m=order_derivative_ontop+1-n
 print*,"n =",n
 print*,"m =",m
!n=2
!m=1
!print*,"n =",n
!print*,"m =",m
 double precision :: denom,nume,x_loc
 double precision, allocatable :: matrix_alpha(:,:)
 double precision, allocatable :: vect_b_temp(:),vect_a(:),vect_b(:),c_pade(:)
 integer, allocatable :: IPIV(:)
 allocate(matrix_alpha(m,m),vect_b_temp(m),IPIV(m),vect_a(0:n),vect_b(0:m+1),c_pade(0:n+m+1))

 do istate = 1, n_states
 matrix_alpha= 0.d0
 vect_b_temp=0.d0

  c_pade = 0.d0
  do i = 1,order_derivative_ontop+1
   c_pade(i)=coeff_H_exp(i-1,istate)
  enddo 


  do j =1,m
   vect_b_temp(j) = - c_pade(n+j)
   do i = 1,m
   !matrix_alpha(i,j)= 0.d0
    matrix_alpha(i,j)=c_pade(n+i-j)
   !print*,i,j,matrix_alpha(i,j)
   !print*,j,vect_b_temp(j)
   enddo
  enddo
 
  call dgesv(m,1,matrix_alpha,m,IPIV,vect_b_temp,m,INFO) 
 
  print*,'***********************'  
  print*,'INFO=',INFO  
  print*,'***********************'  

  vect_b = 0.d0
  vect_b(0) = 1.d0
  do i=1,m
   vect_b(i)=vect_b_temp(i)
  enddo

  vect_a(0) = 0.d0
  do i = 1,n
   do j = 0,i
    vect_a(i) += c_pade(i-j)*vect_b(j)
   enddo
  enddo
 
  x_loc=1/(mu_erf**2.d0) 
 
  nume = 0.d0
  do i= 0,n
   nume += vect_a(i)*x_loc**(dble(i))
  enddo
  
  denom = 0.d0
  do i= 0,m
   denom += vect_b(i)*x_loc**(dble(i))
  enddo
 
 !e_hx_pade_diag_coef(istate) = 0.d0 
  e_h_pade_diag_coef(istate) = nume/denom

 
 enddo  

 END_PROVIDER
