
 BEGIN_PROVIDER[double precision, multipolar_integrals,(0:order_derivative_ontop,0:order_derivative_ontop,0:order_derivative_ontop,0:order_derivative_ontop,n_states)]
 implicit none
 BEGIN_DOC
 !Compute \int dr n(r)*r^{2a} * x^b * y^c * z^d
 END_DOC
 integer :: a,b,c,d,n,istate,r
 double precision :: dm,r2,r1(3)
 double precision :: tmp1,tmp2,tmp3
 double precision, allocatable :: x(:),y(:),z(:),r2n(:) 
 allocate(x(0:order_derivative_ontop),y(0:order_derivative_ontop),z(0:order_derivative_ontop),r2n(0:order_derivative_ontop))

 multipolar_integrals = 0.d0 

 do istate = 1, n_states
  do r = 1, n_points_final_grid
   r1(1) = final_grid_points(1,r)
   r1(2) = final_grid_points(2,r)
   r1(3) = final_grid_points(3,r)

   r2 = r1(1)*r1(1) + r1(2)*r1(2) + r1(3)*r1(3)

   dm = one_e_dm_alpha_at_r(r,istate)+one_e_dm_beta_at_r(r,istate)

   do n = 0, order_derivative_ontop
    x(n)  = r1(1)**dble(n) 
    y(n)  = r1(2)**dble(n) 
    z(n)  = r1(3)**dble(n) 
    r2n(n) = r2**dble(n) 
   enddo 

   do d = 0, order_derivative_ontop      
    tmp1 = z(d) * final_weight_at_r_vector(r) 
    do c = 0, order_derivative_ontop      
     tmp2 = tmp1 * y(c)
     do b = 0, order_derivative_ontop      
      tmp3 = tmp2 * x(b)
      do a = 0, order_derivative_ontop      
       multipolar_integrals(a,b,c,d,istate) += dm * r2n(a) * tmp3 
      enddo
     enddo
    enddo
   enddo

  enddo
 enddo

 END_PROVIDER



 BEGIN_PROVIDER[double precision, double_multipolar_integral,(0:order_derivative_ontop,n_states)]
 implicit none
 BEGIN_DOC
 !Compute \int \int dr1 dr2 n(r1)*n(r2) * r_{12}^{2n}
 END_DOC
 integer :: k,l,i,j,n,istate
 double precision :: tmp1,tmp2,tmp3
 double precision :: binom_func

 double_multipolar_integral = 0.d0 

 do istate = 1, n_states
  do n = 0, order_derivative_ontop
   do k = 0, n 
    tmp1 = binom_func(n,k) * (-2.d0)**(dble(n-k)) 
    do l = 0, k 
     tmp2 = tmp1 * binom_func(k,l)
     do i = 0, n-k 
      tmp3 = tmp2 * binom_func(n-k,i) 
      do j = 0, i 
       double_multipolar_integral(n,istate) +=  tmp3 * binom_func(i,j) * multipolar_integrals(l,j,i-j,n-k-i,istate)  * multipolar_integrals(k-l,j,i-j,n-k-i,istate) 
      enddo
     enddo
    enddo
   enddo

  enddo
 enddo

 END_PROVIDER



 BEGIN_PROVIDER[double precision, sr_hatree_energy_small_mu_coef,(0:order_derivative_ontop,n_states)]
&BEGIN_PROVIDER[double precision, sr_hatree_energy_small_mu_contrib,(0:order_derivative_ontop,n_states)]
 implicit none
 BEGIN_DOC
 !sr_hatree_energy_small_mu_coef(n) = -\frac{1}{\sqrt{pi}}*\frac{-1^n * a_n }{n!} * \int \int dr1 dr2 n(r1)*n(r2) * r_{12}^{2n}
 !sr_hatree_energy_small_mu_contrib(n) = -\frac{1}{\sqrt{pi}}*\frac{-1^n * a_n }{n!}* \mu^{2n+1}  * \int \int dr1 dr2 n(r1)*n(r2) * r_{12}^{2n}
 END_DOC
 include 'utils/constants.include.F'
 integer :: n,istate
 double precision :: a_n
 double precision :: fact 

 sr_hatree_energy_small_mu_coef = 0.d0 
 sr_hatree_energy_small_mu_contrib = 0.d0

 do istate = 1, n_states
  do n = 0, order_derivative_ontop
   a_n = 1.d0/(2.d0*dble(n)+1) 
   sr_hatree_energy_small_mu_coef(n,istate) = (((-1)**(dble(n)))*a_n/fact(n)) * double_multipolar_integral(n,istate)   
   sr_hatree_energy_small_mu_contrib(n,istate) = (mu_erf**(2.d0*dble(n)+1.d0)) * sr_hatree_energy_small_mu_coef(n,istate) 
  enddo
 enddo

 sr_hatree_energy_small_mu_coef = -1/sqrt(pi) * sr_hatree_energy_small_mu_coef
 sr_hatree_energy_small_mu_contrib = -1/sqrt(pi) * sr_hatree_energy_small_mu_contrib


 END_PROVIDER



 BEGIN_PROVIDER[double precision, sr_hatree_energy_small_mu_tot,(n_states)]
 implicit none
 BEGIN_DOC
 !Compute short range Hartree energy like eq B3 in of Toulouse et al PHYSICAL REVIEW A 70, 062505(2004)
 END_DOC
 integer :: n,istate

 sr_hatree_energy_small_mu_tot = 0.d0 

 do istate = 1, n_states
  do n = 0, order_derivative_ontop
   sr_hatree_energy_small_mu_tot(istate) +=  sr_hatree_energy_small_mu_contrib(n,istate) 
  enddo
  sr_hatree_energy_small_mu_tot(istate) = regular_range_Hartree(istate) + sr_hatree_energy_small_mu_tot(istate)
 enddo

 END_PROVIDER


 BEGIN_PROVIDER[double precision, e_h_pade_mu_0_diag_coef,(n_states)]
 implicit none
 BEGIN_DOC
 ! blablablabla
 END_DOC
 integer :: n,m,istate,INFO,i,j
 n=(2*order_derivative_ontop+1)/2  
 m=2*order_derivative_ontop+1-n
!n=(order_derivative_ontop)/2  
!m=order_derivative_ontop-n
 print*,"n =",n
 print*,"m =",m

 double precision :: denom,nume,x_loc
 double precision, allocatable :: matrix_alpha(:,:)
 double precision, allocatable :: vect_b_temp(:),vect_a(:),vect_b(:),c_pade(:)
 integer, allocatable :: IPIV(:)
 allocate(matrix_alpha(m,m),vect_b_temp(m),IPIV(m),vect_a(0:n),vect_b(0:m+1),c_pade(0:2*order_derivative_ontop+1))

 do istate = 1, n_states
 matrix_alpha= 0.d0
 vect_b_temp=0.d0

  c_pade = 0.d0
  do i = 0, order_derivative_ontop
   c_pade(2*i+1)=sr_hatree_energy_small_mu_coef(i,istate)
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
 
  !x_loc=1/(mu_erf**2.d0) 
 
  nume = 0.d0
  do i= 0,n
   !nume += vect_a(i)*x_loc**(dble(i))
   nume += vect_a(i)* mu_erf**(dble(i))
  enddo
  
  denom = 0.d0
  do i= 0,m
  !denom += vect_b(i)*x_loc**(dble(i))
   denom += vect_b(i)* mu_erf**(dble(i))
  enddo
 
 !e_hx_pade_diag_coef(istate) = 0.d0 
 e_h_pade_mu_0_diag_coef(istate) = nume/denom

 
 enddo  

 END_PROVIDER



