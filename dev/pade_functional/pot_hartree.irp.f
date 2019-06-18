program pade_functional
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  read_wf = .true.
  touch read_wf
  io_mo_one_e_integrals = "None"
  touch io_mo_one_e_integrals
  io_mo_two_e_integrals = "None"
  touch io_mo_two_e_integrals
  io_ao_two_e_integrals = "None"
  touch io_ao_two_e_integrals
 
  io_mo_integrals_e_n = "None"
  touch io_mo_integrals_e_n
  io_mo_integrals_kinetic = "None"
  touch io_mo_integrals_kinetic 
  io_ao_integrals_e_n = "None"
  touch io_ao_integrals_e_n 
  io_ao_integrals_kinetic = "None"
  touch io_ao_integrals_kinetic 

  call pot_hartree_trace
end

 subroutine give_v_h_n_at_r(r,v_exact,v_exp,v_pade)
 implicit none
 BEGIN_DOC
! input      : r(1) ==> r(1) = x, r(2) = y, r(3) = z
!
! input      : n_delta => number of application of the laplacian operator on the ao product 
!
! output     : delta_n_at_r ==> Delta^{n_delta}[ao(ao_num_local_1)ao(ao_num_local_2)] at r, where Delta is the laplacian operator 
!
!            : 
 END_DOC
 double precision, intent(in) :: r(3)
 double precision, intent(out) :: v_exact(n_states)
 double precision, intent(out) :: v_exp(0:order_derivative_ontop,n_states) 
 double precision, intent(out) :: v_pade(0:order_derivative_ontop,n_states)

 include 'utils/constants.include.F'
 double precision :: fact

 integer :: k,l,istate,n,der_loc_pad
 double precision :: NAI_pol_mult_erf_ao,n_order_delta 
 double precision :: coef_exp(0:order_derivative_ontop,n_states)

 integer :: m,INFO,i,j
 double precision :: denom,nume,x_loc
 double precision, allocatable :: matrix_alpha(:,:)
 double precision, allocatable :: vect_b_temp(:),vect_a(:),vect_b(:),c_pade(:)
 integer, allocatable :: IPIV(:)


 v_exact= 0.d0
 v_exp = 0.d0
 v_pade =0.d0
 coef_exp = 0.d0
 
 do istate = 1, n_states

  do l = 1, ao_num       
   do k = 1, ao_num        
    v_exact(istate) += ((one_e_dm_alpha_ao_for_dft(k,l,istate)+one_e_dm_beta_ao_for_dft(k,l,istate)) *(NAI_pol_mult_erf_ao(k,l,10000.d0,r) - NAI_pol_mult_erf_ao(k,l,mu_erf,r)))
   enddo
  enddo
!!!!!!! Expension classique

  do n = 0, order_derivative_ontop      

   do l = 1, ao_num       
    do k = 1, ao_num        
     call n_order_delta_ao_product(r,l,k,n,n_order_delta)
     coef_exp(n,istate) += (1.d0/(2.d0*dble(n)+1.d0))*((one_e_dm_alpha_ao_for_dft(k,l,istate)+one_e_dm_beta_ao_for_dft(k,l,istate)) * n_order_delta)
    enddo
   enddo
   coef_exp(n,istate) =  2.d0*sqrt(pi) * gamma((2.d0*dble(n)+3.d0)/2.d0) / ( fact(2*n) *(2.d0*dble(n)+2)) * coef_exp(n,istate)
   v_exp(n,istate)= coef_exp(n,istate)/(mu_erf**(2.d0*dble(n)+2.d0))
  enddo

!!!!!!!!!!!!!!!!!!!!! Pade
 !do der_loc_pad = 2, order_derivative_ontop 
  do der_loc_pad = 5,5 
   n=(der_loc_pad+1)/2  
   m=der_loc_pad+1-n
   print*,"n =",n
   print*,"m =",m
  !n=4
  !m=4
  !print*,"n =",n
  !print*,"m =",m
   allocate(matrix_alpha(m,m),vect_b_temp(m),IPIV(m),vect_a(0:n),vect_b(0:m+1),c_pade(0:n+m+1))

   matrix_alpha= 0.d0
   vect_b_temp=0.d0
   
   c_pade = 0.d0
   do i = 1,der_loc_pad+1
    c_pade(i)=coef_exp(i-1,istate)
   enddo 
   
   
   do j =1,m
    vect_b_temp(j) = - c_pade(n+j)
    do i = 1,m
     matrix_alpha(i,j)=c_pade(n+i-j)
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
  
   x_loc=1.d0/(mu_erf**2.d0) 
  
   nume = 0.d0
   do i= 0,n
    nume += vect_a(i)*x_loc**(dble(i))
   enddo
   
   denom = 0.d0
   do i= 0,m
    denom += vect_b(i)*x_loc**(dble(i))
   enddo
   v_pade(der_loc_pad,istate) = nume/denom

   deallocate(matrix_alpha,vect_b_temp,IPIV,vect_a,vect_b,c_pade) 
  enddo

 enddo  

 end

 subroutine pot_hartree_trace
 implicit none
 integer :: i,j,nr
 double precision :: r12, dr12,r12max,mu
 double precision :: r1(3),r2(3),v_exact(n_states),v_exp(0:order_derivative_ontop,n_states),v_pade(0:order_derivative_ontop,n_states)


 nr = 100 
 r12max = 2.d0
 dr12 = r12max/dble(nr)
 r1(1) = 0.0d0
 r1(2) = 0.0d0
 r1(3) = 0.0d0
 print*,'r1 = ',r1
 !!!!! tracer pot 

 do i = 1, nr
  r1(1) += dr12
  call give_v_h_n_at_r(r1,v_exact,v_exp,v_pade)
  write(33,'(10(F16.10,X))')r1(1),v_exact
  write(44,'(10(F16.10,X))')r1(1),v_exp
  write(55,'(10(F16.10,X))')r1(1),v_pade
 enddo


 end




