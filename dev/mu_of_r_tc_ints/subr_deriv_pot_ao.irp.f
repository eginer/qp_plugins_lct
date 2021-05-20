
BEGIN_PROVIDER [ double precision, deriv_mu_r_pot_chemist_ao, (ao_num, ao_num, ao_num, ao_num)]
 implicit none
 BEGIN_DOC
! deriv_mu_r_pot_chemist_ao(i,k,j,l) = \int dr1 \int dr2 \phi_i(r1) \phi_k(r1) deriv_op \phi_j(r2) \phi_l(r2)
!
! ALSO :: it is non hermitian !! BY convention the differential operators are applied to i(r1) ,j(r2)  
!
! NOTICE THAT : because of mu(r1) and the non symmetric form of the Jastrow factor, the integrals ARE NOT SYMMETRIC in r1, r2 
!
! deriv_mu_r_pot_chemist_ao(i,k,j,l) NOT EQUAL TO deriv_mu_r_pot_chemist_ao(j,l,i,k) for instance 
!
 END_DOC

 print*,''
 print*,'providing deriv_mu_r_pot_chemist_ao ...'
 print*,''
 print*,''
 call wall_time(wall0)
 double precision :: wall0,wall1


 deriv_mu_r_pot_chemist_ao = 0.d0

 call lapl_gamm_r1(deriv_mu_r_pot_chemist_ao)
! call gamma_nabla_r1(deriv_mu_r_pot_chemist_ao)

 call wall_time(wall1)
 print*,''
 print*,''
 print*,'wall time for deriv_mu_r_pot_chemist_ao ', wall1 - wall0   
 print*,''
 print*,''
END_PROVIDER 



subroutine lapl_gamm_r1(big_mat)
 implicit none
 BEGIN_DOC
! enters with the array big_mat and add the following integrals 
!
!  big_mat(i,k,j,l) += - 1/2 \int dr1 \int dr2 \phi_i(r1) \phi_k(r1) \int dr2 (\nabla_1 \gamma_r1 ) \phi_j(r2) \phi_l(r2)
!
!  =  1/2 \int dr1 (\phi_i(r1) [d/dx1 \phi_k(r1)] + \phi_k(r1) [d/dx1 \phi_k(r1)]) . d/dx1 \gamma_r1 \int dr2 \phi_j(r2) \phi_l(r2)
!
! with gamma_r1(r1,r2) = 1/(2 \sqpi (\mu(r1))^2) e^{-(\mu(r1) r12)^2} \nabla_1 \mu(r1)
 END_DOC
 include 'constants.include.F'
 double precision, intent(inout) :: big_mat(ao_num,ao_num,ao_num,ao_num)
 integer :: i,k,ipoint,m
 double precision :: weight,mu
 double precision, allocatable :: a_mat(:,:,:)
 double precision :: cst

 print*,'computing lapl_gamm_r1 ...'
 call wall_time(wall0)
 double precision :: wall0,wall1
 cst = 0.25d0 * inv_sq_pi 

 allocate(a_mat(n_points_final_grid,ao_num,ao_num))

 do m = 1, 3
   do i = 1, ao_num
    do k = 1, ao_num
     do ipoint = 1, n_points_final_grid
      mu = mu_of_r_for_ints(ipoint,1)
      weight = final_weight_at_r_vector(ipoint)
      a_mat(ipoint,k,i) = ( aos_grad_in_r_array_transp_bis(ipoint,k,m) * aos_in_r_array_transp(ipoint,i)   & 
                          + aos_grad_in_r_array_transp_bis(ipoint,i,m) * aos_in_r_array_transp(ipoint,k) ) & 
                          * weight * cst * inv_2_mu_of_r_for_ints(ipoint,1)     * grad_mu_of_r_transp_for_ints(ipoint,1,m)
     enddo
    enddo
   enddo
!   call dgemm("N","N",ao_num*ao_num,ao_num*ao_num,n_points_final_grid,1.d0,gauss_ij_rk(1,1,1),ao_num*ao_num & 
!                     ,a_mat(1,1,1),n_points_final_grid,1.d0,big_mat,ao_num*ao_num)
   call dgemm("T","T",ao_num*ao_num,ao_num*ao_num,n_points_final_grid,1.d0,a_mat(1,1,1),n_points_final_grid & 
                     ,gauss_ij_rk(1,1,1),ao_num*ao_num,1.d0,big_mat,ao_num*ao_num)

 enddo
  call wall_time(wall1)
  print*,'wall time for lapl_gamm_r1 ',wall1 - wall0   
end


subroutine gamma_nabla_r1(big_mat)
 implicit none
 BEGIN_DOC
! enters with the array big_mat and add the following integrals 
!
!  big_mat(i,k,j,l) += - \int dr1 \int dr2 \phi_j(r2) \phi_l(r2) \phi_k(r1) (\nabla_1 \gamma_r1 ) . \nabla_1 \phi_i(r1)
!
! with gamma_r1(r1,r2) = 1/(2 \sqpi (\mu(r1))^2) e^{-(\mu(r1) r12)^2} \nabla_1 \mu(r1)
!
! WARNING : this is not hermitian and the nabla_1 acts on \phi_i(r1) 
 END_DOC
 include 'constants.include.F'
 double precision, intent(inout) :: big_mat(ao_num,ao_num,ao_num,ao_num)
 integer :: i,k,ipoint,m
 double precision :: weight,mu
 double precision, allocatable :: a_mat(:,:,:)
 double precision :: cst

 print*,'computing gamma_nabla_r1 ...'
 call wall_time(wall0)
 double precision :: wall0,wall1
 cst = -1.d0 

 allocate(a_mat(n_points_final_grid,ao_num,ao_num))

 do m = 1, 3
   do k = 1, ao_num
    do i = 1, ao_num
     do ipoint = 1, n_points_final_grid
      mu = mu_of_r_for_ints(ipoint,1)
      weight = final_weight_at_r_vector(ipoint)
      a_mat(ipoint,i,k) =  aos_grad_in_r_array_transp_bis(ipoint,i,m) * aos_in_r_array_transp(ipoint,k)   & 
                        * weight * cst * inv_2_mu_of_r_for_ints(ipoint,1) * grad_mu_of_r_transp_for_ints(ipoint,1,m)
     enddo
    enddo
   enddo
   call dgemm("T","T",ao_num*ao_num,ao_num*ao_num,n_points_final_grid,1.d0,a_mat(1,1,1),n_points_final_grid & 
                     ,gauss_ij_rk(1,1,1),ao_num*ao_num,1.d0,big_mat,ao_num*ao_num)

 enddo
  call wall_time(wall1)
  print*,'wall time for gamma_nabla_r1 ',wall1 - wall0   
end


subroutine test_lapl_gamm_r1
 implicit none
 include 'constants.include.F'
 integer :: ipoint,i,j,k,l,m,jpoint
 double precision :: r1(3),r2(3),weight1,weight2,weight_tot,ao_prod_r2,ao_prod_r1
 double precision :: lapl_gamma,cst_lapl_gamma,func_gamma_nabla_r1
 double precision :: thr, sq_thr,gauss_r12_mu_r1,erf_mur1, cst_gamma_nabla_r1 
 double precision, allocatable :: accu(:,:,:,:)
 allocate(accu(ao_num, ao_num, ao_num, ao_num))
 thr = 1.d-15
 sq_thr = dsqrt(thr)
 cst_lapl_gamma = 0.25d0 * inv_sq_pi 
 cst_gamma_nabla_r1 = -1.d0 
 accu = 0.d0

 do jpoint = 1, n_points_final_grid ! r2
  weight2 = final_weight_at_r_vector(jpoint)
  do l = 1, 1
   do j = l, l
    ao_prod_r2 = aos_in_r_array(j,jpoint) * aos_in_r_array(l,jpoint) * weight2
    if(dabs(ao_prod_r2).lt.sq_thr)cycle
    do ipoint = 1, n_points_final_grid ! r1 
     weight1 = final_weight_at_r_vector(ipoint)
     weight_tot = weight1 * weight2
     if(dabs(weight_tot).lt.thr)cycle
     do k = l, l
      do i = 1, 2
       ao_prod_r1 = aos_in_r_array(i,ipoint) * aos_in_r_array(k,ipoint) * weight1
       if(dabs(ao_prod_r1).lt.sq_thr)cycle
        accu(i,k,j,l) += lapl_gamma(ipoint,jpoint,i,k,cst_lapl_gamma) * weight1 * ao_prod_r2
!        accu(i,k,j,l) += func_gamma_nabla_r1(ipoint,jpoint,i,k,cst_gamma_nabla_r1) * weight1 * ao_prod_r2
                         
      enddo
     enddo
    enddo

   enddo
  enddo

 enddo

 double precision :: num_int,contrib,accu_naive,exact
 accu_naive = 0.d0
 do l = 1, 1
  do j = l, l
   do k = l, l
    do i = 1, 2
     num_int = accu(i,k,j,l)
     exact = deriv_mu_r_pot_chemist_ao(i,k,j,l)
     contrib = dabs(num_int - exact )
     accu_naive += contrib
     if(contrib .gt. 1.d-10)then
      print*,''
      print*,'i,k,j,l',i,k,j,l
      print*,contrib,num_int, exact
     endif
     
    enddo
   enddo
  enddo
 enddo

 print*,'accu naive  = ',accu_naive/dble(ao_num**1)


end

double precision function lapl_gamma(ipoint,jpoint,i,k,cst)
 implicit none 
 integer, intent(in) :: ipoint
 integer, intent(in) :: jpoint
 integer, intent(in) :: i,k
 double precision, intent(in) :: cst
 include 'constants.include.F'
 double precision :: r12,r12_sq,r12_vec(3),r1(3),r2(3),mu,ao_grad_prod(3)
 integer :: m

 lapl_gamma = 0.d0
 r1(:) = final_grid_points(:,ipoint)

 r2(:) = final_grid_points(:,jpoint)

 r12_vec(1) = (r1(1) - r2(1))
 r12_vec(2) = (r1(2) - r2(2))
 r12_vec(3) = (r1(3) - r2(3))
 r12_sq = r12_vec(1)*r12_vec(1) + r12_vec(2)*r12_vec(2) + r12_vec(3)*r12_vec(3)
 r12 = dsqrt(r12_sq)
 do m = 1, 3
  ao_grad_prod(m) = aos_grad_in_r_array_transp_bis(ipoint,k,m) * aos_in_r_array_transp(ipoint,i) & 
                  + aos_grad_in_r_array_transp_bis(ipoint,i,m) * aos_in_r_array_transp(ipoint,k) 
 enddo

 mu = mu_of_r_for_ints(ipoint,1)
 lapl_gamma = cst * inv_2_mu_of_r_for_ints(ipoint,1) * dexp(- mu * mu * r12_sq) & 
             * ( grad_mu_of_r_for_ints(1,ipoint,1) * ao_grad_prod(1) & 
               + grad_mu_of_r_for_ints(2,ipoint,1) * ao_grad_prod(2) &
               + grad_mu_of_r_for_ints(3,ipoint,1) * ao_grad_prod(3) )   
end

double precision function func_gamma_nabla_r1(ipoint,jpoint,i,k,cst)
 implicit none 
 integer, intent(in) :: ipoint
 integer, intent(in) :: jpoint
 integer, intent(in) :: i,k
 double precision, intent(in) :: cst
 include 'constants.include.F'
 double precision :: r12,r12_sq,r12_vec(3),r1(3),r2(3),mu,ao_grad_prod(3)
 integer :: m

 func_gamma_nabla_r1 = 0.d0
 r1(:) = final_grid_points(:,ipoint)

 r2(:) = final_grid_points(:,jpoint)

 r12_vec(1) = (r1(1) - r2(1))
 r12_vec(2) = (r1(2) - r2(2))
 r12_vec(3) = (r1(3) - r2(3))
 r12_sq = r12_vec(1)*r12_vec(1) + r12_vec(2)*r12_vec(2) + r12_vec(3)*r12_vec(3)
 r12 = dsqrt(r12_sq)
 do m = 1, 3
  ao_grad_prod(m) = aos_grad_in_r_array_transp_bis(ipoint,i,m) * aos_in_r_array_transp(ipoint,k) 
 enddo

 mu = mu_of_r_for_ints(ipoint,1)
 func_gamma_nabla_r1 = cst * inv_2_mu_of_r_for_ints(ipoint,1) * dexp(- mu * mu * r12_sq) & 
             * ( grad_mu_of_r_for_ints(1,ipoint,1) * ao_grad_prod(1) & 
               + grad_mu_of_r_for_ints(2,ipoint,1) * ao_grad_prod(2) &
               + grad_mu_of_r_for_ints(3,ipoint,1) * ao_grad_prod(3) )   
end
