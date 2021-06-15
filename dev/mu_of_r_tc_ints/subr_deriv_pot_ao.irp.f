
BEGIN_PROVIDER [ double precision, deriv_mu_r_pot_chemist_ao, (ao_num, ao_num, ao_num, ao_num)]
 implicit none
 BEGIN_DOC
! deriv_mu_r_pot_chemist_ao(i,k,j,l) = \int dr1 \int dr2  \phi_l(r2) \phi_k(r1) deriv_op [\phi_j(r2) \phi_i(r1)]
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
 call gamma_nabla_r1(deriv_mu_r_pot_chemist_ao)
 call non_hermit_r1(deriv_mu_r_pot_chemist_ao)
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

 allocate(a_mat(ao_num,ao_num,n_points_final_grid))

 do m = 1, 3
   do ipoint = 1, n_points_final_grid
    mu = mu_of_r_for_ints(ipoint,1)
    weight = final_weight_at_r_vector(ipoint)
    do k = 1, ao_num
     do i = 1, ao_num
      a_mat(i,k,ipoint) = ( aos_grad_in_r_array(k,ipoint,m) * aos_in_r_array(i,ipoint)   & 
                          + aos_grad_in_r_array(i,ipoint,m) * aos_in_r_array(k,ipoint) ) & 
                          * weight * cst * inv_2_mu_of_r_for_ints(ipoint,1)     * grad_mu_of_r_transp_for_ints(ipoint,1,m)
     enddo
    enddo
   enddo
   call dgemm("N","N",ao_num*ao_num,ao_num*ao_num,n_points_final_grid,1.d0,a_mat(1,1,1),ao_num*ao_num & 
                     ,gauss_ij_rk(1,1,1),n_points_final_grid,1.d0,big_mat,ao_num*ao_num)

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
 cst = -0.5d0 * inv_sq_pi 

 allocate(a_mat(ao_num,ao_num,n_points_final_grid))

 do m = 1, 3
   do ipoint = 1, n_points_final_grid
    mu = mu_of_r_for_ints(ipoint,1)
    weight = final_weight_at_r_vector(ipoint)
    do k = 1, ao_num
     do i = 1, ao_num
      a_mat(i,k,ipoint) =  aos_grad_in_r_array(i,ipoint,m) * aos_in_r_array(k,ipoint)   & 
                        * weight * cst * inv_2_mu_of_r_for_ints(ipoint,1) * grad_mu_of_r_transp_for_ints(ipoint,1,m)
     enddo
    enddo
   enddo
   call dgemm("N","N",ao_num*ao_num,ao_num*ao_num,n_points_final_grid,1.d0,a_mat(1,1,1),ao_num*ao_num & 
                     ,gauss_ij_rk(1,1,1),n_points_final_grid,1.d0,big_mat,ao_num*ao_num)

 enddo
  call wall_time(wall1)
  print*,'wall time for gamma_nabla_r1 ',wall1 - wall0   
end



subroutine non_hermit_r1(big_mat)
 implicit none
 BEGIN_DOC
! enters with the array big_mat and add the following integrals 
!
!  big_mat(i,k,j,l) += \int dr1 \int dr2 \phi_k(r1) \phi_l(r2) (erf (mu * r12) - 1) d/ dr12 \phi_i(r1) \phi_j(r2)
!
! WARNING : this is not hermitian and the nabla_1 acts on \phi_i(r1) 
 END_DOC
 include 'constants.include.F'
 double precision, intent(inout) :: big_mat(ao_num,ao_num,ao_num,ao_num)
 integer :: i,k,ipoint,m
 double precision :: weight,mu
 double precision, allocatable :: a_mat(:,:,:)
 double precision :: cst,r1(3)

 print*,'computing non_hermit_r1 ...'
 call wall_time(wall0)
 double precision :: wall0,wall1
 cst = 0.5d0

 allocate(a_mat(ao_num,ao_num,n_points_final_grid))

 do m = 1, 3
   ! x1 phi_k(r1) d/dx1 phi_i(r1) \int dr2 (erf(mu*r12) -1)/r12 phi_j(r2) phi_l(r2)
   do ipoint = 1, n_points_final_grid
    mu = mu_of_r_for_ints(ipoint,1)
    r1(:) = final_grid_points(:,ipoint)
    weight = final_weight_at_r_vector(ipoint)
    do k = 1, ao_num
     do i = 1, ao_num
      a_mat(i,k,ipoint) =  r1(m) * aos_grad_in_r_array(i,ipoint,m) * aos_in_r_array(k,ipoint) * weight * cst 
     enddo
    enddo
   enddo
   call dgemm("N","N",ao_num*ao_num,ao_num*ao_num,n_points_final_grid,1.d0,a_mat(1,1,1),ao_num*ao_num & 
                     ,v_ij_erf_rk_transp(1,1,1),n_points_final_grid,1.d0,big_mat,ao_num*ao_num)

   ! phi_k(r1) phi_i(r1) \int dr2 (erf(mu*r12) -1)/r12 x2 (d/dx2 phi_j(r2)) phi_l(r2)
   do ipoint = 1, n_points_final_grid
    mu = mu_of_r_for_ints(ipoint,1)
    weight = final_weight_at_r_vector(ipoint)
    do k = 1, ao_num
     do i = 1, ao_num
      a_mat(i,k,ipoint) =  aos_in_r_array(i,ipoint) * aos_in_r_array(k,ipoint) * weight * cst 
     enddo
    enddo
   enddo
   call dgemm("N","N",ao_num*ao_num,ao_num*ao_num,n_points_final_grid,1.d0,a_mat(1,1,1),ao_num*ao_num & 
                     ,x_d_dx_v_ij_erf_rk(1,1,1,m),n_points_final_grid,1.d0,big_mat,ao_num*ao_num)

   ! - x1 phi_k(r1) phi_i(r1) \int dr2 (erf(mu*r12) -1)/r12 (d/dx2 phi_j(r2)) phi_l(r2)
   do ipoint = 1, n_points_final_grid
    mu = mu_of_r_for_ints(ipoint,1)
    r1(:) = final_grid_points(:,ipoint)
    weight = final_weight_at_r_vector(ipoint)
    do k = 1, ao_num
     do i = 1, ao_num
      a_mat(i,k,ipoint) =  -r1(m) * aos_in_r_array(i,ipoint) * aos_in_r_array(k,ipoint) * weight * cst 
     enddo
    enddo
   enddo
   call dgemm("N","N",ao_num*ao_num,ao_num*ao_num,n_points_final_grid,1.d0,a_mat(1,1,1),ao_num*ao_num & 
                     ,d_dx_v_ij_erf_rk(1,1,1,m),n_points_final_grid,1.d0,big_mat,ao_num*ao_num)
   ! - d/d x1 phi_k(r1) phi_i(r1) \int dr2 (erf(mu*r12) -1)/r12 x2 phi_j(r2) phi_l(r2)
   do ipoint = 1, n_points_final_grid
    mu = mu_of_r_for_ints(ipoint,1)
    r1(:) = final_grid_points(:,ipoint)
    weight = final_weight_at_r_vector(ipoint)
    do k = 1, ao_num
     do i = 1, ao_num
      a_mat(i,k,ipoint) =  -aos_grad_in_r_array(i,ipoint,m) * aos_in_r_array(k,ipoint) * weight * cst 
     enddo
    enddo
   enddo
   call dgemm("N","N",ao_num*ao_num,ao_num*ao_num,n_points_final_grid,1.d0,a_mat(1,1,1),ao_num*ao_num & 
                     ,x_v_ij_erf_rk_transp_bis(1,1,1,m),n_points_final_grid,1.d0,big_mat,ao_num*ao_num)
 enddo
  call wall_time(wall1)
  print*,'wall time for non_hermit_r1 ',wall1 - wall0   
end



subroutine test_num_deriv_pot 
 implicit none
 include 'constants.include.F'
 integer :: ipoint,i,j,k,l,m,jpoint
 double precision :: r1(3),r2(3),weight1,weight2,weight_tot,ao_prod_r2,ao_prod_r1
 double precision :: lapl_gamma,cst_lapl_gamma,func_gamma_nabla_r1
 double precision :: thr, sq_thr,gauss_r12_mu_r1,erf_mur1, cst_gamma_nabla_r1 , func_non_hermit_at_r1
 double precision, allocatable :: accu(:,:,:,:)
 allocate(accu(ao_num, ao_num, ao_num, ao_num))
 thr = 1.d-15
 sq_thr = dsqrt(thr)
 cst_lapl_gamma = 0.25d0 * inv_sq_pi 
 cst_gamma_nabla_r1 = -1.d0 
 accu = 0.d0

 do jpoint = 1, n_points_final_grid ! r2
  weight2 = final_weight_at_r_vector(jpoint)
  do l = 1, ao_num
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
        accu(i,k,j,l) += func_gamma_nabla_r1(ipoint,jpoint,i,k,cst_gamma_nabla_r1) * weight1 * ao_prod_r2
        accu(i,k,j,l) += func_non_hermit_at_r1(ipoint,jpoint,i,k,j,l) * weight1 * weight2
                         
      enddo
     enddo
    enddo

   enddo
  enddo

 enddo

 double precision :: num_int,contrib,accu_naive,exact
 accu_naive = 0.d0
 do l = 1, ao_num
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
      print*,num_int, exact, contrib
     endif
     
    enddo
   enddo
  enddo
 enddo

 print*,'accu naive  = ',accu_naive/dble(ao_num**1)


end

subroutine test_big_array_ao_deriv
 implicit none
 integer :: i,j,k,l
 double precision :: num_int,contrib,accu_naive,exact
 accu_naive = 0.d0
 do l = 1, ao_num
  do j = 1, ao_num
   do k = 1, ao_num
    do i = 1, ao_num
     num_int = ao_two_e_eff_dr12_pot_array_new(k,i,l,j)
     exact = deriv_mu_r_pot_chemist_ao(i,k,j,l)
     contrib = dabs(num_int - exact )
     accu_naive += contrib
     if(contrib .gt. 1.d-10)then
      print*,''
      print*,'i,k,j,l',i,k,j,l
      print*,num_int, exact,contrib
     endif
     
    enddo
   enddo
  enddo
 enddo

 print*,'accu = ',accu_naive/dble(ao_num**4)


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

double precision function func_non_hermit_at_r1(ipoint,jpoint,i,k,j,l)
 BEGIN_DOC
! (erf(mu(r1)) - 1) d/dr12
 END_DOC
 implicit none 
 integer, intent(in) :: ipoint
 integer, intent(in) :: jpoint
 integer, intent(in) :: i,k,j,l
 double precision :: cst
 include 'constants.include.F'
 double precision :: r12,r12_sq,r12_vec(3),r1(3),r2(3),mu,ao_grad_prod(3,4)
 integer :: m,kk
 cst = 0.5d0
 func_non_hermit_at_r1 = 0.d0
 r1(:) = final_grid_points(:,ipoint)

 r2(:) = final_grid_points(:,jpoint)

 r12_vec(1) = (r1(1) - r2(1))
 r12_vec(2) = (r1(2) - r2(2))
 r12_vec(3) = (r1(3) - r2(3))
 r12_sq = r12_vec(1)*r12_vec(1) + r12_vec(2)*r12_vec(2) + r12_vec(3)*r12_vec(3)
 r12 = dsqrt(r12_sq)
 if(dabs(r12).lt.1.d-10)return
 ao_grad_prod = 0.d0
 do m = 1, 3
  ao_grad_prod(m,1) = r1(m) * aos_grad_in_r_array_transp_bis(ipoint,i,m) * aos_in_r_array_transp(ipoint,k) & 
                    * aos_in_r_array_transp(jpoint,j) * aos_in_r_array_transp(jpoint,l)
  ao_grad_prod(m,2) = r2(m) * aos_grad_in_r_array_transp_bis(jpoint,j,m) * aos_in_r_array_transp(jpoint,l) & 
                    * aos_in_r_array_transp(ipoint,i) * aos_in_r_array_transp(ipoint,k)
  ao_grad_prod(m,3) = -r1(m) * aos_in_r_array_transp(ipoint,i) * aos_in_r_array_transp(ipoint,k) & 
                    * aos_grad_in_r_array_transp_bis(jpoint,j,m) * aos_in_r_array_transp(jpoint,l)
  ao_grad_prod(m,4) = -r2(m) * aos_in_r_array_transp(jpoint,j) * aos_in_r_array_transp(jpoint,l) & 
                    * aos_grad_in_r_array_transp_bis(ipoint,i,m) * aos_in_r_array_transp(ipoint,k)
 enddo

 mu = mu_of_r_for_ints(ipoint,1)
 double precision :: func
 func = (derf(mu*r12) - 1.d0)/r12
 do kk = 1, 4
  do m = 1, 3
   func_non_hermit_at_r1 += cst * ao_grad_prod(m,kk) * func
  enddo
 enddo
end

double precision function gamma_r1(ipoint,jpoint,cst)
 implicit none
 integer, intent(in) :: ipoint, jpoint
 double precision, intent(in) :: cst

 double precision :: mu
 mu = mu_of_r_for_ints(ipoint,1)
 ! TODO : that function requires to evaluate mu(r1) and its gradient at any points ...


end

double precision function func_non_hermit_at_r1_bis(ipoint,jpoint,i,k,j,l,mu)
 BEGIN_DOC
! (erf(mu(r1)) - 1) d/dr12
 END_DOC
 implicit none 
 integer, intent(in) :: ipoint
 integer, intent(in) :: jpoint
 integer, intent(in) :: i,k,j,l
 double precision, intent(in) :: mu
 double precision :: cst
 include 'constants.include.F'
 double precision :: r12,r12_sq,r12_vec(3),r1(3),r2(3),ao_grad_prod(3,4)
 integer :: m,kk
 cst = 0.5d0
 func_non_hermit_at_r1_bis = 0.d0
 r1(:) = final_grid_points(:,ipoint)

 r2(:) = final_grid_points(:,jpoint)

 r12_vec(1) = (r1(1) - r2(1))
 r12_vec(2) = (r1(2) - r2(2))
 r12_vec(3) = (r1(3) - r2(3))
 r12_sq = r12_vec(1)*r12_vec(1) + r12_vec(2)*r12_vec(2) + r12_vec(3)*r12_vec(3)
 r12 = dsqrt(r12_sq)
 if(dabs(r12).lt.1.d-10)return
 ao_grad_prod = 0.d0
 do m = 1, 3
  ao_grad_prod(m,1) = r1(m) * aos_grad_in_r_array_transp_bis(ipoint,i,m) * aos_in_r_array_transp(ipoint,k) & 
                    * aos_in_r_array_transp(jpoint,j) * aos_in_r_array_transp(jpoint,l)
  ao_grad_prod(m,2) = r2(m) * aos_grad_in_r_array_transp_bis(jpoint,j,m) * aos_in_r_array_transp(jpoint,l) & 
                    * aos_in_r_array_transp(ipoint,i) * aos_in_r_array_transp(ipoint,k)
  ao_grad_prod(m,3) = -r1(m) * aos_in_r_array_transp(ipoint,i) * aos_in_r_array_transp(ipoint,k) & 
                    * aos_grad_in_r_array_transp_bis(jpoint,j,m) * aos_in_r_array_transp(jpoint,l)
  ao_grad_prod(m,4) = -r2(m) * aos_in_r_array_transp(jpoint,j) * aos_in_r_array_transp(jpoint,l) & 
                    * aos_grad_in_r_array_transp_bis(ipoint,i,m) * aos_in_r_array_transp(ipoint,k)
 enddo

 double precision :: func
 func = (derf(mu*r12) - 1.d0)/r12
 do kk = 1, 4
  do m = 1, 3
   func_non_hermit_at_r1_bis += cst * ao_grad_prod(m,kk) * func
  enddo
 enddo
end

