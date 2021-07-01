
BEGIN_PROVIDER [ double precision, deriv_mu_r_pot_chemist_mo, (mo_num, mo_num, mo_num, mo_num)]
 implicit none
 BEGIN_DOC
! deriv_mu_r_pot_chemist_mo(i,k,j,l) = \int dr1 \int dr2 \phi_l(r2) \phi_k(r1) deriv_op [\phi_j(r2) \phi_i(r1)]
!
! ALSO :: it is non hermitian !! BY convention the differential operators are applied to i(r1) ,j(r2)  
!
! NOTICE THAT : because of mu(r1) and the non symmetric form of the Jastrow factor, the integrals ARE NOT SYMMETRIC in r1, r2 
!
! deriv_mu_r_pot_chemist_mo(i,k,j,l) NOT EQUAL TO deriv_mu_r_pot_chemist_mo(j,l,i,k) for instance 
 END_DOC
 integer :: i,j,k,l,m,n,p,q
 double precision, allocatable :: mo_tmp_1(:,:,:,:),mo_tmp_2(:,:,:,:),mo_tmp_3(:,:,:,:)
 double precision :: wall0, wall1
 call wall_time(wall0)
 provide deriv_mu_r_pot_chemist_ao
 print*,'providing deriv_mu_r_pot_chemist_mo ...'
 allocate(mo_tmp_1(mo_num,ao_num,ao_num,ao_num))
 mo_tmp_1 = 0.d0
 do m = 1, ao_num
  do p = 1, ao_num
   do n = 1, ao_num
    do q = 1, ao_num
     do k = 1, mo_num                                                                                                                        
      !       (k n|p m)    = sum_q c_qk * (q n|p m)
      mo_tmp_1(k,n,p,m) += mo_coef_transp(k,q) * deriv_mu_r_pot_chemist_ao(q,n,p,m)
     enddo
    enddo
   enddo
  enddo
 enddo

 free deriv_mu_r_pot_chemist_ao
 allocate(mo_tmp_2(mo_num,mo_num,ao_num,ao_num))
 mo_tmp_2 = 0.d0
 do m = 1, ao_num
  do p = 1, ao_num
   do n = 1, ao_num
    do i = 1, mo_num
     do k = 1, mo_num
      !       (k i|p m) = sum_n c_ni * (k n|p m)
      mo_tmp_2(k,i,p,m) += mo_coef_transp(i,n) * mo_tmp_1(k,n,p,m)
     enddo
    enddo
   enddo
  enddo
 enddo
 deallocate(mo_tmp_1)
 allocate(mo_tmp_1(mo_num,mo_num,mo_num,ao_num))
 mo_tmp_1 = 0.d0
 do m = 1, ao_num
  do p = 1, ao_num
   do l = 1, mo_num
    do i = 1, mo_num
     do k = 1, mo_num
      mo_tmp_1(k,i,l,m) += mo_coef_transp(l,p) * mo_tmp_2(k,i,p,m)
     enddo
    enddo
   enddo
  enddo
 enddo
 deallocate(mo_tmp_2)
 deriv_mu_r_pot_chemist_mo = 0.d0
 do m = 1, ao_num
  do j = 1, mo_num
   do l = 1, mo_num
    do i = 1, mo_num
     do k = 1, mo_num
      deriv_mu_r_pot_chemist_mo(k,i,l,j) += mo_coef_transp(j,m)  * mo_tmp_1(k,i,l,m)
     enddo
    enddo
   enddo
  enddo
 enddo
 print*,'Time to perform 4 indx transformation = ',wall1 - wall0
 call wall_time(wall1)
END_PROVIDER 


BEGIN_PROVIDER [ double precision, deriv_mu_r_pot_physicist_mo, (mo_num, mo_num, mo_num, mo_num)]
 implicit none
 BEGIN_DOC
! deriv_mu_r_pot_physicist_mo(l,k,j,i) = <lk|ji> = \int dr1 dr2 \phi_l(r2) \phi_k(r1) \tilde{W}_ee(r1,r2) \phi_j(r2) \phi_i(r1)
!
! !!!! WARNING !!!! If constant_mu == .False. then a \mu(r1) is used and therefore IT IS NOT SYMMETRIC ANYMORE IN (r1, r2)
 END_DOC
 integer :: i,j,k,l
 
 if(constant_mu)then
  PROVIDE mo_two_e_eff_dr12_pot_array_physicist
  do i = 1, mo_num 
   do j = 1, mo_num 
    do k = 1, mo_num 
     do l = 1, mo_num 
      !                           2 1 2 1 
      deriv_mu_r_pot_physicist_mo(l,k,j,i) = mo_two_e_eff_dr12_pot_array_physicist(l,k,j,i)
     enddo
    enddo
   enddo
  enddo
 else
  PROVIDE deriv_mu_r_pot_chemist_mo
  do i = 1, mo_num 
   do j = 1, mo_num 
    do k = 1, mo_num 
     do l = 1, mo_num 
      !                         2 1 2 1 
      if(symmetrized_tc_h)then
       deriv_mu_r_pot_physicist_mo(l,k,j,i) = 0.5d0 * (deriv_mu_r_pot_chemist_mo(i,k,j,l) + deriv_mu_r_pot_chemist_mo(j,l,i,k))
      else
       deriv_mu_r_pot_physicist_mo(l,k,j,i) = deriv_mu_r_pot_chemist_mo(i,k,j,l)
      endif
     enddo
    enddo
   enddo
  enddo
  FREE deriv_mu_r_pot_chemist_mo 
 endif

END_PROVIDER 

subroutine test_num_deriv_pot_mo
 implicit none
 include 'constants.include.F'
 integer :: ipoint,i,j,k,l,m,jpoint
 double precision :: r1(3),r2(3),weight1,weight2,weight_tot,mo_prod_r2,mo_prod_r1
 double precision :: lapl_gamma_mo,cst_lapl_gamma,func_gamma_nabla_r1_mo
 double precision :: thr, sq_thr,gauss_r12_mu_r1,erf_mur1, cst_gamma_nabla_r1 , func_non_hermit_at_r1_mo, cst_non_hermit
 double precision, allocatable :: accu(:,:,:,:)
 allocate(accu(mo_num, mo_num, mo_num, mo_num))
 thr = 1.d-15
 sq_thr = dsqrt(thr)
 cst_lapl_gamma = 0.25d0 * inv_sq_pi 
 cst_gamma_nabla_r1 = -1.d0 
 cst_non_hermit = 0.5d0
 accu = 0.d0

 do jpoint = 1, n_points_final_grid ! r2
  weight2 = final_weight_at_r_vector(jpoint)
  do l = 1, mo_num
   do j = l, l
    mo_prod_r2 = mos_in_r_array(j,jpoint) * mos_in_r_array(l,jpoint) * weight2
    if(dabs(mo_prod_r2).lt.sq_thr)cycle
    do ipoint = 1, n_points_final_grid ! r1 
     weight1 = final_weight_at_r_vector(ipoint)
     weight_tot = weight1 * weight2
     if(dabs(weight_tot).lt.thr)cycle
     do k = l, l
      do i = 1, 2
       mo_prod_r1 = mos_in_r_array(i,ipoint) * mos_in_r_array(k,ipoint) * weight1
       if(dabs(mo_prod_r1).lt.sq_thr)cycle
        accu(i,k,j,l) += lapl_gamma_mo(ipoint,jpoint,i,k,cst_lapl_gamma) * weight1 * mo_prod_r2
        accu(i,k,j,l) += func_gamma_nabla_r1_mo(ipoint,jpoint,i,k,cst_gamma_nabla_r1) * weight1 * mo_prod_r2
        accu(i,k,j,l) += func_non_hermit_at_r1_mo(ipoint,jpoint,i,k,j,l) * weight1 * weight2
                         
      enddo
     enddo
    enddo

   enddo
  enddo

 enddo

 double precision :: num_int,contrib,accu_naive,exact
 accu_naive = 0.d0
 do l = 1, mo_num
  do j = l, l
   do k = l, l
    do i = 1, 2
     num_int = accu(i,k,j,l)
     exact = deriv_mu_r_pot_chemist_mo(i,k,j,l)
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
 print*,'accu naive  = ',accu_naive/dble(mo_num**1)
end

subroutine test_big_array_mo_deriv
 implicit none
 integer :: i,j,k,l
 double precision :: num_int,contrib,accu_naive,exact
 double precision :: accu_n,accu_relat,contrib_relat
 accu_naive = 0.d0
 accu_relat = 0.d0
 accu_n = 0.d0
 print*,''
 print*,''
 print*,'testing the mixed numerical with the analytical for the non hermitian part'
 print*,''
 do l = 1, mo_num
  do j = 1, mo_num
   do k = 1, mo_num
    do i = 1, mo_num
!     num_int = mo_two_e_eff_dr12_pot_array_new_3(k,i,l,j)
!     num_int = mo_two_e_eff_dr12_pot_array_new_bis(k,i,l,j)
     num_int = deriv_mu_r_pot_chemist_mo(i,k,j,l)
     exact = mo_two_e_eff_dr12_pot_array_chemist(k,i,l,j)
     contrib = dabs(num_int - exact )
     if(dabs(exact).gt.1.d-12)then
      accu_naive += contrib
      accu_n += 1.d0
      contrib_relat = contrib/dabs(exact)
      accu_relat += contrib_relat
      if(contrib_relat .gt. 1.d-3)then
       print*,''
       print*,'i,k,j,l',i,k,j,l
       print*,num_int, exact
       print*,contrib, contrib_relat
      endif
     endif
     
    enddo
   enddo
  enddo
 enddo

 print*,'accu       = ',accu_naive/accu_n
 print*,'accu relat = ',accu_relat/accu_n


end

double precision function lapl_gamma_mo(ipoint,jpoint,i,k,cst)
 implicit none 
 integer, intent(in) :: ipoint
 integer, intent(in) :: jpoint
 integer, intent(in) :: i,k
 double precision, intent(in) :: cst
 include 'constants.include.F'
 double precision :: r12,r12_sq,r12_vec(3),r1(3),r2(3),mu,mo_grad_prod(3)
 integer :: m

 lapl_gamma_mo = 0.d0
 r1(:) = final_grid_points(:,ipoint)

 r2(:) = final_grid_points(:,jpoint)

 r12_vec(1) = (r1(1) - r2(1))
 r12_vec(2) = (r1(2) - r2(2))
 r12_vec(3) = (r1(3) - r2(3))
 r12_sq = r12_vec(1)*r12_vec(1) + r12_vec(2)*r12_vec(2) + r12_vec(3)*r12_vec(3)
 r12 = dsqrt(r12_sq)
 do m = 1, 3
  mo_grad_prod(m) = mos_grad_in_r_array_transp_bis(ipoint,k,m) * mos_in_r_array_transp(ipoint,i) & 
                  + mos_grad_in_r_array_transp_bis(ipoint,i,m) * mos_in_r_array_transp(ipoint,k) 
 enddo

 mu = mu_of_r_for_ints(ipoint,1)
 lapl_gamma_mo = cst * inv_2_mu_of_r_for_ints(ipoint,1) * dexp(- mu * mu * r12_sq) & 
             * ( grad_mu_of_r_for_ints(1,ipoint,1) * mo_grad_prod(1) & 
               + grad_mu_of_r_for_ints(2,ipoint,1) * mo_grad_prod(2) &
               + grad_mu_of_r_for_ints(3,ipoint,1) * mo_grad_prod(3) )   
end

double precision function func_gamma_nabla_r1_mo(ipoint,jpoint,i,k,cst)
 implicit none 
 integer, intent(in) :: ipoint
 integer, intent(in) :: jpoint
 integer, intent(in) :: i,k
 double precision, intent(in) :: cst
 include 'constants.include.F'
 double precision :: r12,r12_sq,r12_vec(3),r1(3),r2(3),mu,mo_grad_prod(3)
 integer :: m

 func_gamma_nabla_r1_mo = 0.d0
 r1(:) = final_grid_points(:,ipoint)

 r2(:) = final_grid_points(:,jpoint)

 r12_vec(1) = (r1(1) - r2(1))
 r12_vec(2) = (r1(2) - r2(2))
 r12_vec(3) = (r1(3) - r2(3))
 r12_sq = r12_vec(1)*r12_vec(1) + r12_vec(2)*r12_vec(2) + r12_vec(3)*r12_vec(3)
 r12 = dsqrt(r12_sq)
 do m = 1, 3
  mo_grad_prod(m) = mos_grad_in_r_array_transp_bis(ipoint,i,m) * mos_in_r_array_transp(ipoint,k) 
 enddo

 mu = mu_of_r_for_ints(ipoint,1)
 func_gamma_nabla_r1_mo = cst * inv_2_mu_of_r_for_ints(ipoint,1) * dexp(- mu * mu * r12_sq) & 
             * ( grad_mu_of_r_for_ints(1,ipoint,1) * mo_grad_prod(1) & 
               + grad_mu_of_r_for_ints(2,ipoint,1) * mo_grad_prod(2) &
               + grad_mu_of_r_for_ints(3,ipoint,1) * mo_grad_prod(3) )   
end

!double precision function func_non_hermit_at_r1_mo(ipoint,jpoint,i,k,j,l,cst)
! implicit none 
! integer, intent(in) :: ipoint
! integer, intent(in) :: jpoint
! integer, intent(in) :: i,k,j,l
! double precision, intent(in) :: cst
! include 'constants.include.F'
! double precision :: r12,r12_sq,r12_vec(3),r1(3),r2(3),mu,mo_grad_prod(3,4)
! integer :: m,kk
!
! func_non_hermit_at_r1_mo = 0.d0
! r1(:) = final_grid_points(:,ipoint)
!
! r2(:) = final_grid_points(:,jpoint)
!
! r12_vec(1) = (r1(1) - r2(1))
! r12_vec(2) = (r1(2) - r2(2))
! r12_vec(3) = (r1(3) - r2(3))
! r12_sq = r12_vec(1)*r12_vec(1) + r12_vec(2)*r12_vec(2) + r12_vec(3)*r12_vec(3)
! r12 = dsqrt(r12_sq)
! if(dabs(r12).lt.1.d-10)return
! mo_grad_prod = 0.d0
! do m = 1, 3
!  mo_grad_prod(m,1) = r1(m) * mos_grad_in_r_array_transp_bis(ipoint,i,m) * mos_in_r_array_transp(ipoint,k) & 
!                    * mos_in_r_array_transp(jpoint,j) * mos_in_r_array_transp(jpoint,l)
!  mo_grad_prod(m,2) = r2(m) * mos_grad_in_r_array_transp_bis(jpoint,j,m) * mos_in_r_array_transp(jpoint,l) & 
!                    * mos_in_r_array_transp(ipoint,i) * mos_in_r_array_transp(ipoint,k)
!  mo_grad_prod(m,3) = -r1(m) * mos_in_r_array_transp(ipoint,i) * mos_in_r_array_transp(ipoint,k) & 
!                    * mos_grad_in_r_array_transp_bis(jpoint,j,m) * mos_in_r_array_transp(jpoint,l)
!  mo_grad_prod(m,4) = -r2(m) * mos_in_r_array_transp(jpoint,j) * mos_in_r_array_transp(jpoint,l) & 
!                    * mos_grad_in_r_array_transp_bis(ipoint,i,m) * mos_in_r_array_transp(ipoint,k)
! enddo
!
! mu = mu_of_r_for_ints(ipoint,1)
! double precision :: func
! func = (derf(mu*r12) - 1.d0)/r12
! do kk = 1, 4
!  do m = 1, 3
!   func_non_hermit_at_r1_mo += cst * mo_grad_prod(m,kk) * func
!  enddo
! enddo
!end
!
