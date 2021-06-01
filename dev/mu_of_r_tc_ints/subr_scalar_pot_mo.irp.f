BEGIN_PROVIDER [ double precision, scalar_mu_r_pot_chemist_mo, (mo_num, mo_num, mo_num, mo_num)]
 implicit none
 BEGIN_DOC
! scalar_mu_r_pot_chemist_mo(i,k,j,l) = \int dr1 \int dr2 \phi_i(r1) \phi_k(r1) W_ee(r12,\mu(r1)) \phi_j(r2) \phi_l(r2)
!
! NOTICE THAT : because of mu(r1) and the non symmetric form of the Jastrow factor, the integrals ARE NOT SYMMETRIC in r1, r2 
!
! scalar_mu_r_pot_chemist_mo(i,k,j,l) NOT EQUAL TO scalar_mu_r_pot_chemist_mo(j,l,i,k) for instance 
 END_DOC
 integer :: i,j,k,l,m,n,p,q
 double precision, allocatable :: mo_tmp_1(:,:,:,:),mo_tmp_2(:,:,:,:),mo_tmp_3(:,:,:,:)

 allocate(mo_tmp_1(mo_num,ao_num,ao_num,ao_num))
 mo_tmp_1 = 0.d0
 do m = 1, ao_num
  do p = 1, ao_num
   do n = 1, ao_num
    do q = 1, ao_num
     do k = 1, mo_num                                                                                                                        
      !       (k n|p m)    = sum_q c_qk * (q n|p m)
      mo_tmp_1(k,n,p,m) += mo_coef_transp(k,q) * scalar_mu_r_pot_chemist_ao(q,n,p,m)
     enddo
    enddo
   enddo
  enddo
 enddo

 free scalar_mu_r_pot_chemist_ao
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
 scalar_mu_r_pot_chemist_mo = 0.d0
 do m = 1, ao_num
  do j = 1, mo_num
   do l = 1, mo_num
    do i = 1, mo_num
     do k = 1, mo_num
      scalar_mu_r_pot_chemist_mo(k,i,l,j) += mo_coef_transp(j,m)  * mo_tmp_1(k,i,l,m)
     enddo
    enddo
   enddo
  enddo
 enddo

END_PROVIDER 

BEGIN_PROVIDER [ double precision, scalar_mu_r_pot_physicist_mo, (mo_num, mo_num, mo_num, mo_num)]
 implicit none
 BEGIN_DOC
! scalar_mu_r_pot_physicist_mo(l,k,j,i) = <lk|ji> = \int dr1 dr2 \phi_l(r2) \phi_k(r1) \tilde{W}_ee(r1,r2) \phi_j(r2) \phi_i(r1)
!
! !!!! WARNING !!!! If constant_mu == .False. then a \mu(r1) is used and therefore IT IS NOT SYMMETRIC ANYMORE IN (r1, r2)
 END_DOC
 integer :: i,j,k,l
 double precision :: get_mo_two_e_integral_erf,mo_two_e_integral_eff_pot
 if(constant_mu)then
  PROVIDE mo_two_e_integrals_in_map mo_integrals_map big_array_exchange_integrals 
  PROVIDE mo_two_e_integrals_eff_pot_in_map mo_two_e_integrals_erf_in_map
  do i = 1, mo_num 
   do j = 1, mo_num 
    do k = 1, mo_num 
     do l = 1, mo_num 
      !                            2 1 2 1 
      scalar_mu_r_pot_physicist_mo(l,k,j,i) = get_mo_two_e_integral_erf(l,k,j,i,mo_integrals_erf_map) & 
                                            + mo_two_e_integral_eff_pot(l,k,j,i)
     enddo
    enddo
   enddo
  enddo
 else
  do i = 1, mo_num 
   do j = 1, mo_num 
    do k = 1, mo_num 
     do l = 1, mo_num 
      !                          2 1 2 1 
      scalar_mu_r_pot_physicist_mo(l,k,j,i) = scalar_mu_r_pot_chemist_mo(i,k,j,l)
     enddo
    enddo
   enddo
  enddo
  FREE scalar_mu_r_pot_chemist_mo 
 endif

END_PROVIDER 


subroutine test_num_scal_pot_mo
 implicit none
 include 'constants.include.F'
 integer :: ipoint,i,j,k,l,m,jpoint
 double precision :: r1(3),r2(3),weight1,weight2,weight_tot
 double precision :: nabla_sq_term,cst_nabla,ao_prod_r1,ao_prod_r2,nabla_r12_1,erf_mu_sq,nabla_r12_2
 double precision :: thr, sq_thr,cst_nabla_r12_1,cst_gauss_r12,gauss_r12_mu_r1,erf_mur1,cst_nabla_r12_2
 double precision, allocatable :: accu(:,:,:,:)
 allocate(accu(mo_num, mo_num, mo_num, mo_num))
 thr = 1.d-15
 sq_thr = dsqrt(thr)
 cst_nabla = -1.d0/(8.d0 * pi)
 cst_nabla_r12_1 = 0.5d0 * inv_sq_pi
 cst_nabla_r12_2 = -0.25d0 * inv_sq_pi
 cst_gauss_r12 = inv_sq_pi
 accu = 0.d0

 do jpoint = 1, n_points_final_grid ! r2
  weight2 = final_weight_at_r_vector(jpoint)
  do l = 1, mo_num
!  do l = 1, 1
!   do j = 1, mo_num
   do j = l, l
    ao_prod_r2 = mos_in_r_array(j,jpoint) * mos_in_r_array(l,jpoint) * weight2
    if(dabs(ao_prod_r2).lt.sq_thr)cycle
    do ipoint = 1, n_points_final_grid ! r1 
     weight1 = final_weight_at_r_vector(ipoint)
     weight_tot = weight1 * weight2
     if(dabs(weight_tot).lt.thr)cycle
!     do i = 1, mo_num
!      do k = 1, mo_num
     do k = l, l
      do i = j, j
       ao_prod_r1 = mos_in_r_array(i,ipoint) * mos_in_r_array(k,ipoint) * weight1
       if(dabs(ao_prod_r1).lt.sq_thr)cycle
        accu(i,k,j,l) += erf_mur1(ipoint,jpoint) * ao_prod_r1 * ao_prod_r2
        accu(i,k,j,l) += gauss_r12_mu_r1(ipoint,jpoint,cst_gauss_r12) * ao_prod_r1 * ao_prod_r2
        accu(i,k,j,l) += erf_mu_sq(ipoint,jpoint) * ao_prod_r1 * ao_prod_r2
        accu(i,k,j,l) += nabla_sq_term(ipoint,jpoint,cst_nabla) * ao_prod_r1 * ao_prod_r2
        accu(i,k,j,l) += nabla_r12_1(ipoint,jpoint,cst_nabla_r12_1) * ao_prod_r1 * ao_prod_r2
        accu(i,k,j,l) += nabla_r12_2(ipoint,jpoint,cst_nabla_r12_2) * ao_prod_r1 * ao_prod_r2
      enddo
     enddo
    enddo

   enddo
  enddo

 enddo

 double precision :: num_int,contrib,accu_naive
 accu_naive = 0.d0
! do l = 1, 1
 do l = 1, mo_num
!  do j = 1, mo_num
  do j = l, l
   do k = l, l
    do i = j, j
     num_int = accu(i,k,j,l)
     contrib = dabs(num_int - scalar_mu_r_pot_chemist_mo(i,k,j,l) )
     accu_naive += contrib
     if(contrib .gt. 1.d-10)then
      print*,''
      print*,'i,k,j,l',i,k,j,l
      print*,contrib,num_int, scalar_mu_r_pot_chemist_mo(i,k,j,l)
     endif
     
    enddo
   enddo
  enddo
 enddo

 print*,'accu naive  = ',accu_naive/dble(mo_num**1)
! print*,'accu dgemm  = ',accu2/dble(mo_num**4)
end

subroutine test_big_array_mo_scal
 implicit none
 integer :: i,j,k,l
 double precision :: num_int,contrib,accu_naive,exact,get_mo_two_e_integral_erf,mo_two_e_integral_eff_pot
 double precision :: accu_n,accu_relat
 accu_n = 0.d0
 accu_relat = 0.d0
 accu_naive = 0.d0
 print*,''
 print*,''
 print*,'testing the mixed numerical with the analytical for the hermitian part'
 print*,''
 do l = 1, mo_num
  do j = 1, mo_num
   do k = 1, mo_num
    do i = 1, mo_num
!     num_int = mo_two_e_eff_dr12_pot_array_new_3(k,i,l,j)
!     num_int = mo_two_e_eff_dr12_pot_array_new_bis(k,i,l,j)
     exact = get_mo_two_e_integral_erf(i,j,k,l,mo_integrals_erf_map)
     exact += mo_two_e_integral_eff_pot(i,j,k,l)
     num_int = scalar_mu_r_pot_chemist_mo(i,k,j,l)
     contrib = dabs(num_int - exact )
     if(dabs(exact).gt.1.d-15)then
      accu_naive += contrib
      accu_n += 1.d0
      accu_relat += contrib/dabs(exact)
     endif
     if(contrib .gt. 1.d-10)then
      print*,''
      print*,'i,k,j,l',i,k,j,l
      print*,'num, exact, delta '
      print*,num_int, exact,contrib
     endif
     
    enddo
   enddo
  enddo
 enddo

 print*,'accu       = ',accu_naive/accu_n**4
 print*,'accu relat = ',accu_relat/accu_n**4


end

