BEGIN_PROVIDER [double precision, ao_two_e_eff_dr12_pot_array, (ao_num,ao_num,ao_num,ao_num)]
 implicit none
 BEGIN_DOC
!   1 2                                1 2 
! < k l | [erf( mu r12) - 1] d/d_r12 | i j > on the AO basis
 END_DOC
 integer :: i,j,k,l,m
 double precision :: mu_in,d_dr12(3),d_dr12_large(3),accu
 mu_in = mu_erf
 !$OMP PARALLEL & 
 !$OMP DEFAULT(NONE) &
 !$OMP PRIVATE(i,j,k,l,m,d_dr12,d_dr12_large,accu) & 
 !$OMP SHARED(ao_two_e_eff_dr12_pot_array,mu_in,ao_num) 
 !$OMP DO SCHEDULE (dynamic)
 do j = 1, ao_num ! r2
  do i = 1, ao_num ! r1
   do l = 1, ao_num ! r2 
    do k = 1, ao_num ! r1 
      ! <kl|ij>
      !    the d/dr12 op acts on 1   2    
      !                            1   2
      call ao_two_e_eff_dr12_pot(i,k,j,l,mu_in,d_dr12,d_dr12_large)
      accu  = 0.d0
      do m = 1, 3
       accu += d_dr12(m) - d_dr12_large(m) 
      enddo
      !                           1 2 1 2 
      ao_two_e_eff_dr12_pot_array(k,l,i,j) = accu
    enddo
   enddo
  enddo
 enddo
 !$OMP END DO
 !$OMP END PARALLEL
END_PROVIDER 


BEGIN_PROVIDER [double precision, mo_two_e_eff_dr12_pot_array, (mo_num,mo_num,mo_num,mo_num)]
 implicit none
 BEGIN_DOC
! < k l | [erf( mu r12) - 1] d/d_r12 | i j > on the MO basis
 END_DOC
 integer :: i,j,k,l,m,n,p,q
 double precision, allocatable :: mo_tmp_1(:,:,:,:),mo_tmp_2(:,:,:,:),mo_tmp_3(:,:,:,:)
 mo_two_e_eff_dr12_pot_array = 0.d0
 allocate(mo_tmp_1(ao_num,ao_num,ao_num,mo_num))
 mo_tmp_1 = 0.d0
 do i = 1, mo_num ! r1
  do m = 1, ao_num
   do n = 1, ao_num
    do p = 1, ao_num
     do q = 1, ao_num
      ! <j p
      mo_tmp_1(p,n,m,i) += mo_coef(q,i) * ao_two_e_eff_dr12_pot_array(q,p,n,m)
!      double precision :: get_ao_two_e_integral
!      mo_tmp_1(p,n,m,i) += mo_coef(q,i) * get_ao_two_e_integral(q,p,n,m,ao_integrals_map)
     enddo
    enddo
   enddo
  enddo
 enddo
 FREE ao_two_e_eff_dr12_pot_array
 allocate(mo_tmp_2(ao_num,ao_num,mo_num,mo_num))
 mo_tmp_2 = 0.d0
 do i = 1, mo_num
  do j = 1, mo_num
    do m = 1, ao_num
     do n = 1, ao_num
      do p = 1, ao_num
       mo_tmp_2(n,m,j,i) += mo_tmp_1(p,n,m,i) * mo_coef(p,j)
      enddo
     enddo
    enddo
  enddo
 enddo
 deallocate(mo_tmp_1)
 allocate(mo_tmp_3(ao_num,mo_num,mo_num,mo_num))
 mo_tmp_3 = 0.d0
 do i = 1, mo_num
  do j = 1, mo_num
   do k = 1, mo_num
    do m = 1, ao_num
     do n = 1, ao_num
      mo_tmp_3(m,k,j,i) += mo_tmp_2(n,m,j,i) * mo_coef(n,k)
     enddo
    enddo
   enddo
  enddo
 enddo
 deallocate(mo_tmp_2)
 do i = 1, mo_num
  do j = 1, mo_num
   do k = 1, mo_num
    do l = 1, mo_num
     do m = 1, ao_num
      mo_two_e_eff_dr12_pot_array(l,k,j,i) += mo_coef(m,l) * mo_tmp_3(m,k,j,i)
     enddo
    enddo
   enddo
  enddo
 enddo
 deallocate(mo_tmp_3)

END_PROVIDER 
