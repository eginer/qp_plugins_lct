BEGIN_PROVIDER [double precision, ao_two_e_eff_dr12_pot_array, (ao_num,ao_num,ao_num,ao_num)]
 implicit none
 BEGIN_DOC
                                       !   1 2                                1 2 
! ao_two_e_eff_dr12_pot_array(k,l,i,j) = < k l | [erf( mu r12) - 1] d/d_r12 | i j > on the AO basis
!
! WARNING !!!!!! ONLY WORKS FOR ATOMIC SYSTEMS, NOT MOLECULAR SYSTEMS 
 END_DOC
 integer :: i,j,k,l,m
 double precision :: mu_in,d_dr12(3),d_dr12_large(3),accu
 double precision, allocatable :: ao_ints(:)
 mu_in = mu_erf
 double precision :: wall0,wall1
 call wall_time(wall0)
 print*,'Providing ao_two_e_eff_dr12_pot_array ...'
 PROVIDE ao_two_e_integrals_in_map ao_integrals_map
 !$OMP PARALLEL & 
 !$OMP DEFAULT(NONE) &
 !$OMP PRIVATE(i,j,k,l,m,d_dr12,d_dr12_large,accu, ao_ints) & 
 !$OMP SHARED(ao_two_e_eff_dr12_pot_array,mu_in,ao_num,ao_integrals_threshold)   
 allocate(ao_ints(ao_num))
 !$OMP DO SCHEDULE (dynamic)
 do j = 1, ao_num ! r2
  do i = 1, ao_num ! r1
   do l = 1, ao_num ! r2 
     call get_ao_two_e_integrals(i,j,l,ao_num,ao_ints)
    do k = 1, ao_num ! r1 
      if(dabs(ao_ints(k)).lt.ao_integrals_threshold)then
       ao_two_e_eff_dr12_pot_array(k,l,i,j)  = 0.d0
       cycle
      endif
      ! <kl|ij>
      !    the d/dr12 op acts on 1   2    
      !                            1   2
      call ao_two_e_eff_dr12_pot(i,k,j,l,mu_in,d_dr12,d_dr12_large)
      accu  = 0.d0
      do m = 1, 3
       accu += d_dr12(m) - d_dr12_large(m) 
      enddo
      !                        < k l | d/dr12 | i j >
      !                           1 2 1 2 
      ao_two_e_eff_dr12_pot_array(k,l,i,j) = accu
    enddo
   enddo
  enddo
 enddo
 !$OMP END DO
 !$OMP END PARALLEL
 call wall_time(wall1)
 print*,'ao_two_e_eff_dr12_pot_array provided in ',wall1 - wall0
END_PROVIDER 

BEGIN_PROVIDER [double precision, ao_two_e_eff_dr12_pot_array_no_cycle, (ao_num,ao_num,ao_num,ao_num)]
 implicit none
 BEGIN_DOC
!   1 2                                1 2 
! < k l | [erf( mu r12) - 1] d/d_r12 | i j > on the AO basis
!
! WARNING !!!!!! ONLY WORKS FOR ATOMIC SYSTEMS, NOT MOLECULAR SYSTEMS 
 END_DOC
 integer :: i,j,k,l,m
 double precision :: mu_in,d_dr12(3),d_dr12_large(3),accu
 mu_in = mu_erf
 !$OMP PARALLEL & 
 !$OMP DEFAULT(NONE) &
 !$OMP PRIVATE(i,j,k,l,m,d_dr12,d_dr12_large,accu) & 
 !$OMP SHARED(ao_two_e_eff_dr12_pot_array_no_cycle,mu_in,ao_num) 
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
      ao_two_e_eff_dr12_pot_array_no_cycle(k,l,i,j) = accu
    enddo
   enddo
  enddo
 enddo
 !$OMP END DO
 !$OMP END PARALLEL
END_PROVIDER 


BEGIN_PROVIDER [double precision, mo_two_e_eff_dr12_pot_array_old, (mo_num,mo_num,mo_num,mo_num)]
 implicit none
 BEGIN_DOC
! mo_two_e_eff_dr12_pot_array_old(k,l,i,j) = < i j | [erf( mu r12) - 1] d/d_r12 | k l > on the MO basis
!
!
! WARNING !!!!!! ONLY WORKS FOR ATOMIC SYSTEMS, NOT MOLECULAR SYSTEMS 
 END_DOC
 integer :: i,j,k,l,m,n,p,q
 double precision, allocatable :: mo_tmp_1(:,:,:,:),mo_tmp_2(:,:,:,:),mo_tmp_3(:,:,:,:)
 double precision :: wall0,wall1
 provide ao_two_e_eff_dr12_pot_array
 call wall_time(wall0)
 print*,'Providing mo_two_e_eff_dr12_pot_array ...'
 mo_two_e_eff_dr12_pot_array_old = 0.d0
 allocate(mo_tmp_1(ao_num,ao_num,ao_num,mo_num))
 mo_tmp_1 = 0.d0
 do i = 1, mo_num ! r1
  do m = 1, ao_num
   do n = 1, ao_num
    do p = 1, ao_num
     do q = 1, ao_num
      ! <j p
      mo_tmp_1(p,n,m,i) += mo_coef(q,i) * ao_two_e_eff_dr12_pot_array(q,p,n,m)
      !  1 2          1 2 
      !mo_tmp_1(p,n,m,i) = \sum_q c_qi < q p |d/dr12| n m > = < i p | d/dr12 | n m >
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
      !mo_tmp_2(n,m,j,i) = \sum_p c_pj < i p | d/dr12 | n m > = < i j | d/dr12 | n m >
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
     !mo_tmp_2(m,k,j,i) = \sum_n c_nk < i j | d/dr12 | n m > = < i j | d/dr12 | k m >
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
     !mo_tmp_2(m,k,j,i) = \sum_m c_ml < i j | d/dr12 | k m > = < i j | d/dr12 | k l >
      mo_two_e_eff_dr12_pot_array_old(l,k,j,i) += mo_coef(m,l) * mo_tmp_3(m,k,j,i)
     enddo
    enddo
   enddo
  enddo
 enddo
 deallocate(mo_tmp_3)
 call wall_time(wall1)
 print*,'mo_two_e_eff_dr12_pot_array provided in ',wall1 - wall0

END_PROVIDER 

BEGIN_PROVIDER [double precision, mo_two_e_eff_dr12_pot_array, (mo_num,mo_num,mo_num,mo_num)]
 implicit none
 BEGIN_DOC
! mo_two_e_eff_dr12_pot_array(k,l,i,j) = < i j | [erf( mu r12) - 1] d/d_r12 | k l > on the MO basis
!
!
! WARNING !!!!!! ONLY WORKS FOR ATOMIC SYSTEMS, NOT MOLECULAR SYSTEMS 
 END_DOC
 integer :: i,j,k,l,m,n,p,q
 mo_two_e_eff_dr12_pot_array = 0.d0
 do i = 1, mo_num
  do j = 1, mo_num
   do l = 1, mo_num
    do k = 1, mo_num
     mo_two_e_eff_dr12_pot_array(i,j,k,l) = mo_two_e_eff_dr12_pot_array_old(k,l,i,j)
    enddo
   enddo
  enddo
 enddo
 FREE   mo_two_e_eff_dr12_pot_array_old 

END_PROVIDER 
