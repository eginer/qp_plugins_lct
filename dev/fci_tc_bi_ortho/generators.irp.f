use bitmasks

 BEGIN_PROVIDER [ double precision, psi_coef_sorted_gen, (N_det,N_states)]
&BEGIN_PROVIDER [ integer(bit_kind), psi_det_sorted_gen, (N_int,2,N_det) ]
&BEGIN_PROVIDER [ integer, psi_det_sorted_gen_order,     (N_det)  ]

 implicit none
 integer :: i
 integer, allocatable :: iorder(:)
 double precision :: accu
 double precision, allocatable :: coef(:)
 allocate(iorder(N_det), coef(N_det))
 accu = 0.d0
 do i = 1, N_det
!  accu += dabs(psi_r_coef_bi_ortho(i,1) * psi_r_coef_bi_ortho(i,1))
  accu += dabs(psi_l_coef_bi_ortho(i,1) * psi_r_coef_bi_ortho(i,1))
 enddo
 print*,'accu = ',accu
 accu = 1.d0/dsqrt(accu)
 do i = 1, N_det
!  coef(i) = -dsqrt(dabs(psi_r_coef_bi_ortho(i,1) * psi_r_coef_bi_ortho(i,1))) * accu
  coef(i) = -dsqrt(dabs(psi_l_coef_bi_ortho(i,1) * psi_r_coef_bi_ortho(i,1))) * accu
  iorder(i) = i
 enddo
 accu = 0.D0
 do i = 1, N_det
  accu += coef(i)*coef(i)
 enddo
 print*,'accu = ',accu
 call dsort(coef,iorder,N_det)
 do i = 1, N_det
  psi_det_sorted_gen_order(i) = iorder(i)
  psi_coef_sorted_gen(i,1) = -coef(i)
  psi_det_sorted_gen(1:N_int,1:2,i) = psi_det(1:N_int,1:2,iorder(i))
 enddo
END_PROVIDER 

 BEGIN_PROVIDER [ integer(bit_kind), psi_det_generators, (N_int,2,psi_det_size) ]
&BEGIN_PROVIDER [ double precision, psi_coef_generators, (psi_det_size,N_states) ]
 implicit none
 BEGIN_DOC
 ! For Single reference wave functions, the generator is the
 ! Hartree-Fock determinant
 END_DOC
 psi_det_generators(1:N_int,1:2,1:N_det) = psi_det_sorted_gen(1:N_int,1:2,1:N_det)
 psi_coef_generators(1:N_det,1:N_states) = psi_coef_sorted_gen(1:N_det,1:N_states)

END_PROVIDER                                                                          

BEGIN_PROVIDER [ integer, N_det_generators ]                                                                              
 implicit none
 BEGIN_DOC
 ! For Single reference wave functions, the number of generators is 1 : the
 ! Hartree-Fock determinant
 END_DOC
 integer :: i
 double precision :: norm
 call write_time(6)
 norm = 1.d0
 N_det_generators = N_det
! do i=1,N_det
!   norm = norm - psi_average_norm_contrib_sorted(i)
!   if (norm - 1.d-10 < 1.d0 - threshold_generators) then
!     N_det_generators = i
!     exit
!   endif
! enddo
! N_det_generators = max(N_det_generators,1)
 call write_int(6,N_det_generators,'Number of generators')
END_PROVIDER


 BEGIN_PROVIDER [double precision, reigvec_tc_bi_orth_sorted, (psi_det_size, N_states)]
&BEGIN_PROVIDER [double precision, leigvec_tc_bi_orth_sorted, (psi_det_size, N_states)]

   implicit none
   integer                       :: i, j, k
   reigvec_tc_bi_orth_sorted = 0.d0
   leigvec_tc_bi_orth_sorted = 0.d0
   do i = 1, N_det
    reigvec_tc_bi_orth_sorted(i,1) = psi_r_coef_bi_ortho(psi_det_sorted_gen_order(i),1)
    leigvec_tc_bi_orth_sorted(i,1) = psi_l_coef_bi_ortho(psi_det_sorted_gen_order(i),1)
   enddo
!   integer, allocatable          :: iorder(:)
!   double precision, allocatable :: tmp_sort(:)
!
!   allocate ( tmp_sort(N_det), iorder(N_det) )
!
!   do i = 1, N_det
!     tmp_sort(i) = -dabs(reigvec_tc_bi_orth(i,1) * leigvec_tc_bi_orth(i,1))
!     iorder(i) = i
!   enddo
!   call dsort(tmp_sort, iorder, N_det)
!   deallocate(tmp_sort)
!
!   do k = 1, N_states
!     do i = 1, N_det
!       j = iorder(i)
!       reigvec_tc_bi_orth_sorted(i,k) = reigvec_tc_bi_orth(j,k)
!       leigvec_tc_bi_orth_sorted(i,k) = leigvec_tc_bi_orth(j,k)
!     enddo
!   enddo
!
!   reigvec_tc_bi_orth_sorted(N_det+1:psi_det_size,:) = 0.d0
!   leigvec_tc_bi_orth_sorted(N_det+1:psi_det_size,:) = 0.d0
!
!   deallocate(iorder)

END_PROVIDER


