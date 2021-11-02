 BEGIN_PROVIDER [ double precision, two_e_tc_hermit_integral_alpha, (ao_num, ao_num) ]
&BEGIN_PROVIDER [ double precision, two_e_tc_hermit_integral_beta ,  (ao_num, ao_num) ]
 use map_module
 implicit none
 BEGIN_DOC
 ! Alpha and Beta Fock matrices in AO basis set
 END_DOC

 integer                        :: i,j,k,l,k1,r,s
 integer                        :: i0,j0,k0,l0
 integer*8                      :: p,q
 double precision               :: integral, c0, c1, c2
 double precision               :: ao_two_e_integral, local_threshold
 double precision, allocatable  :: two_e_tc_hermit_integral_alpha_tmp(:,:)
 double precision, allocatable  :: two_e_tc_hermit_integral_beta_tmp(:,:)

 two_e_tc_hermit_integral_alpha = 0.d0
 two_e_tc_hermit_integral_beta  = 0.d0
 PROVIDE ao_tc_sym_two_e_pot_in_map

   integer(omp_lock_kind) :: lck(ao_num)
   integer(map_size_kind)     :: i8
   integer                        :: ii(8), jj(8), kk(8), ll(8), k2
   integer(cache_map_size_kind)   :: n_elements_max, n_elements
   integer(key_kind), allocatable :: keys(:)
   double precision, allocatable  :: values(:)

   !$OMP PARALLEL DEFAULT(NONE)                                      &
       !$OMP PRIVATE(i,j,l,k1,k,integral,ii,jj,kk,ll,i8,keys,values,n_elements_max, &
       !$OMP  n_elements,two_e_tc_hermit_integral_alpha_tmp,two_e_tc_hermit_integral_beta_tmp)&
       !$OMP SHARED(ao_num,SCF_density_matrix_ao_alpha,SCF_density_matrix_ao_beta,&
       !$OMP  ao_tc_sym_two_e_pot_map, two_e_tc_hermit_integral_alpha, two_e_tc_hermit_integral_beta)

   call get_cache_map_n_elements_max(ao_tc_sym_two_e_pot_map,n_elements_max)
   allocate(keys(n_elements_max), values(n_elements_max))
   allocate(two_e_tc_hermit_integral_alpha_tmp(ao_num,ao_num), &
            two_e_tc_hermit_integral_beta_tmp(ao_num,ao_num))
   two_e_tc_hermit_integral_alpha_tmp = 0.d0
   two_e_tc_hermit_integral_beta_tmp  = 0.d0

   !$OMP DO SCHEDULE(static,1)
   do i8=0_8,ao_tc_sym_two_e_pot_map%map_size
     n_elements = n_elements_max
     call get_cache_map(ao_tc_sym_two_e_pot_map,i8,keys,values,n_elements)
     do k1=1,n_elements
       call two_e_integrals_index_reverse(kk,ii,ll,jj,keys(k1))

       do k2=1,8
         if (kk(k2)==0) then
           cycle
         endif
         i = ii(k2)
         j = jj(k2)
         k = kk(k2)
         l = ll(k2)
         integral = (SCF_density_matrix_ao_alpha(k,l)+SCF_density_matrix_ao_beta(k,l)) * values(k1)
         two_e_tc_hermit_integral_alpha_tmp(i,j) += integral
         two_e_tc_hermit_integral_beta_tmp (i,j) += integral
         integral = values(k1)
         two_e_tc_hermit_integral_alpha_tmp(l,j) -= SCF_density_matrix_ao_alpha(k,i) * integral
         two_e_tc_hermit_integral_beta_tmp (l,j) -= SCF_density_matrix_ao_beta (k,i) * integral
       enddo
     enddo
   enddo
   !$OMP END DO NOWAIT
   !$OMP CRITICAL
   two_e_tc_hermit_integral_alpha += two_e_tc_hermit_integral_alpha_tmp
   two_e_tc_hermit_integral_beta  += two_e_tc_hermit_integral_beta_tmp
   !$OMP END CRITICAL
   deallocate(keys,values,two_e_tc_hermit_integral_alpha_tmp,two_e_tc_hermit_integral_beta_tmp)
   !$OMP END PARALLEL


END_PROVIDER


 BEGIN_PROVIDER [ double precision, two_e_tc_non_hermit_integral_alpha, (ao_num, ao_num)]
&BEGIN_PROVIDER [ double precision, two_e_tc_non_hermit_integral_beta, (ao_num, ao_num)]
 implicit none
 integer :: i,j,k,l
 double precision :: density, density_a, density_b
 two_e_tc_non_hermit_integral_alpha = 0.d0
 two_e_tc_non_hermit_integral_beta  = 0.d0
 do i = 1, ao_num
  do k = 1, ao_num
   do j = 1, ao_num
    do l = 1, ao_num
     density_a = SCF_density_matrix_ao_alpha(l,j)
     density_b = SCF_density_matrix_ao_beta(l,j)
     density   = density_a + density_b
     two_e_tc_non_hermit_integral_alpha(k,i) += density   * ao_non_hermit_term_chemist(l,j,k,i)
     two_e_tc_non_hermit_integral_beta(k,i)  += density   * ao_non_hermit_term_chemist(l,j,k,i)
     two_e_tc_non_hermit_integral_alpha(k,i) -= density_a   * ao_non_hermit_term_chemist(k,j,l,i)
     two_e_tc_non_hermit_integral_beta(k,i)  -= density_b   * ao_non_hermit_term_chemist(k,j,l,i)
    enddo
   enddo
  enddo
 enddo

END_PROVIDER 

BEGIN_PROVIDER [ double precision, Fock_matrix_tc_ao_alpha, (ao_num, ao_num)]
 implicit none
 Fock_matrix_tc_ao_alpha = two_e_tc_non_hermit_integral_alpha + two_e_tc_hermit_integral_alpha + ao_one_e_integrals
END_PROVIDER 

BEGIN_PROVIDER [ double precision, Fock_matrix_tc_mo_alpha, (mo_num,mo_num) ]
   implicit none
   BEGIN_DOC
   ! Fock matrix on the MO basis
   END_DOC
   call ao_to_mo(Fock_matrix_tc_ao_alpha,size(Fock_matrix_tc_ao_alpha,1), &
                 Fock_matrix_tc_mo_alpha,size(Fock_matrix_tc_mo_alpha,1))
END_PROVIDER

BEGIN_PROVIDER [ double precision, Fock_matrix_tc_ao_beta, (ao_num, ao_num)]
 implicit none
 Fock_matrix_tc_ao_beta = two_e_tc_non_hermit_integral_beta + two_e_tc_hermit_integral_beta + ao_one_e_integrals
END_PROVIDER 

BEGIN_PROVIDER [ double precision, Fock_matrix_tc_mo_beta, (mo_num,mo_num) ]
   implicit none
   BEGIN_DOC
   ! Fock matrix on the MO basis
   END_DOC
   call ao_to_mo(Fock_matrix_tc_ao_beta,size(Fock_matrix_tc_ao_beta,1), &
                 Fock_matrix_tc_mo_beta,size(Fock_matrix_tc_mo_beta,1))
END_PROVIDER
