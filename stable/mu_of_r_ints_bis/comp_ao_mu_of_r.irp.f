
subroutine compute_all_ijkl_ao_for_jl_erf_mu_of_r_int(j,l,ao_integrals)
 implicit none
 integer, intent(in) :: j,l
 double precision, intent(out) :: ao_integrals(ao_num,ao_num)
 integer :: i,k
 do k = 1, ao_num
  do i = 1, ao_num
   ao_integrals(i,k) = ao_ints_erf_mu_of_r_chemist(i,k,j,l)
  enddo
 enddo
end

subroutine compute_ao_erf_mu_of_r_integrals_jl(j,l,n_integrals,buffer_i,buffer_value)
  implicit none
  use map_module
  BEGIN_DOC
  !  Parallel client for AO integrals
  END_DOC

  integer, intent(in)             :: j,l
  integer,intent(out)             :: n_integrals
  integer(key_kind),intent(out)   :: buffer_i(ao_num*ao_num)
  real(integral_kind),intent(out) :: buffer_value(ao_num*ao_num)
  logical, external               :: ao_two_e_integral_zero

  integer                         :: i,k
  double precision                :: ao_two_e_integral,cpu_1,cpu_2, wall_1, wall_2
  double precision                :: integral, wall_0
  double precision                :: thr
  integer                         :: kk, m, j1, i1
  double precision, allocatable :: ao_integrals(:,:)
  allocate( ao_integrals(ao_num,ao_num) )

  call compute_all_ijkl_ao_for_jl_erf_mu_of_r_int(j,l,ao_integrals)
  thr = ao_integrals_threshold

  n_integrals = 0

  j1 = j+shiftr(l*l-l,1)
  do k = 1, ao_num           ! r1
    i1 = shiftr(k*k-k,1)
    if (i1 > j1) then
      exit
    endif
    do i = 1, k
      i1 += 1
      if (i1 > j1) then
        exit
      endif
      if (ao_two_e_integral_zero(i,j,k,l)) then
        cycle
      endif
      integral = ao_integrals(i,k)
      !DIR$ FORCEINLINE
      if (abs(integral) < thr) then
        cycle
      endif
      n_integrals += 1
      !DIR$ FORCEINLINE
      call two_e_integrals_index(i,j,k,l,buffer_i(n_integrals))
      buffer_value(n_integrals) =  integral 
    enddo
  enddo

end

BEGIN_PROVIDER [ double precision, ao_ints_erf_mu_of_r_no_sym, (ao_num, ao_num, ao_num, ao_num)]
 implicit none
! BEGIN_DOC
! enters with the array ao_ints_erf_mu_of_r_no_sym and add the following integrals 
!
! ao_ints_erf_mu_of_r_no_sym(i,k,j,l) += \int dr1 \phi_i(r1) \phi_k(r1) \int dr2 erf(\mu(r1) r12)/r12 \phi_j(r2) \phi_l(r2)
! END_DOC
 include 'constants.include.F'
 integer :: i,k,ipoint,m
 double precision :: weight,mu
 double precision, allocatable :: a_mat(:,:,:)
! double precision, intent(inout) :: ao_ints_erf_mu_of_r_no_sym(ao_num,ao_num,ao_num,ao_num)
 

 print*,'computing lr_int_mu_r1 ...'
 call wall_time(wall0)
 double precision :: wall0,wall1

 allocate(a_mat(ao_num,ao_num,n_points_final_grid))

 ao_ints_erf_mu_of_r_no_sym = 0.d0
 do ipoint = 1, n_points_final_grid
  weight = final_weight_at_r_vector(ipoint)
  do k = 1, ao_num
   do i = 1, ao_num
    a_mat(i,k, ipoint) = aos_in_r_array(k,ipoint) * aos_in_r_array(i,ipoint) * weight
   enddo
  enddo
 enddo

 ! erf(mu(r) * r12)/r12
 call dgemm("N","N",ao_num*ao_num,ao_num*ao_num,n_points_final_grid,1.d0,a_mat(1,1,1),ao_num*ao_num &                                  
                   ,erf_mu_of_r_ao_ao_transp(1,1,1),n_points_final_grid,1.d0,ao_ints_erf_mu_of_r_no_sym,ao_num*ao_num)

 call wall_time(wall1)
 print*,'wall time for lr_int_mu_r1 ',wall1 - wall0

END_PROVIDER 

BEGIN_PROVIDER [ double precision, ao_ints_erf_mu_of_r_chemist, (ao_num, ao_num, ao_num, ao_num)]
 implicit none
 integer :: i,j,k,l
 do l = 1, ao_num
  do j = 1, ao_num
   do k = 1, ao_num
    do i = 1, ao_num
     ao_ints_erf_mu_of_r_chemist(i,k,j,l) = 0.5d0 * (ao_ints_erf_mu_of_r_no_sym(i,k,j,l) + ao_ints_erf_mu_of_r_no_sym(j,l,i,k))
    enddo
   enddo
  enddo
 enddo

END_PROVIDER 

BEGIN_PROVIDER [ double precision, ao_ints_erf_mu_of_r_physicist, (ao_num, ao_num, ao_num, ao_num)]
 implicit none
 integer :: i,j,k,l
 do l = 1, ao_num
  do k = 1, ao_num
   do j = 1, ao_num
    do i = 1, ao_num
     ao_ints_erf_mu_of_r_physicist(i,j,k,l) = 0.5d0 * (ao_ints_erf_mu_of_r_no_sym(i,k,j,l) + ao_ints_erf_mu_of_r_no_sym(j,l,i,k))
    enddo
   enddo
  enddo
 enddo

END_PROVIDER 

BEGIN_PROVIDER [double precision, erf_mu_of_r_ao_ao,( ao_num, ao_num,n_points_final_grid)]
 implicit none
 BEGIN_DOC
! erf_mu_of_r_ao_ao(i,j,r) = \int dr' phi_i(r') phi_j(r') erf(mu(r) |r - r'|)/|r-r'|
 END_DOC
 double precision :: mu,r(3),int_mu,delta,wall0,wall1,phi_j_erf_mu_r_phi
 integer :: i,j,ipoint
  provide mu_of_r_dft
 call wall_time(wall0)
  !$OMP PARALLEL                  &
  !$OMP DEFAULT (NONE)            &
  !$OMP PRIVATE (i,j,ipoint,mu,r,int_mu,delta) & 
  !$OMP SHARED (ao_num,n_points_final_grid,mu_of_r_dft,erf_mu_of_r_ao_ao,final_grid_points)
  !$OMP DO SCHEDULE (dynamic)
  do ipoint = 1, n_points_final_grid
   do i = 1, ao_num
    do j = i, ao_num
     mu = mu_of_r_dft(ipoint)
     r(1) = final_grid_points(1,ipoint)
     r(2) = final_grid_points(2,ipoint)
     r(3) = final_grid_points(3,ipoint)
     int_mu = phi_j_erf_mu_r_phi(i,j,mu, r)
     erf_mu_of_r_ao_ao(j,i,ipoint) = int_mu
    enddo
   enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

 do ipoint = 1, n_points_final_grid
  do i = 1, ao_num
   do j = 1, i-1
    erf_mu_of_r_ao_ao(j,i,ipoint)= erf_mu_of_r_ao_ao(i,j,ipoint)
   enddo
  enddo
 enddo                                                                                                                                                                 

 call wall_time(wall1)
 print*,'wall time for erf_mu_of_r_ao_ao  ',wall1 - wall0

END_PROVIDER

BEGIN_PROVIDER [double precision, erf_mu_of_r_ao_ao_transp,( n_points_final_grid,ao_num, ao_num)]
 implicit none
 integer :: ipoint, i, j
 do ipoint = 1, n_points_final_grid
  do i = 1, ao_num
   do j = 1, ao_num
    erf_mu_of_r_ao_ao_transp(ipoint,j,i) = erf_mu_of_r_ao_ao(j,i,ipoint)
   enddo
  enddo
 enddo

END_PROVIDER 
