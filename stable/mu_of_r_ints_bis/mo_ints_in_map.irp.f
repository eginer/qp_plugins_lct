BEGIN_PROVIDER [ logical, mo_int_erf_mu_of_r_in_map ]
  use map_module
  implicit none
  BEGIN_DOC
  ! If True, the map of MO two-electron integrals is provided
  END_DOC
  integer(bit_kind)              :: mask_ijkl(N_int,4)
  integer(bit_kind)              :: mask_ijk(N_int,3)
  double precision               :: cpu_1, cpu_2, wall_1, wall_2

  PROVIDE mo_class

  mo_int_erf_mu_of_r_in_map = .True.
!  if (read_mo_int_erf_mu_of_r) then
!    print*,'Reading the MO integrals'
!    call map_load_from_disk(trim(ezfio_filename)//'/work/mo_ints',mo_integrals_map)
!    print*, 'MO integrals provided'
!    return
!  else
    PROVIDE ao_int_erf_mu_of_r_in_map
!  endif

  print *,  ''
  print *,  'AO -> MO integrals transformation'
  print *,  '---------------------------------'
  print *,  ''

  call wall_time(wall_1)
  call cpu_time(cpu_1)

  call add_integrals_to_mo_int_erf_mu_of_r_map(full_ijkl_bitmask_4)

  call wall_time(wall_2)
  call cpu_time(cpu_2)

  integer*8                      :: get_mo_erf_mu_of_r_map_size, mo_map_size
  mo_map_size = get_mo_erf_mu_of_r_map_size ()

  double precision, external     :: map_mb
  print*,'Molecular mo_int_erf_mu_of_r_map provided :'
  print*,' Size of MO map           ', map_mb(mo_int_erf_mu_of_r_map) ,'MB'
  print*,' Number of MO mo_int_erf_mu_of_r_map integrals : ',  mo_map_size
  print*,' cpu  time :',cpu_2 - cpu_1, 's'
  print*,' wall time :',wall_2 - wall_1, 's  ( x ', (cpu_2-cpu_1)/(wall_2-wall_1), ')'

!  if (write_mo_int_erf_mu_of_r.and.mpi_master) then
!    call ezfio_set_work_empty(.False.)
!    call map_save_to_disk(trim(ezfio_filename)//'/work/mo_ints',mo_integrals_map)
!    call ezfio_set_mo_two_e_ints_io_mo_int_erf_mu_of_r('Read')
!  endif

END_PROVIDER


subroutine add_integrals_to_mo_int_erf_mu_of_r_map(mask_ijkl)
  use bitmasks
  implicit none

  BEGIN_DOC
  ! Adds integrals to tha MO map according to some bitmask
  END_DOC

  integer(bit_kind), intent(in)  :: mask_ijkl(N_int,4)

  integer                        :: i,j,k,l
  integer                        :: i0,j0,k0,l0
  double precision               :: c, cpu_1, cpu_2, wall_1, wall_2, wall_0

  integer, allocatable           :: list_ijkl(:,:)
  integer                        :: n_i, n_j, n_k, n_l
  integer, allocatable           :: two_e_tmp_0_idx(:)
  real(integral_kind), allocatable :: two_e_tmp_0(:,:)
  double precision, allocatable  :: two_e_tmp_1(:)
  double precision, allocatable  :: two_e_tmp_2(:,:)
  double precision, allocatable  :: two_e_tmp_3(:,:,:)
  !DIR$ ATTRIBUTES ALIGN : 64    :: two_e_tmp_1, two_e_tmp_2, two_e_tmp_3

  integer                        :: n_integrals
  integer                        :: size_buffer
  integer(key_kind),allocatable  :: buffer_i(:)
  real(integral_kind),allocatable :: buffer_value(:)
  double precision, external     :: map_mb

  integer                        :: i1,j1,k1,l1, ii1, kmax, thread_num
  integer                        :: i2,i3,i4
  double precision,parameter     :: thr_coef = 1.d-10

  PROVIDE ao_int_erf_mu_of_r_in_map  mo_coef

  !Get list of MOs for i,j,k and l
  !-------------------------------

  allocate(list_ijkl(mo_num,4))
  call bitstring_to_list( mask_ijkl(1,1), list_ijkl(1,1), n_i, N_int )
  call bitstring_to_list( mask_ijkl(1,2), list_ijkl(1,2), n_j, N_int )
  call bitstring_to_list( mask_ijkl(1,3), list_ijkl(1,3), n_k, N_int )
  call bitstring_to_list( mask_ijkl(1,4), list_ijkl(1,4), n_l, N_int )
  j = 0
  do i = 1, N_int
    j += popcnt(mask_ijkl(i,1))
  enddo
  if(j==0)then
    return
  endif

  j = 0
  do i = 1, N_int
    j += popcnt(mask_ijkl(i,2))
  enddo
  if(j==0)then
    return
  endif

  j = 0
  do i = 1, N_int
    j += popcnt(mask_ijkl(i,3))
  enddo
  if(j==0)then
    return
  endif

  j = 0
  do i = 1, N_int
    j += popcnt(mask_ijkl(i,4))
  enddo
  if(j==0)then
    return
  endif

  size_buffer = min(ao_num*ao_num*ao_num,16000000)
  print*, 'Buffers : ', 8.*(mo_num*(n_j)*(n_k+1) + mo_num+&
      ao_num+ao_num*ao_num+ size_buffer*3)/(1024*1024), 'MB / core'

  double precision               :: accu_bis
  accu_bis = 0.d0
  call wall_time(wall_1)

  !$OMP PARALLEL PRIVATE(l1,k1,j1,i1,i2,i3,i4,i,j,k,l,c, ii1,kmax,   &
      !$OMP  two_e_tmp_0_idx, two_e_tmp_0, two_e_tmp_1,two_e_tmp_2,two_e_tmp_3,&
      !$OMP  buffer_i,buffer_value,n_integrals,wall_2,i0,j0,k0,l0,   &
      !$OMP  wall_0,thread_num,accu_bis)                             &
      !$OMP  DEFAULT(NONE)                                           &
      !$OMP  SHARED(size_buffer,ao_num,mo_num,n_i,n_j,n_k,n_l,   &
      !$OMP  mo_coef_transp,                                         &
      !$OMP  mo_coef_transp_is_built, list_ijkl,                     &
      !$OMP  mo_coef_is_built, wall_1,                               &
      !$OMP  mo_coef,mo_integrals_threshold,mo_int_erf_mu_of_r_map)
  n_integrals = 0
  wall_0 = wall_1
  allocate(two_e_tmp_3(mo_num, n_j, n_k),                 &
      two_e_tmp_1(mo_num),                                &
      two_e_tmp_0(ao_num,ao_num),                                   &
      two_e_tmp_0_idx(ao_num),                                      &
      two_e_tmp_2(mo_num, n_j),                           &
      buffer_i(size_buffer),                                         &
      buffer_value(size_buffer) )

  thread_num = 0
  !$  thread_num = omp_get_thread_num()
  !$OMP DO SCHEDULE(guided)
  do l1 = 1,ao_num
    two_e_tmp_3 = 0.d0
    do k1 = 1,ao_num
      two_e_tmp_2 = 0.d0
      do j1 = 1,ao_num
        call get_ao_two_e_ints_erf_mu_of_r(j1,k1,l1,ao_num,two_e_tmp_0(1,j1))
      enddo
      do j1 = 1,ao_num
        kmax = 0
        do i1 = 1,ao_num
          c = two_e_tmp_0(i1,j1)
          if (c == 0.d0) then
            cycle
          endif
          kmax += 1
          two_e_tmp_0(kmax,j1) = c
          two_e_tmp_0_idx(kmax) = i1
        enddo

        if (kmax==0) then
          cycle
        endif

        two_e_tmp_1 = 0.d0
        ii1=1
        do ii1 = 1,kmax-4,4
          i1 = two_e_tmp_0_idx(ii1)
          i2 = two_e_tmp_0_idx(ii1+1)
          i3 = two_e_tmp_0_idx(ii1+2)
          i4 = two_e_tmp_0_idx(ii1+3)
          do i = list_ijkl(1,1), list_ijkl(n_i,1)
            two_e_tmp_1(i)  =  two_e_tmp_1(i) +                    &
                mo_coef_transp(i,i1) * two_e_tmp_0(ii1,j1) +        &
                mo_coef_transp(i,i2) * two_e_tmp_0(ii1+1,j1) +      &
                mo_coef_transp(i,i3) * two_e_tmp_0(ii1+2,j1) +      &
                mo_coef_transp(i,i4) * two_e_tmp_0(ii1+3,j1)
          enddo ! i
        enddo  ! ii1

        i2 = ii1
        do ii1 = i2,kmax
          i1 = two_e_tmp_0_idx(ii1)
          do i = list_ijkl(1,1), list_ijkl(n_i,1)
            two_e_tmp_1(i) = two_e_tmp_1(i) + mo_coef_transp(i,i1) * two_e_tmp_0(ii1,j1)
          enddo ! i
        enddo  ! ii1
        c = 0.d0

        do i = list_ijkl(1,1), list_ijkl(n_i,1)
          c = max(c,abs(two_e_tmp_1(i)))
          if (c>mo_integrals_threshold) exit
        enddo
        if ( c < mo_integrals_threshold ) then
          cycle
        endif

        do j0 = 1, n_j
          j = list_ijkl(j0,2)
          c = mo_coef_transp(j,j1)
          if (abs(c) < thr_coef) then
            cycle
          endif
          do i = list_ijkl(1,1), list_ijkl(n_i,1)
            two_e_tmp_2(i,j0)  = two_e_tmp_2(i,j0) + c * two_e_tmp_1(i)
          enddo ! i
        enddo  ! j
      enddo !j1
      if ( maxval(abs(two_e_tmp_2)) < mo_integrals_threshold ) then
        cycle
      endif


      do k0 = 1, n_k
        k = list_ijkl(k0,3)
        c = mo_coef_transp(k,k1)
        if (abs(c) < thr_coef) then
          cycle
        endif

        do j0 = 1, n_j
          j = list_ijkl(j0,2)
          do i = list_ijkl(1,1), k
            two_e_tmp_3(i,j0,k0) = two_e_tmp_3(i,j0,k0) + c* two_e_tmp_2(i,j0)
          enddo!i
        enddo !j

      enddo  !k
    enddo   !k1



    do l0 = 1,n_l
      l = list_ijkl(l0,4)
      c = mo_coef_transp(l,l1)
      if (abs(c) < thr_coef) then
        cycle
      endif
      j1 = shiftr((l*l-l),1)
      do j0 = 1, n_j
        j = list_ijkl(j0,2)
        if (j > l)  then
          exit
        endif
        j1 += 1
        do k0 = 1, n_k
          k = list_ijkl(k0,3)
          i1 = shiftr((k*k-k),1)
          if (i1<=j1) then
            continue
          else
            exit
          endif
          two_e_tmp_1 = 0.d0
          do i0 = 1, n_i
            i = list_ijkl(i0,1)
            if (i>k) then
              exit
            endif
            two_e_tmp_1(i) = c*two_e_tmp_3(i,j0,k0)
            !           i1+=1
          enddo

          do i0 = 1, n_i
            i = list_ijkl(i0,1)
            if(i> min(k,j1-i1+list_ijkl(1,1)-1))then
              exit
            endif
            if (abs(two_e_tmp_1(i)) < mo_integrals_threshold) then
              cycle
            endif
            n_integrals += 1
            buffer_value(n_integrals) = two_e_tmp_1(i)
            !DIR$ FORCEINLINE
            call mo_two_e_integrals_index(i,j,k,l,buffer_i(n_integrals))
            if (n_integrals == size_buffer) then
              call insert_into_mo_int_erf_mu_of_r_map(n_integrals,buffer_i,buffer_value,&
                  real(mo_integrals_threshold,integral_kind))
              n_integrals = 0
            endif
          enddo
        enddo
      enddo
    enddo

    call wall_time(wall_2)
    if (thread_num == 0) then
      if (wall_2 - wall_0 > 1.d0) then
        wall_0 = wall_2
        print*, 100.*float(l1)/float(ao_num), '% in ',               &
            wall_2-wall_1, 's', map_mb(mo_int_erf_mu_of_r_map) ,'MB'
      endif
    endif
  enddo
  !$OMP END DO NOWAIT
  deallocate (two_e_tmp_1,two_e_tmp_2,two_e_tmp_3)

  integer                        :: index_needed

  call insert_into_mo_int_erf_mu_of_r_map(n_integrals,buffer_i,buffer_value,&
      real(mo_integrals_threshold,integral_kind))
  deallocate(buffer_i, buffer_value)
  !$OMP END PARALLEL
  call map_merge(mo_int_erf_mu_of_r_map)

  call wall_time(wall_2)
  call cpu_time(cpu_2)
  integer*8                      :: get_mo_erf_mu_of_r_map_size, mo_map_size
  mo_map_size = get_mo_erf_mu_of_r_map_size()

  deallocate(list_ijkl)


end


