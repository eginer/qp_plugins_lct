
BEGIN_PROVIDER [ logical, mo_tc_sym_two_e_pot_dag_in_map ]

  BEGIN_DOC
  ! If True, the map of MO two-electron integrals is provided
  END_DOC

  use map_module

  implicit none
  integer*8 :: get_mo_tc_sym_two_e_pot_dag_map_size, mo_eff_pot_dag_map_size

  PROVIDE ao_tc_sym_two_e_pot_dag_in_map

  mo_tc_sym_two_e_pot_dag_in_map = .True.

  call add_integrals_to_mo_tc_sym_two_e_pot_dag_map(full_ijkl_bitmask_4)
  mo_eff_pot_dag_map_size = get_mo_tc_sym_two_e_pot_dag_map_size()

END_PROVIDER



subroutine clear_mo_eff_pot_dag_map()
  BEGIN_DOC
  ! Frees the memory of the MO map
  END_DOC
  implicit none
  call map_deinit(mo_tc_sym_two_e_pot_dag_map)
  FREE mo_tc_sym_two_e_pot_dag_in_map
end subroutine clear_mo_eff_pot_dag_map



subroutine add_integrals_to_mo_tc_sym_two_e_pot_dag_map(mask_ijkl)

  BEGIN_DOC
  ! Adds integrals to tha MO map according to some bitmask
  END_DOC

  use bitmasks

  implicit none

  integer(bit_kind), intent(in)    :: mask_ijkl(N_int,4)

  double precision, parameter      :: thr_coef = 1.d-10

  character*(2048)                 :: output(1)
  integer*8                        :: get_mo_tc_sym_two_e_pot_dag_map_size, mo_map_size
  integer                          :: index_needed
  integer                          :: i, j, k, l, i0, j0, k0, l0, i2, i3, i4
  integer                          :: n_i, n_j, n_k, n_l
  integer                          :: i1, j1, k1, l1, ii1, kmax, thread_num
  integer                          :: n_integrals, size_buffer
  double precision                 :: c, cpu_1, cpu_2, wall_1, wall_2, wall_0
  double precision                 :: map_mb, accu_bis
  integer,             allocatable :: list_ijkl(:,:)
  integer,             allocatable :: two_e_tmp_0_idx(:)
  double precision,    allocatable :: two_e_tmp_1(:)
  double precision,    allocatable :: two_e_tmp_2(:,:)
  double precision,    allocatable :: two_e_tmp_3(:,:,:)
  integer(key_kind),   allocatable :: buffer_i(:)
  real(integral_kind), allocatable :: buffer_value(:)
  real(integral_kind), allocatable :: two_e_tmp_0(:,:)

  !DIR$ ATTRIBUTES ALIGN : 64      :: two_e_tmp_1, two_e_tmp_2, two_e_tmp_3

  PROVIDE ao_two_e_integrals_in_map  mo_coef

  !Get list of MOs for i,j,k and l
  !-------------------------------

  allocate( list_ijkl(mo_num,4) )
  call bitstring_to_list( mask_ijkl(1,1), list_ijkl(1,1), n_i, N_int )
  call bitstring_to_list( mask_ijkl(1,2), list_ijkl(1,2), n_j, N_int )
  call bitstring_to_list( mask_ijkl(1,3), list_ijkl(1,3), n_k, N_int )
  call bitstring_to_list( mask_ijkl(1,4), list_ijkl(1,4), n_l, N_int )

  print*, 'i'
  call bitstring_to_str( output(1), mask_ijkl(1,1), N_int )
  print *,  trim(output(1))
  j = 0
  do i = 1, N_int
    j += popcnt(mask_ijkl(i,1))
  enddo
  if(j==0) return
  print*, 'j'
  call bitstring_to_str( output(1), mask_ijkl(1,2), N_int )
  print *,  trim(output(1))
  j = 0
  do i = 1, N_int
    j += popcnt(mask_ijkl(i,2))
  enddo
  if(j==0) return
  print*, 'k'
  call bitstring_to_str( output(1), mask_ijkl(1,3), N_int )
  print *,  trim(output(1))
  j = 0
  do i = 1, N_int
    j += popcnt(mask_ijkl(i,3))
  enddo
  if(j==0) return
  print*, 'l'
  call bitstring_to_str( output(1), mask_ijkl(1,4), N_int )
  print *,  trim(output(1))
  j = 0
  do i = 1, N_int
    j += popcnt(mask_ijkl(i,4))
  enddo
  if(j==0) return

  size_buffer = min(ao_num*ao_num*ao_num,16000000)
  print*, 'Providing the eff_pot_dag molecular integrals '
  print*, 'Buffers : ' &
        , 8.*(mo_num*(n_j)*(n_k+1) + mo_num+ao_num+ao_num*ao_num+ size_buffer*3)/(1024*1024), 'MB / core'

  call wall_time(wall_1)
  call cpu_time(cpu_1)

  accu_bis = 0.d0

  !$OMP PARALLEL PRIVATE(l1,k1,j1,i1,i2,i3,i4,i,j,k,l,c, ii1,kmax,   &
      !$OMP  two_e_tmp_0_idx, two_e_tmp_0, two_e_tmp_1,two_e_tmp_2,two_e_tmp_3,&
      !$OMP  buffer_i,buffer_value,n_integrals,wall_2,i0,j0,k0,l0,   &
      !$OMP  wall_0,thread_num,accu_bis)                             &
      !$OMP  DEFAULT(NONE)                                           &
      !$OMP  SHARED(size_buffer,ao_num,mo_num,n_i,n_j,n_k,n_l,   &
      !$OMP  mo_coef_transp,                                         &
      !$OMP  mo_coef_transp_is_built, list_ijkl,                     &
      !$OMP  mo_coef_is_built, wall_1,                               &
      !$OMP  mo_coef,mo_integrals_threshold,mo_tc_sym_two_e_pot_dag_map)
  n_integrals = 0
  wall_0 = wall_1
  allocate(two_e_tmp_3(mo_num, n_j, n_k), &
      two_e_tmp_1(mo_num),                &
      two_e_tmp_0(ao_num,ao_num),         &
      two_e_tmp_0_idx(ao_num),            &
      two_e_tmp_2(mo_num, n_j),           &
      buffer_i(size_buffer),              &
      buffer_value(size_buffer) )

  thread_num = 0
  !$  thread_num = omp_get_thread_num()
  !$OMP DO SCHEDULE(guided)
  do l1 = 1, ao_num
    two_e_tmp_3 = 0.d0
    do k1 = 1,ao_num
      two_e_tmp_2 = 0.d0
      do j1 = 1, ao_num
        call get_many_ao_tc_sym_two_e_pot_dag(j1, k1, l1, ao_num, two_e_tmp_0(1,j1)) ! all integrals for a given l1, k1
      enddo
      do j1 = 1, ao_num
        kmax = 0
        do i1 = 1, ao_num
          c = two_e_tmp_0(i1,j1)
          if(c == 0.d0) cycle
          kmax += 1
          two_e_tmp_0(kmax,j1)  = c
          two_e_tmp_0_idx(kmax) = i1
        enddo

        if(kmax==0) cycle

        two_e_tmp_1 = 0.d0
        ii1 = 1
        ! sum_m c_m^i (m)
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
      j1 = ishft((l*l-l),-1)
      do j0 = 1, n_j
        j = list_ijkl(j0,2)
        if (j > l)  then
          exit
        endif
        j1 += 1
        do k0 = 1, n_k
          k = list_ijkl(k0,3)
          i1 = ishft((k*k-k),-1)
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
              call insert_into_mo_tc_sym_two_e_pot_dag_map(n_integrals,buffer_i,buffer_value,&
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
            wall_2-wall_1, 's', map_mb(mo_tc_sym_two_e_pot_dag_map) ,'MB'
      endif
    endif
  enddo
  !$OMP END DO NOWAIT
  deallocate (two_e_tmp_1,two_e_tmp_2,two_e_tmp_3)


  call insert_into_mo_tc_sym_two_e_pot_dag_map(n_integrals,buffer_i,buffer_value,&
      real(mo_integrals_threshold,integral_kind))
  deallocate(buffer_i, buffer_value)
  !$OMP END PARALLEL
  call map_merge(mo_tc_sym_two_e_pot_dag_map)

  call wall_time(wall_2)
  call cpu_time(cpu_2)
  mo_map_size = get_mo_tc_sym_two_e_pot_dag_map_size()

  deallocate(list_ijkl)


  print*,'Molecular eff_pot_dag integrals provided:'
  print*,' Size of MO eff_pot_dag map           ', map_mb(mo_tc_sym_two_e_pot_dag_map) ,'MB'
  print*,' Number of MO eff_pot_dag integrals: ',  mo_map_size
  print*,' cpu  time :',cpu_2 - cpu_1, 's'
  print*,' wall time :',wall_2 - wall_1, 's  ( x ', (cpu_2-cpu_1)/(wall_2-wall_1), ')'

end


