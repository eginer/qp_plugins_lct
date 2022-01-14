use map_module


integer function load_mo_tc_sym_two_e_pot_dag(filename)

  BEGIN_DOC
  ! Read from disk the |MO| eff_pot_dag integrals
  END_DOC

  implicit none

  character*(*), intent(in)        :: filename

  integer                          :: iknd, kknd
  integer*8                        :: i, n, j
  integer(cache_key_kind), pointer :: key(:)
  real(integral_kind),     pointer :: val(:)

  load_mo_tc_sym_two_e_pot_dag = 1

  open(unit = 66, file = filename, FORM = 'unformatted', STATUS = 'UNKNOWN')
    read(66, err=98, end=98) iknd, kknd
    if(iknd /= integral_kind) then
      print *,  'Wrong integrals kind in file :', iknd
      stop 1
    endif
    if(kknd /= key_kind) then
      print *,  'Wrong key kind in file :', kknd
      stop 1
    endif
    read(66, err=98, end=98) mo_tc_sym_two_e_pot_dag_map%sorted, mo_tc_sym_two_e_pot_dag_map%map_size &
                           , mo_tc_sym_two_e_pot_dag_map%n_elements
    do i = 0_8, mo_tc_sym_two_e_pot_dag_map%map_size
      read(66, err=99, end=99) mo_tc_sym_two_e_pot_dag_map%map(i)%sorted, mo_tc_sym_two_e_pot_dag_map%map(i)%map_size &
                             , mo_tc_sym_two_e_pot_dag_map%map(i)%n_elements
      call cache_map_reallocate(mo_tc_sym_two_e_pot_dag_map%map(i),mo_tc_sym_two_e_pot_dag_map%map(i)%map_size)
    enddo
    do i = 0_8, mo_tc_sym_two_e_pot_dag_map%map_size
      key => mo_tc_sym_two_e_pot_dag_map%map(i)%key
      val => mo_tc_sym_two_e_pot_dag_map%map(i)%value
      n = mo_tc_sym_two_e_pot_dag_map%map(i)%n_elements
      read(66,err=99,end=99) (key(j), j=1,n), (val(j), j=1,n)
    enddo
    call map_sort(mo_tc_sym_two_e_pot_dag_map)
    load_mo_tc_sym_two_e_pot_dag = 0
    return
  99 continue
    call map_deinit(mo_tc_sym_two_e_pot_dag_map)
  98 continue
    stop 'Problem reading mo_tc_sym_two_e_pot_dag_map file in work/'

end function load_mo_tc_sym_two_e_pot_dag




BEGIN_PROVIDER [ type(map_type), mo_tc_sym_two_e_pot_dag_map ]

  BEGIN_DOC
  ! |MO| integrals
  END_DOC

  implicit none
  integer(key_kind)      :: key_max
  integer(map_size_kind) :: sze

  call two_e_integrals_index(mo_num, mo_num, mo_num, mo_num, key_max)
  sze = key_max
  call map_init(mo_tc_sym_two_e_pot_dag_map, sze)
  print*, 'MO mo_tc_sym_two_e_pot_dag_map  initialized'
END_PROVIDER



subroutine insert_into_mo_tc_sym_two_e_pot_dag_map(n_integrals, buffer_i, buffer_values, thr)

  BEGIN_DOC
  ! Create new entry into |MO| map, or accumulate in an existing entry
  END_DOC

  use map_module

  implicit none
  integer,             intent(in)    :: n_integrals
  real(integral_kind), intent(in)    :: thr
  integer(key_kind),   intent(inout) :: buffer_i(n_integrals)
  real(integral_kind), intent(inout) :: buffer_values(n_integrals)

  call map_update(mo_tc_sym_two_e_pot_dag_map, buffer_i, buffer_values, n_integrals, thr)

end subroutine insert_into_mo_tc_sym_two_e_pot_dag_map



BEGIN_PROVIDER [ double precision, mo_tc_sym_two_e_pot_dag_cache, (0:64*64*64*64) ]

  BEGIN_DOC
  ! Cache of |MO| integrals for fast access
  END_DOC

  implicit none

  integer                        :: i, j, k, l, ii
  integer(key_kind)              :: idx
  real(integral_kind)            :: integral

  PROVIDE mo_tc_sym_two_e_pot_dag_in_map
  FREE ao_tc_sym_two_e_pot_dag_cache

 !$OMP PARALLEL DO PRIVATE (i,j,k,l,idx,ii,integral)
  do l = mo_tc_sym_two_e_pot_cache_min, mo_tc_sym_two_e_pot_cache_max
    do k = mo_tc_sym_two_e_pot_cache_min,mo_tc_sym_two_e_pot_cache_max
      do j = mo_tc_sym_two_e_pot_cache_min, mo_tc_sym_two_e_pot_cache_max
        do i = mo_tc_sym_two_e_pot_cache_min, mo_tc_sym_two_e_pot_cache_max
          !DIR$ FORCEINLINE
          call two_e_integrals_index(i,j,k,l,idx)
          !DIR$ FORCEINLINE
          call map_get(mo_tc_sym_two_e_pot_dag_map, idx, integral)
          ii = l-mo_tc_sym_two_e_pot_cache_min
          ii = ior( ishft(ii,6), k-mo_tc_sym_two_e_pot_cache_min)
          ii = ior( ishft(ii,6), j-mo_tc_sym_two_e_pot_cache_min)
          ii = ior( ishft(ii,6), i-mo_tc_sym_two_e_pot_cache_min)
          mo_tc_sym_two_e_pot_dag_cache(ii) = integral
        enddo
      enddo
    enddo
  enddo
 !$OMP END PARALLEL DO

END_PROVIDER


double precision function get_mo_tc_sym_two_e_pot_dag(i, j, k, l, map)

  BEGIN_DOC
  ! Returns one integral $\langle ij|kl \rangle$ in the |MO| basis
  END_DOC

  use map_module

  implicit none

  integer, intent(in)            :: i, j, k, l
  type(map_type), intent(inout)  :: map

  integer                        :: ii
  integer(key_kind)              :: idx
  real(integral_kind)            :: tmp

  PROVIDE mo_tc_sym_two_e_pot_dag_in_map mo_tc_sym_two_e_pot_dag_cache

  ii = l-mo_tc_sym_two_e_pot_cache_min
  ii = ior(ii, k-mo_tc_sym_two_e_pot_cache_min)
  ii = ior(ii, j-mo_tc_sym_two_e_pot_cache_min)
  ii = ior(ii, i-mo_tc_sym_two_e_pot_cache_min)
  if (iand(ii, -64) /= 0) then
    !DIR$ FORCEINLINE
    call two_e_integrals_index(i, j, k, l, idx)
    !DIR$ FORCEINLINE
    call map_get(map, idx, tmp)
    get_mo_tc_sym_two_e_pot_dag = dble(tmp)
  else
    ii = l-mo_tc_sym_two_e_pot_cache_min
    ii = ior( ishft(ii,6), k-mo_tc_sym_two_e_pot_cache_min)
    ii = ior( ishft(ii,6), j-mo_tc_sym_two_e_pot_cache_min)
    ii = ior( ishft(ii,6), i-mo_tc_sym_two_e_pot_cache_min)
    get_mo_tc_sym_two_e_pot_dag = mo_tc_sym_two_e_pot_dag_cache(ii)
  endif

end function get_mo_tc_sym_two_e_pot_dag



double precision function mo_tc_sym_two_e_pot_dag(i,j,k,l)

  BEGIN_DOC
  ! Returns one integral $\langle ij|kl \rangle$ in the |MO| basis
  END_DOC

  implicit none
  integer, intent(in) :: i, j, k, l
  double precision    :: get_mo_tc_sym_two_e_pot_dag

  PROVIDE mo_tc_sym_two_e_pot_dag_in_map mo_tc_sym_two_e_pot_dag_cache
  !DIR$ FORCEINLINE
  mo_tc_sym_two_e_pot_dag = get_mo_tc_sym_two_e_pot_dag(i, j, k, l, mo_tc_sym_two_e_pot_dag_map)

  return
end function mo_tc_sym_two_e_pot_dag



subroutine get_many_mo_tc_sym_two_e_pot_dag(j, k, l, sze, out_val,map)

  BEGIN_DOC
  ! Returns multiple integrals $\langle ij|kl \rangle$ in the |MO| basis, all
  ! i for j,k,l fixed.
  END_DOC

  use map_module

  implicit none

  integer, intent(in)            :: j, k, l, sze
  double precision, intent(out)  :: out_val(sze)
  type(map_type), intent(inout)  :: map

  integer                        :: i
  integer(key_kind)              :: hash(sze)
  real(integral_kind)            :: tmp_val(sze)

  PROVIDE mo_tc_sym_two_e_pot_dag_in_map

  do i = 1, sze
    !DIR$ FORCEINLINE
    call two_e_integrals_index(i, j, k, l, hash(i))
  enddo

  if (key_kind == 8) then
    call map_get_many(map, hash, out_val, sze)
  else
    call map_get_many(map, hash, tmp_val, sze)
    do i = 1, sze
      out_val(i) = dble(tmp_val(i))
    enddo
  endif

end subroutine get_many_mo_tc_sym_two_e_pot_dag



subroutine get_many_mo_tc_sym_two_e_pot_dag_ij(k, l, sze, out_array, map)

  BEGIN_DOC
  ! Returns multiple integrals $\langle ij|kl \rangle$ in the |MO| basis, all
  ! $\int i(1)j(2) \frac{1}{r_{12}} k(1)l(2)$
  ! i, j for k,l fixed.
  END_DOC

  use map_module

  implicit none

  integer, intent(in)              :: k,l, sze
  double precision, intent(out)    :: out_array(sze,sze)
  type(map_type), intent(inout)    :: map

  integer                          :: i,j,kk,ll,m
  integer,             allocatable :: pairs(:,:), iorder(:)
  integer(key_kind),   allocatable :: hash(:)
  real(integral_kind), allocatable :: tmp_val(:)

  PROVIDE mo_tc_sym_two_e_pot_dag_in_map

  allocate( hash(sze*sze), pairs(2,sze*sze), iorder(sze*sze), tmp_val(sze*sze) )

  kk = 0
  out_array = 0.d0
  do j = 1, sze
    do i = 1, sze
      kk += 1
      !DIR$ FORCEINLINE
      call two_e_integrals_index(i, j, k, l, hash(kk))
      pairs(1,kk) = i
      pairs(2,kk) = j
      iorder(kk) = kk
    enddo
  enddo

  if (key_kind == 8) then
    call i8radix_sort(hash,iorder,kk,-1)
  else if (key_kind == 4) then
    call iradix_sort(hash,iorder,kk,-1)
  else if (key_kind == 2) then
    call i2radix_sort(hash,iorder,kk,-1)
  endif

  call map_get_many(mo_tc_sym_two_e_pot_dag_map, hash, tmp_val, kk)

  do ll = 1, kk
    m = iorder(ll)
    i = pairs(1,m)
    j = pairs(2,m)
    out_array(i,j) = tmp_val(ll)
  enddo

  deallocate(pairs,hash,iorder,tmp_val)

end subroutine get_many_mo_tc_sym_two_e_pot_dag_ij




integer*8 function get_mo_tc_sym_two_e_pot_dag_map_size()
  BEGIN_DOC
  ! Returns the number of elements in the |MO| map
  END_DOC
  implicit none
  get_mo_tc_sym_two_e_pot_dag_map_size = mo_tc_sym_two_e_pot_dag_map%n_elements
end function get_mo_tc_sym_two_e_pot_dag_map_size
