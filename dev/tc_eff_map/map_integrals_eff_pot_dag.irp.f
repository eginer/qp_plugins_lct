use map_module

!! mo Map
!! ======

BEGIN_PROVIDER [ type(map_type), mo_integrals_tcdag_int_map ]

  BEGIN_DOC
  ! |mo| integrals
  END_DOC

  implicit none
  integer(key_kind)      :: key_max
  integer(map_size_kind) :: sze

  call two_e_integrals_index(mo_num,mo_num,mo_num,mo_num,key_max)
  sze = key_max
  call map_init(mo_integrals_tcdag_int_map, sze)
  print*,  'mo tcdag_int map initialized : ', sze

END_PROVIDER



BEGIN_PROVIDER [ double precision, mo_integrals_tcdag_int_cache, (0:64*64*64*64) ]

  BEGIN_DOC
  ! Cache of |mo| integrals for fast access
  END_DOC

  use map_module

  implicit none
  integer             :: i, j, k, l, ii
  integer(key_kind)   :: idx
  real(integral_kind) :: integral

  PROVIDE mo_two_e_integrals_tcdag_int_in_map

 !$OMP PARALLEL DO PRIVATE (i,j,k,l,idx,ii,integral)
  do l = mo_integrals_tc_int_cache_min, mo_integrals_tc_int_cache_max
    do k = mo_integrals_tc_int_cache_min, mo_integrals_tc_int_cache_max
      do j = mo_integrals_tc_int_cache_min, mo_integrals_tc_int_cache_max
        do i = mo_integrals_tc_int_cache_min, mo_integrals_tc_int_cache_max
          !DIR$ FORCEINLINE
          call two_e_integrals_index(i, j, k, l, idx)
          !DIR$ FORCEINLINE
          call map_get(mo_integrals_tcdag_int_map, idx, integral)
          ii = l - mo_integrals_tc_int_cache_min
          ii = ior( ishft(ii,6), k-mo_integrals_tc_int_cache_min)
          ii = ior( ishft(ii,6), j-mo_integrals_tc_int_cache_min)
          ii = ior( ishft(ii,6), i-mo_integrals_tc_int_cache_min)
          mo_integrals_tcdag_int_cache(ii) = integral
        enddo
      enddo
    enddo
  enddo
 !$OMP END PARALLEL DO

END_PROVIDER



subroutine insert_into_mo_integrals_tcdag_int_map(n_integrals, buffer_i, buffer_values)

  BEGIN_DOC
  ! Create new entry into |mo| map
  END_DOC

  use map_module

  implicit none
  integer, intent(in)                :: n_integrals
  integer(key_kind),   intent(inout) :: buffer_i(n_integrals)
  real(integral_kind), intent(inout) :: buffer_values(n_integrals)

  call map_append(mo_integrals_tcdag_int_map, buffer_i, buffer_values, n_integrals)

end subroutine insert_into_mo_integrals_tcdag_int_map



double precision function get_mo_two_e_integral_tcdag_int(i, j, k, l, map) result(result)

  BEGIN_DOC
  ! Gets one |mo| two-electron integral from the |mo| map
  END_DOC

  use map_module

  implicit none
  integer, intent(in)            :: i, j, k, l
  type(map_type), intent(inout)  :: map
  integer                        :: ii
  integer(key_kind)              :: idx
  real(integral_kind)            :: tmp

  PROVIDE mo_two_e_integrals_tcdag_int_in_map mo_integrals_tcdag_int_cache mo_integrals_tc_int_cache_min

  !DIR$ FORCEINLINE
  ii = l-mo_integrals_tc_int_cache_min
  ii = ior(ii, k-mo_integrals_tc_int_cache_min)
  ii = ior(ii, j-mo_integrals_tc_int_cache_min)
  ii = ior(ii, i-mo_integrals_tc_int_cache_min)
  if (iand(ii, -64) /= 0) then
    !DIR$ FORCEINLINE
    call two_e_integrals_index(i, j, k, l, idx)
    !DIR$ FORCEINLINE
    call map_get(map,idx,tmp)
    tmp = tmp
  else
    ii = l-mo_integrals_tc_int_cache_min
    ii = ior( ishft(ii,6), k-mo_integrals_tc_int_cache_min)
    ii = ior( ishft(ii,6), j-mo_integrals_tc_int_cache_min)
    ii = ior( ishft(ii,6), i-mo_integrals_tc_int_cache_min)
    tmp = mo_integrals_tcdag_int_cache(ii)
  endif
  result = tmp

end function get_mo_two_e_integral_tcdag_int



subroutine get_mo_two_e_integrals_tcdag_int(j, k, l, sze, out_val)

  BEGIN_DOC
  ! Gets multiple |mo| two-electron integral from the |mo| map .
  ! All i are retrieved for j,k,l fixed.
  END_DOC

  use map_module

  implicit none
  integer, intent(in)              :: j, k, l, sze
  real(integral_kind), intent(out) :: out_val(sze)
  integer                          :: i
  integer(key_kind)                :: hash
  double precision                 :: thresh
  double precision                 :: get_mo_two_e_integral_tcdag_int

  PROVIDE mo_two_e_integrals_tcdag_int_in_map mo_integrals_tcdag_int_map

  thresh = mo_integrals_threshold

  if (.False.) then
    out_val = 0.d0
    return
  endif

  do i=1,sze
    out_val(i) = get_mo_two_e_integral_tcdag_int(i, j, k, l, mo_integrals_tcdag_int_map)
  enddo

end subroutine get_mo_two_e_integrals_tcdag_int



subroutine get_mo_two_e_integrals_tcdag_int_non_zero(j, k, l, sze, out_val, out_val_index, non_zero_int)

  BEGIN_DOC
  ! Gets multiple |mo| two-electron integrals from the |mo| map .
  ! All non-zero i are retrieved for j,k,l fixed.
  END_DOC

  use map_module

  implicit none
  integer, intent(in)              :: j, k, l, sze
  real(integral_kind), intent(out) :: out_val(sze)
  integer,             intent(out) :: out_val_index(sze), non_zero_int

  integer                        :: i
  integer(key_kind)              :: hash
  double precision               :: thresh, tmp

  PROVIDE mo_two_e_integrals_tcdag_int_in_map

  thresh = mo_integrals_threshold

  non_zero_int = 0
  if (.False.) then
    out_val = 0.d0
    return
  endif

  non_zero_int = 0
  do i=1,sze
    !DIR$ FORCEINLINE
    call two_e_integrals_index(i, j, k, l, hash)
    call map_get(mo_integrals_tcdag_int_map, hash,tmp)
    if (dabs(tmp) < thresh ) cycle
    non_zero_int = non_zero_int+1
    out_val_index(non_zero_int) = i
    out_val(non_zero_int) = tmp
  enddo

end subroutine get_mo_two_e_integrals_tcdag_int_non_zero



function get_mo_tcdag_int_map_size()
  BEGIN_DOC
  ! Returns the number of elements in the |mo| map
  END_DOC
  implicit none
  integer (map_size_kind) :: get_mo_tcdag_int_map_size
  get_mo_tcdag_int_map_size = mo_integrals_tcdag_int_map%n_elements
end function get_mo_tcdag_int_map_size



subroutine clear_mo_tcdag_int_map()
  BEGIN_DOC
  ! Frees the memory of the |mo| map
  END_DOC
  implicit none
  call map_deinit(mo_integrals_tcdag_int_map)
  FREE mo_integrals_tcdag_int_map
end subroutine clear_mo_tcdag_int_map



subroutine dump_mo_integrals_tcdag_int(filename)

  BEGIN_DOC
  ! Save to disk the |mo| tcdag_int integrals
  END_DOC

  use map_module

  implicit none
  character*(*), intent(in)        :: filename
  integer*8                        :: i, j, n
  integer(cache_key_kind), pointer :: key(:)
  real(integral_kind),     pointer :: val(:)

  call ezfio_set_work_empty(.False.)

  open(unit=66, file=filename, FORM='unformatted')
    write(66) integral_kind, key_kind
    write(66) mo_integrals_tcdag_int_map%sorted, mo_integrals_tcdag_int_map%map_size, mo_integrals_tcdag_int_map%n_elements
    do i = 0_8, mo_integrals_tcdag_int_map%map_size
      write(66) mo_integrals_tcdag_int_map%map(i)%sorted, mo_integrals_tcdag_int_map%map(i)%map_size &
              , mo_integrals_tcdag_int_map%map(i)%n_elements
    enddo
    do i = 0_8, mo_integrals_tcdag_int_map%map_size
      key => mo_integrals_tcdag_int_map%map(i)%key
      val => mo_integrals_tcdag_int_map%map(i)%value
      n = mo_integrals_tcdag_int_map%map(i)%n_elements
      write(66) (key(j), j=1,n), (val(j), j=1,n)
    enddo
  close(66)

end subroutine dump_mo_integrals_tcdag_int



integer function load_mo_integrals_tcdag_int(filename)

  BEGIN_DOC
  ! Read from disk the |mo| tcdag_int integrals
  END_DOC

  implicit none
  character*(*), intent(in)        :: filename
  integer*8                        :: i, j, n
  integer                          :: iknd, kknd
  integer(cache_key_kind), pointer :: key(:)
  real(integral_kind),     pointer :: val(:)

  load_mo_integrals_tcdag_int = 1

  open(unit=66, file=filename, FORM='unformatted', STATUS='UNKNOWN')
    read(66,err=98,end=98) iknd, kknd
    if (iknd /= integral_kind) then
      print *,  'Wrong integrals kind in file :', iknd
      stop 1
    endif
    if (kknd /= key_kind) then
      print *,  'Wrong key kind in file :', kknd
      stop 1
    endif
    read(66,err=98,end=98) mo_integrals_tcdag_int_map%sorted, mo_integrals_tcdag_int_map%map_size &
                         , mo_integrals_tcdag_int_map%n_elements
    do i=0_8, mo_integrals_tcdag_int_map%map_size
      read(66,err=99,end=99) mo_integrals_tcdag_int_map%map(i)%sorted, mo_integrals_tcdag_int_map%map(i)%map_size &
                           , mo_integrals_tcdag_int_map%map(i)%n_elements
      call cache_map_reallocate(mo_integrals_tcdag_int_map%map(i),mo_integrals_tcdag_int_map%map(i)%map_size)
    enddo
    do i = 0_8, mo_integrals_tcdag_int_map%map_size
      key => mo_integrals_tcdag_int_map%map(i)%key
      val => mo_integrals_tcdag_int_map%map(i)%value
      n = mo_integrals_tcdag_int_map%map(i)%n_elements
      read(66,err=99,end=99) (key(j), j=1,n), (val(j), j=1,n)
    enddo
    call map_sort(mo_integrals_tcdag_int_map)
    load_mo_integrals_tcdag_int = 0
    return
  99 continue
    call map_deinit(mo_integrals_tcdag_int_map)
  98 continue
    stop 'Problem reading mo_integrals_tcdag_int_map file in work/'

end function load_mo_integrals_tcdag_int

