use map_module

!! mo Map
!! ======

BEGIN_PROVIDER [ type(map_type), mo_integrals_tc_int_map ]
  implicit none
  BEGIN_DOC
  ! |mo| integrals
  END_DOC
  return
END_PROVIDER

 BEGIN_PROVIDER [ integer, mo_integrals_tc_int_cache_min ]
&BEGIN_PROVIDER [ integer, mo_integrals_tc_int_cache_max ]
 implicit none
 BEGIN_DOC
 ! Min and max values of the mos for which the integrals are in the cache
 END_DOC
 mo_integrals_tc_int_cache_min = max(1,mo_num - 63)
 mo_integrals_tc_int_cache_max = mo_num

END_PROVIDER

BEGIN_PROVIDER [ double precision, mo_integrals_tc_int_cache, (0:64*64*64*64) ]
  use map_module
 implicit none
 BEGIN_DOC
 ! Cache of |mo| integrals for fast access
 END_DOC
 PROVIDE mo_two_e_integrals_tc_int_in_map
 integer                        :: i,j,k,l,ii
 integer(key_kind)              :: idx
 real(integral_kind)            :: integral
 !$OMP PARALLEL DO PRIVATE (i,j,k,l,idx,ii,integral)
 do l=mo_integrals_tc_int_cache_min,mo_integrals_tc_int_cache_max
   do k=mo_integrals_tc_int_cache_min,mo_integrals_tc_int_cache_max
     do j=mo_integrals_tc_int_cache_min,mo_integrals_tc_int_cache_max
       do i=mo_integrals_tc_int_cache_min,mo_integrals_tc_int_cache_max
         !DIR$ FORCEINLINE
         call two_e_integrals_index(i,j,k,l,idx)
         !DIR$ FORCEINLINE
         call map_get(mo_integrals_tc_int_map,idx,integral)
         ii = l-mo_integrals_tc_int_cache_min
         ii = ior( ishft(ii,6), k-mo_integrals_tc_int_cache_min)
         ii = ior( ishft(ii,6), j-mo_integrals_tc_int_cache_min)
         ii = ior( ishft(ii,6), i-mo_integrals_tc_int_cache_min)
         mo_integrals_tc_int_cache(ii) = integral
       enddo
     enddo
   enddo
 enddo
 !$OMP END PARALLEL DO

END_PROVIDER


subroutine insert_into_mo_integrals_tc_int_map(n_integrals,buffer_i, buffer_values)
  use map_module
  implicit none
  BEGIN_DOC
  ! Create new entry into |mo| map
  END_DOC

  integer, intent(in)                :: n_integrals
  integer(key_kind), intent(inout)   :: buffer_i(n_integrals)
  real(integral_kind), intent(inout) :: buffer_values(n_integrals)

  call map_append(mo_integrals_tc_int_map, buffer_i, buffer_values, n_integrals)
end

double precision function get_mo_two_e_integral_tc_int(i,j,k,l,map) result(result)
  use map_module
  implicit none
  BEGIN_DOC
  ! Gets one |mo| two-electron integral from the |mo| map
  END_DOC
  integer, intent(in)            :: i,j,k,l
  integer(key_kind)              :: idx
  type(map_type), intent(inout)  :: map
  integer                        :: ii
  real(integral_kind)            :: tmp
  logical, external              :: mo_two_e_integral_zero
  PROVIDE mo_two_e_integrals_tc_int_in_map mo_integrals_tc_int_cache mo_integrals_tc_int_cache_min
  !DIR$ FORCEINLINE
!  if (mo_two_e_integral_zero(i,j,k,l)) then
!  if (.False.) then
!    tmp = 0.d0
!  else if (mo_two_e_integral_tc_int_schwartz(i,k)*mo_two_e_integral_tc_int_schwartz(j,l) < mo_integrals_threshold) then
!    tmp = 0.d0
!  else
    ii = l-mo_integrals_tc_int_cache_min
    ii = ior(ii, k-mo_integrals_tc_int_cache_min)
    ii = ior(ii, j-mo_integrals_tc_int_cache_min)
    ii = ior(ii, i-mo_integrals_tc_int_cache_min)
    if (iand(ii, -64) /= 0) then
      !DIR$ FORCEINLINE
      call two_e_integrals_index(i,j,k,l,idx)
      !DIR$ FORCEINLINE
      call map_get(map,idx,tmp)
      tmp = tmp
    else
      ii = l-mo_integrals_tc_int_cache_min
      ii = ior( ishft(ii,6), k-mo_integrals_tc_int_cache_min)
      ii = ior( ishft(ii,6), j-mo_integrals_tc_int_cache_min)
      ii = ior( ishft(ii,6), i-mo_integrals_tc_int_cache_min)
      tmp = mo_integrals_tc_int_cache(ii)
    endif
!  endif
  result = tmp
end


subroutine get_mo_two_e_integrals_tc_int(j,k,l,sze,out_val)
  use map_module
  BEGIN_DOC
  ! Gets multiple |mo| two-electron integral from the |mo| map .
  ! All i are retrieved for j,k,l fixed.
  END_DOC
  implicit none
  integer, intent(in)            :: j,k,l, sze
  real(integral_kind), intent(out) :: out_val(sze)

  integer                        :: i
  integer(key_kind)              :: hash
  double precision               :: thresh
! logical, external              :: mo_one_e_integral_zero
  PROVIDE mo_two_e_integrals_tc_int_in_map mo_integrals_tc_int_map
  thresh = mo_integrals_threshold

! if (mo_one_e_integral_zero(j,l)) then
  if (.False.) then
    out_val = 0.d0
    return
  endif

  double precision :: get_mo_two_e_integral_tc_int
  do i=1,sze
    out_val(i) = get_mo_two_e_integral_tc_int(i,j,k,l,mo_integrals_tc_int_map)
  enddo

end

subroutine get_mo_two_e_integrals_tc_int_non_zero(j,k,l,sze,out_val,out_val_index,non_zero_int)
  use map_module
  implicit none
  BEGIN_DOC
  ! Gets multiple |mo| two-electron integrals from the |mo| map .
  ! All non-zero i are retrieved for j,k,l fixed.
  END_DOC
  integer, intent(in)            :: j,k,l, sze
  real(integral_kind), intent(out) :: out_val(sze)
  integer, intent(out)           :: out_val_index(sze),non_zero_int

  integer                        :: i
  integer(key_kind)              :: hash
  double precision               :: thresh,tmp
! logical, external              :: mo_one_e_integral_zero
  PROVIDE mo_two_e_integrals_tc_int_in_map
  thresh = mo_integrals_threshold

  non_zero_int = 0
! if (mo_one_e_integral_zero(j,l)) then
  if (.False.) then
    out_val = 0.d0
    return
  endif

  non_zero_int = 0
  do i=1,sze
    integer, external :: mo_l4
!    double precision, external :: mo_two_e_integral_tc_int
    !DIR$ FORCEINLINE
!    if (mo_two_e_integral_tc_int_schwartz(i,k)*mo_two_e_integral_tc_int_schwartz(j,l) < thresh) then
!      cycle
!    endif
    call two_e_integrals_index(i,j,k,l,hash)
    call map_get(mo_integrals_tc_int_map, hash,tmp)
    if (dabs(tmp) < thresh ) cycle
    non_zero_int = non_zero_int+1
    out_val_index(non_zero_int) = i
    out_val(non_zero_int) = tmp
  enddo

end


function get_mo_tc_int_map_size()
  implicit none
  integer (map_size_kind) :: get_mo_tc_int_map_size
  BEGIN_DOC
  ! Returns the number of elements in the |mo| map
  END_DOC
  get_mo_tc_int_map_size = mo_integrals_tc_int_map % n_elements
end

subroutine clear_mo_tc_int_map
  implicit none
  BEGIN_DOC
  ! Frees the memory of the |mo| map
  END_DOC
  call map_deinit(mo_integrals_tc_int_map)
  FREE mo_integrals_tc_int_map
end



subroutine dump_mo_integrals_tc_int(filename)
  use map_module
  implicit none
  BEGIN_DOC
  ! Save to disk the |mo| tc_int integrals
  END_DOC
  character*(*), intent(in)      :: filename
  integer(cache_key_kind), pointer :: key(:)
  real(integral_kind), pointer   :: val(:)
  integer*8                      :: i,j, n
  call ezfio_set_work_empty(.False.)
  open(unit=66,file=filename,FORM='unformatted')
  write(66) integral_kind, key_kind
  write(66) mo_integrals_tc_int_map%sorted, mo_integrals_tc_int_map%map_size,    &
      mo_integrals_tc_int_map%n_elements
  do i=0_8,mo_integrals_tc_int_map%map_size
    write(66) mo_integrals_tc_int_map%map(i)%sorted, mo_integrals_tc_int_map%map(i)%map_size,&
        mo_integrals_tc_int_map%map(i)%n_elements
  enddo
  do i=0_8,mo_integrals_tc_int_map%map_size
    key => mo_integrals_tc_int_map%map(i)%key
    val => mo_integrals_tc_int_map%map(i)%value
    n = mo_integrals_tc_int_map%map(i)%n_elements
    write(66) (key(j), j=1,n), (val(j), j=1,n)
  enddo
  close(66)

end



integer function load_mo_integrals_tc_int(filename)
  implicit none
  BEGIN_DOC
  ! Read from disk the |mo| tc_int integrals
  END_DOC
  character*(*), intent(in)      :: filename
  integer*8                      :: i
  integer(cache_key_kind), pointer :: key(:)
  real(integral_kind), pointer   :: val(:)
  integer                        :: iknd, kknd
  integer*8                      :: n, j
  load_mo_integrals_tc_int = 1
  open(unit=66,file=filename,FORM='unformatted',STATUS='UNKNOWN')
  read(66,err=98,end=98) iknd, kknd
  if (iknd /= integral_kind) then
    print *,  'Wrong integrals kind in file :', iknd
    stop 1
  endif
  if (kknd /= key_kind) then
    print *,  'Wrong key kind in file :', kknd
    stop 1
  endif
  read(66,err=98,end=98) mo_integrals_tc_int_map%sorted, mo_integrals_tc_int_map%map_size,&
      mo_integrals_tc_int_map%n_elements
  do i=0_8, mo_integrals_tc_int_map%map_size
    read(66,err=99,end=99) mo_integrals_tc_int_map%map(i)%sorted,          &
        mo_integrals_tc_int_map%map(i)%map_size, mo_integrals_tc_int_map%map(i)%n_elements
    call cache_map_reallocate(mo_integrals_tc_int_map%map(i),mo_integrals_tc_int_map%map(i)%map_size)
  enddo
  do i=0_8, mo_integrals_tc_int_map%map_size
    key => mo_integrals_tc_int_map%map(i)%key
    val => mo_integrals_tc_int_map%map(i)%value
    n = mo_integrals_tc_int_map%map(i)%n_elements
    read(66,err=99,end=99) (key(j), j=1,n), (val(j), j=1,n)
  enddo
  call map_sort(mo_integrals_tc_int_map)
  load_mo_integrals_tc_int = 0
  return
  99 continue
  call map_deinit(mo_integrals_tc_int_map)
  98 continue
  stop 'Problem reading mo_integrals_tc_int_map file in work/'

end




