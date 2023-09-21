
use map_module

!! MO Map
!! ======

BEGIN_PROVIDER [ type(map_type), mo_int_erf_mu_of_r_map ]
  implicit none
  BEGIN_DOC
  ! MO integrals
  END_DOC
  integer(key_kind)              :: key_max
  integer(map_size_kind)         :: sze
  call two_e_integrals_index(mo_num,mo_num,mo_num,mo_num,key_max)
  sze = key_max
  call map_init(mo_int_erf_mu_of_r_map,sze)
  print*,  'MO map initialized : ', sze
END_PROVIDER

double precision function get_mo_two_e_int_erf_mu_of_r(i,j,k,l,map) result(result)
  use map_module
  implicit none
  BEGIN_DOC
  ! Gets one MO bi-electronic integral from the MO map
  ! (i,k) = 1, (j,l) = 2
  END_DOC
  integer, intent(in)            :: i,j,k,l
  integer(key_kind)              :: idx
  type(map_type), intent(inout)  :: map
  integer                        :: ii
  real(integral_kind)            :: tmp
  PROVIDE mo_int_erf_mu_of_r_in_map 
  !DIR$ FORCEINLINE
!  if (mo_overlap_abs(i,k)*mo_overlap_abs(j,l) < mo_integrals_threshold ) then
!    tmp = 0.d0
!  else
      !DIR$ FORCEINLINE
      call two_e_integrals_index(i,j,k,l,idx)
      !DIR$ FORCEINLINE
      call map_get(map,idx,tmp)
!  endif
  result = tmp
end



subroutine get_mo_two_e_ints_erf_mu_of_r(j,k,l,sze,out_val)
  use map_module
  BEGIN_DOC
  ! Gets multiple MO bi-electronic integral from the MO map .
  ! All i are retrieved for j,k,l fixed.
  ! (j,l) = 1 ; (k,:) = 2
  END_DOC
  implicit none
  integer, intent(in)            :: j,k,l, sze
  real(integral_kind), intent(out) :: out_val(sze)
  
  integer                        :: i
  integer(key_kind)              :: hash
  double precision               :: thresh
  PROVIDE mo_int_erf_mu_of_r_in_map mo_int_erf_mu_of_r_map
  thresh = mo_integrals_threshold
  
!  if (mo_overlap_abs(j,l) < thresh) then
!    out_val = 0.d0
!    return
!  endif
  
  double precision :: get_mo_two_e_int_erf_mu_of_r
  do i=1,sze
    out_val(i) = get_mo_two_e_int_erf_mu_of_r(i,j,k,l,mo_int_erf_mu_of_r_map) 
  enddo
  
end

subroutine get_mo_two_e_ints_erf_mu_of_r_non_zero(j,k,l,sze,out_val,out_val_index,non_zero_int)
  use map_module
  implicit none
  BEGIN_DOC
  ! Gets multiple MO bi-electronic integral from the MO map .
  ! All non-zero i are retrieved for j,k,l fixed.
  END_DOC
  integer, intent(in)            :: j,k,l, sze
  real(integral_kind), intent(out) :: out_val(sze)
  integer, intent(out)           :: out_val_index(sze),non_zero_int
  
  integer                        :: i
  integer(key_kind)              :: hash
  double precision               :: thresh,tmp
  PROVIDE mo_int_erf_mu_of_r_in_map
  thresh = mo_integrals_threshold
  
  non_zero_int = 0
!  if (mo_overlap_abs(j,l) < thresh) then
!    out_val = 0.d0
!    return
!  endif
 
  non_zero_int = 0
  do i=1,sze
    integer, external :: mo_l4
    double precision, external :: mo_two_e_int_erf_mu_of_r
    !DIR$ FORCEINLINE
!    if (mo_overlap_abs(i,k)*mo_overlap_abs(j,l) < thresh) then
!      cycle
!    endif
    call two_e_integrals_index(i,j,k,l,hash)
    call map_get(mo_int_erf_mu_of_r_map, hash,tmp)
    if (dabs(tmp) < thresh ) cycle
     non_zero_int = non_zero_int+1
     out_val_index(non_zero_int) = i
     out_val(non_zero_int) = tmp
  enddo
  
end


function get_mo_erf_mu_of_r_map_size()
  implicit none
  integer (map_size_kind) :: get_mo_erf_mu_of_r_map_size
  BEGIN_DOC
  ! Returns the number of elements in the MO map
  END_DOC
  get_mo_erf_mu_of_r_map_size = mo_int_erf_mu_of_r_map % n_elements
end

subroutine clear_mo_erf_mu_of_r_map
  implicit none
  BEGIN_DOC
  ! Frees the memory of the MO map
  END_DOC
  call map_deinit(mo_int_erf_mu_of_r_map)
  FREE mo_int_erf_mu_of_r_map
end


subroutine insert_into_mo_int_erf_mu_of_r_map(n_integrals,buffer_i, buffer_values)
  use map_module
  implicit none
  BEGIN_DOC
  ! Create new entry into AO map
  END_DOC
  
  integer, intent(in)                :: n_integrals
  integer(key_kind), intent(inout)   :: buffer_i(n_integrals)
  real(integral_kind), intent(inout) :: buffer_values(n_integrals)
  
  call map_append(mo_int_erf_mu_of_r_map, buffer_i, buffer_values, n_integrals)
end

