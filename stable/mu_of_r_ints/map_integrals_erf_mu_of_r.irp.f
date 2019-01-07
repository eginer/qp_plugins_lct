use map_module


BEGIN_TEMPLATE

subroutine dump_$ao_integrals(filename)
  use map_module
  implicit none
  BEGIN_DOC
  ! Save to disk the $ao integrals
  END_DOC
  character*(*), intent(in)      :: filename
  integer(cache_key_kind), pointer :: key(:)
  real(integral_kind), pointer   :: val(:)
  integer*8                      :: i,j, n
  call ezfio_set_work_empty(.False.)
  open(unit=66,file=filename,FORM='unformatted')
  write(66) integral_kind, key_kind
  write(66) $ao_integrals_map%sorted, $ao_integrals_map%map_size,    &
      $ao_integrals_map%n_elements
  do i=0_8,$ao_integrals_map%map_size
    write(66) $ao_integrals_map%map(i)%sorted, $ao_integrals_map%map(i)%map_size,&
        $ao_integrals_map%map(i)%n_elements
  enddo
  do i=0_8,$ao_integrals_map%map_size
    key => $ao_integrals_map%map(i)%key
    val => $ao_integrals_map%map(i)%value
    n = $ao_integrals_map%map(i)%n_elements
    write(66) (key(j), j=1,n), (val(j), j=1,n)
  enddo
  close(66)
  
end

IRP_IF COARRAY
subroutine communicate_$ao_integrals()
  use map_module
  implicit none
  BEGIN_DOC
  ! Communicate the $ao integrals with co-array
  END_DOC
  integer(cache_key_kind), pointer :: key(:)
  real(integral_kind), pointer   :: val(:)
  integer*8                      :: i,j, k, nmax
  integer*8, save                :: n[*]
  integer                        :: copy_n

  real(integral_kind), allocatable            :: buffer_val(:)[:]
  integer(cache_key_kind), allocatable        :: buffer_key(:)[:]
  real(integral_kind), allocatable            :: copy_val(:)
  integer(key_kind), allocatable              :: copy_key(:)

  n = 0_8
  do i=0_8,$ao_integrals_map%map_size
    n = max(n,$ao_integrals_map%map(i)%n_elements)
  enddo
  sync all
  nmax = 0_8
  do j=1,num_images()
    nmax = max(nmax,n[j])
  enddo
  allocate( buffer_key(nmax)[*], buffer_val(nmax)[*])
  allocate( copy_key(nmax), copy_val(nmax))
  do i=0_8,$ao_integrals_map%map_size
    key => $ao_integrals_map%map(i)%key
    val => $ao_integrals_map%map(i)%value
    n = $ao_integrals_map%map(i)%n_elements
    do j=1,n
      buffer_key(j) = key(j)
      buffer_val(j) = val(j)
    enddo
    sync all
    do j=1,num_images()
      if (j /= this_image()) then
        copy_n = n[j]
        do k=1,copy_n
          copy_val(k) = buffer_val(k)[j]
          copy_key(k) = buffer_key(k)[j]
          copy_key(k) = copy_key(k)+ishft(i,-map_shift)
        enddo
        call map_append($ao_integrals_map, copy_key, copy_val, copy_n )
      endif
    enddo
    sync all
  enddo
  deallocate( buffer_key, buffer_val, copy_val, copy_key)
  
end
IRP_ENDIF 


integer function load_$ao_integrals(filename)
  implicit none
  BEGIN_DOC
  ! Read from disk the $ao integrals
  END_DOC
  character*(*), intent(in)      :: filename
  integer*8                      :: i
  integer(cache_key_kind), pointer :: key(:)
  real(integral_kind), pointer   :: val(:)
  integer                        :: iknd, kknd
  integer*8                      :: n, j
  load_$ao_integrals = 1
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
  read(66,err=98,end=98) $ao_integrals_map%sorted, $ao_integrals_map%map_size,&
      $ao_integrals_map%n_elements
  do i=0_8, $ao_integrals_map%map_size
    read(66,err=99,end=99) $ao_integrals_map%map(i)%sorted,          &
        $ao_integrals_map%map(i)%map_size, $ao_integrals_map%map(i)%n_elements
    call cache_map_reallocate($ao_integrals_map%map(i),$ao_integrals_map%map(i)%map_size)
  enddo
  do i=0_8, $ao_integrals_map%map_size
    key => $ao_integrals_map%map(i)%key
    val => $ao_integrals_map%map(i)%value
    n = $ao_integrals_map%map(i)%n_elements
    read(66,err=99,end=99) (key(j), j=1,n), (val(j), j=1,n)
  enddo
  call map_sort($ao_integrals_map)
  load_$ao_integrals = 0
  return
  99 continue
  call map_deinit($ao_integrals_map)
  98 continue
  stop 'Problem reading $ao_integrals_map file in work/'
  
end

SUBST [ ao_integrals_map, ao_integrals, ao_num ]
ao_integrals_erf_mu_of_r_map ; ao_integrals_erf_mu_of_r ; ao_num ;;
mo_integrals_erf_mu_of_r_map ; mo_integrals_erf_mu_of_r ; mo_tot_num;;
END_TEMPLATE




BEGIN_PROVIDER [ type(map_type), mo_integrals_erf_mu_of_r_map ]
  implicit none
  BEGIN_DOC
  ! MO integrals
  END_DOC
  integer(key_kind)              :: key_max
  integer(map_size_kind)         :: sze
  call bielec_integrals_index(mo_tot_num,mo_tot_num,mo_tot_num,mo_tot_num,key_max)
  sze = key_max
  call map_init(mo_integrals_erf_mu_of_r_map,sze)
  print*, 'MO ERF map initialized'
END_PROVIDER

subroutine insert_into_ao_integrals_erf_mu_of_r_map(n_integrals,buffer_i, buffer_values)
  use map_module
  implicit none
  BEGIN_DOC
  ! Create new entry into AO map
  END_DOC
  
  integer, intent(in)                :: n_integrals
  integer(key_kind), intent(inout)   :: buffer_i(n_integrals)
  real(integral_kind), intent(inout) :: buffer_values(n_integrals)
  
  call map_append(ao_integrals_erf_mu_of_r_map, buffer_i, buffer_values, n_integrals)
end

subroutine insert_into_mo_integrals_erf_mu_of_r_map(n_integrals,                 &
      buffer_i, buffer_values, thr)
  use map_module
  implicit none
  
  BEGIN_DOC
  ! Create new entry into MO map, or accumulate in an existing entry
  END_DOC
  
  integer, intent(in)                :: n_integrals
  integer(key_kind), intent(inout)   :: buffer_i(n_integrals)
  real(integral_kind), intent(inout) :: buffer_values(n_integrals)
  real(integral_kind), intent(in)    :: thr
  call map_update(mo_integrals_erf_mu_of_r_map, buffer_i, buffer_values, n_integrals, thr)
end

 BEGIN_PROVIDER [ integer, mo_integrals_erf_mu_of_r_cache_min ]
&BEGIN_PROVIDER [ integer, mo_integrals_erf_mu_of_r_cache_max ]
 implicit none
 BEGIN_DOC
 ! Min and max values of the MOs for which the integrals are in the cache
 END_DOC
 mo_integrals_erf_mu_of_r_cache_min = max(1,elec_alpha_num - 31)
 mo_integrals_erf_mu_of_r_cache_max = min(mo_tot_num,mo_integrals_erf_mu_of_r_cache_min+63)

END_PROVIDER

BEGIN_PROVIDER [ double precision, mo_integrals_erf_mu_of_r_cache, (0:64*64*64*64) ]
 implicit none
 BEGIN_DOC
 ! Cache of MO integrals for fast access
 END_DOC
 PROVIDE mo_bielec_integrals_erf_mu_of_r_in_map
 integer                        :: i,j,k,l
 integer                        :: ii
 integer(key_kind)              :: idx
 real(integral_kind)            :: integral
!FREE ao_integrals_erf_mu_of_r_cache
 !$OMP PARALLEL DO PRIVATE (i,j,k,l,idx,ii,integral)
 do l=mo_integrals_erf_mu_of_r_cache_min,mo_integrals_erf_mu_of_r_cache_max
   do k=mo_integrals_erf_mu_of_r_cache_min,mo_integrals_erf_mu_of_r_cache_max
     do j=mo_integrals_erf_mu_of_r_cache_min,mo_integrals_erf_mu_of_r_cache_max
       do i=mo_integrals_erf_mu_of_r_cache_min,mo_integrals_erf_mu_of_r_cache_max
         !DIR$ FORCEINLINE
         call bielec_integrals_index(i,j,k,l,idx)
         !DIR$ FORCEINLINE
         call map_get(mo_integrals_erf_mu_of_r_map,idx,integral)
         ii = l-mo_integrals_erf_mu_of_r_cache_min
         ii = ior( ishft(ii,6), k-mo_integrals_erf_mu_of_r_cache_min)
         ii = ior( ishft(ii,6), j-mo_integrals_erf_mu_of_r_cache_min)
         ii = ior( ishft(ii,6), i-mo_integrals_erf_mu_of_r_cache_min)
         mo_integrals_erf_mu_of_r_cache(ii) = integral
       enddo
     enddo
   enddo
 enddo
 !$OMP END PARALLEL DO

END_PROVIDER


double precision function get_mo_bielec_integral_erf_mu_of_r(i,j,k,l,map)
  use map_module
  implicit none
  BEGIN_DOC
  ! Returns one integral <ij|kl> in the MO basis
  END_DOC
  integer, intent(in)            :: i,j,k,l
  integer(key_kind)              :: idx
  integer                        :: ii
  type(map_type), intent(inout)  :: map
  real(integral_kind)            :: tmp
  PROVIDE mo_bielec_integrals_erf_mu_of_r_in_map mo_integrals_erf_mu_of_r_cache
  ii = l-mo_integrals_erf_mu_of_r_cache_min
  ii = ior(ii, k-mo_integrals_erf_mu_of_r_cache_min)
  ii = ior(ii, j-mo_integrals_erf_mu_of_r_cache_min)
  ii = ior(ii, i-mo_integrals_erf_mu_of_r_cache_min)
  if (iand(ii, -64) /= 0) then
    !DIR$ FORCEINLINE
    call bielec_integrals_index(i,j,k,l,idx)
    !DIR$ FORCEINLINE
    call map_get(map,idx,tmp)
    get_mo_bielec_integral_erf_mu_of_r = dble(tmp)
  else
    ii = l-mo_integrals_erf_mu_of_r_cache_min
    ii = ior( ishft(ii,6), k-mo_integrals_erf_mu_of_r_cache_min)
    ii = ior( ishft(ii,6), j-mo_integrals_erf_mu_of_r_cache_min)
    ii = ior( ishft(ii,6), i-mo_integrals_erf_mu_of_r_cache_min)
    get_mo_bielec_integral_erf_mu_of_r = mo_integrals_erf_mu_of_r_cache(ii)
  endif
end


double precision function mo_bielec_integral_erf_mu_of_r(i,j,k,l)
  implicit none
  BEGIN_DOC
  ! Returns one integral <ij|kl> in the MO basis
  END_DOC
  integer, intent(in)            :: i,j,k,l
  double precision               :: get_mo_bielec_integral_erf_mu_of_r
  PROVIDE mo_bielec_integrals_erf_mu_of_r_in_map mo_integrals_erf_mu_of_r_cache
  !DIR$ FORCEINLINE
  PROVIDE mo_bielec_integrals_erf_mu_of_r_in_map
  mo_bielec_integral_erf_mu_of_r = get_mo_bielec_integral_erf_mu_of_r(i,j,k,l,mo_integrals_erf_mu_of_r_map)
  return
end

subroutine get_mo_bielec_integrals_erf_mu_of_r(j,k,l,sze,out_val,map)
  use map_module
  implicit none
  BEGIN_DOC
  ! Returns multiple integrals <ij|kl> in the MO basis, all
  ! i for j,k,l fixed.
  END_DOC
  integer, intent(in)            :: j,k,l, sze
  double precision, intent(out)  :: out_val(sze)
  type(map_type), intent(inout)  :: map
  integer                        :: i
  integer(key_kind)              :: hash(sze)
  real(integral_kind)            :: tmp_val(sze)
  PROVIDE mo_bielec_integrals_erf_mu_of_r_in_map
  
  do i=1,sze
    !DIR$ FORCEINLINE
    call bielec_integrals_index(i,j,k,l,hash(i))
  enddo
  
  if (key_kind == 8) then
    call map_get_many(map, hash, out_val, sze)
  else
    call map_get_many(map, hash, tmp_val, sze)
    ! Conversion to double precision 
    do i=1,sze
      out_val(i) = dble(tmp_val(i))
    enddo
  endif
end

subroutine get_mo_bielec_integrals_erf_mu_of_r_ij(k,l,sze,out_array,map)
  use map_module
  implicit none
  BEGIN_DOC
  ! Returns multiple integrals <ij|kl> in the MO basis, all
  ! i(1)j(2) 1/r12 k(1)l(2)
  ! i, j for k,l fixed.
  END_DOC
  integer, intent(in)            :: k,l, sze
  double precision, intent(out)  :: out_array(sze,sze)
  type(map_type), intent(inout)  :: map
  integer                        :: i,j,kk,ll,m
  integer(key_kind),allocatable  :: hash(:)
  integer  ,allocatable          :: pairs(:,:), iorder(:)
  real(integral_kind), allocatable :: tmp_val(:)

  PROVIDE mo_bielec_integrals_erf_mu_of_r_in_map
  allocate (hash(sze*sze), pairs(2,sze*sze),iorder(sze*sze), &
  tmp_val(sze*sze))
  
  kk=0
  out_array = 0.d0
  do j=1,sze
   do i=1,sze
    kk += 1
    !DIR$ FORCEINLINE
    call bielec_integrals_index(i,j,k,l,hash(kk))
    pairs(1,kk) = i
    pairs(2,kk) = j
    iorder(kk) = kk
   enddo
  enddo

  logical :: integral_is_in_map
  if (key_kind == 8) then
    call i8radix_sort(hash,iorder,kk,-1)
  else if (key_kind == 4) then
    call iradix_sort(hash,iorder,kk,-1)
  else if (key_kind == 2) then
    call i2radix_sort(hash,iorder,kk,-1)
  endif

  call map_get_many(mo_integrals_erf_mu_of_r_map, hash, tmp_val, kk)

  do ll=1,kk
    m = iorder(ll)
    i=pairs(1,m)
    j=pairs(2,m)
    out_array(i,j) = tmp_val(ll)
  enddo  

  deallocate(pairs,hash,iorder,tmp_val)
end

subroutine get_mo_bielec_integrals_erf_mu_of_r_coulomb_ii(k,l,sze,out_val,map)
  use map_module
  implicit none
  BEGIN_DOC
  ! Returns multiple integrals <ki|li> 
  ! k(1)i(2) 1/r12 l(1)i(2) :: out_val(i1)
  ! for k,l fixed.
  END_DOC
  integer, intent(in)            :: k,l, sze
  double precision, intent(out)  :: out_val(sze)
  type(map_type), intent(inout)  :: map
  integer                        :: i
  integer(key_kind)              :: hash(sze)
  real(integral_kind)            :: tmp_val(sze)
  PROVIDE mo_bielec_integrals_erf_mu_of_r_in_map
  
  integer :: kk
  do i=1,sze
    !DIR$ FORCEINLINE
    call bielec_integrals_index(k,i,l,i,hash(i))
  enddo
  
  if (key_kind == 8) then
    call map_get_many(map, hash, out_val, sze)
  else
    call map_get_many(map, hash, tmp_val, sze)
    ! Conversion to double precision 
    do i=1,sze
      out_val(i) = dble(tmp_val(i))
    enddo
  endif
end

subroutine get_mo_bielec_integrals_erf_mu_of_r_exch_ii(k,l,sze,out_val,map)
  use map_module
  implicit none
  BEGIN_DOC
  ! Returns multiple integrals <ki|il> 
  ! k(1)i(2) 1/r12 i(1)l(2) :: out_val(i1)
  ! for k,l fixed.
  END_DOC
  integer, intent(in)            :: k,l, sze
  double precision, intent(out)  :: out_val(sze)
  type(map_type), intent(inout)  :: map
  integer                        :: i
  integer(key_kind)              :: hash(sze)
  real(integral_kind)            :: tmp_val(sze)
  PROVIDE mo_bielec_integrals_erf_mu_of_r_in_map
  
  integer :: kk
  do i=1,sze
    !DIR$ FORCEINLINE
    call bielec_integrals_index(k,i,i,l,hash(i))
  enddo
  
  if (key_kind == 8) then
    call map_get_many(map, hash, out_val, sze)
  else
    call map_get_many(map, hash, tmp_val, sze)
    ! Conversion to double precision 
    do i=1,sze
      out_val(i) = dble(tmp_val(i))
    enddo
  endif
end


integer*8 function get_mo_erf_mu_of_r_map_size()
  implicit none
  BEGIN_DOC
  ! Return the number of elements in the MO map
  END_DOC
  get_mo_erf_mu_of_r_map_size = mo_integrals_erf_mu_of_r_map % n_elements
end
