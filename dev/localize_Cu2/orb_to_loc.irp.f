BEGIN_PROVIDER [integer, n_orb_loc]
 implicit none
 n_orb_loc = n_ref_vec_total
END_PROVIDER 

BEGIN_PROVIDER [integer, index_orb_loc, (n_orb_loc)]
 implicit none
 integer :: i,j,k,l,index_i
 double precision, allocatable :: overlap_orb(:)
 integer , allocatable :: iorder(:)
 logical,  allocatable :: is_ok(:)
 allocate(overlap_orb(n_inact_orb),iorder(n_inact_orb),is_ok(n_inact_orb))
 is_ok = .True.
 index_orb_loc = -100000
 do i = 1, n_ref_vec_total
  print*,'i = ',i
  overlap_orb = 0.d0
  do j = 1, n_inact_orb
   index_i = list_inact(j)
   iorder(j) = j
   if(is_ok(j))then
    do l = 1, ao_num
     do k = 1, ao_num 
      overlap_orb(j) += mo_coef(k,index_i) * ref_vec_total(l,i) * ao_overlap(k,l)
     enddo
    enddo
    overlap_orb(j) = -dabs(overlap_orb(j))
   else
    overlap_orb(j) = 100.d0
   endif
  enddo
  call dsort(overlap_orb,iorder,n_inact_orb)
  do j = 1, 3
   print*,'overlap_orb = ',overlap_orb(1)
  enddo
  index_orb_loc(i) = list_inact(iorder(1))
  is_ok(iorder(1)) = .False.
 enddo
 do i = 1, n_orb_loc
  print*,'index_orb_loc',index_orb_loc(i)
 
 enddo

END_PROVIDER 
