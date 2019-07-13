subroutine index_two_e_no_sym(i,j,k,l,n,i1)
  use map_module
  implicit none
 BEGIN_DOC
! symetric index for (i,k) and (j,l) but not between (i,k) and (j,l)
! n is the maximum value of the indices
 END_DOC
  integer, intent(in)            :: i,j,k,l
  integer, intent(in)            :: n
  integer(key_kind), intent(out) :: i1
  integer(key_kind)              :: p,q,r,s,i2
  p = min(i,k)
  r = max(i,k)
  p = p+ishft(r*r-r,-1)
  q = min(j,l)
  s = max(j,l)
  q = q+ishft(s*s-s,-1)
  i1 = (p-1)*n*n + q
end

subroutine index_reverse_two_e_no_sym(i,j,k,l,n,i1)
  use map_module
  implicit none
 BEGIN_DOC
! symetric index for (i,k) and (j,l) but not between (i,k) and (j,l)
! n is the maximum value of the indices
 END_DOC
  integer, intent(in)            :: n
  integer, intent(out) :: i(4),j(4),k(4),l(4)
  integer(key_kind), intent(in)  :: i1
  integer(key_kind)  :: i2
  integer*8 :: n2, p, q, r, s, p2,q2,rp,sq
  double precision  :: x
  i = 0
  n2 = n * n
  p2 = i1/n2 + 1_8
  q2 = i1 - (p2 - 1_8) * n2
  rp = (-1_8 + int(dsqrt(1.d0+8.d0 * dble(p2)),8))/2_8 
  if(rp + (rp*rp -rp)/2_8 == p2)then
   r = rp
   p = rp
  else 
   r = (1_8 + int(dsqrt(1.d0+8.d0 * dble(p2)),8))/2_8
   p = p2 - r*(r-1_8)/2_8
  endif
  sq = (-1_8 + int(dsqrt(1.d0+8.d0 * dble(q2)),8))/2_8
  if(sq + (sq*sq -sq)/2_8 == q2)then
   s = sq 
   q = sq
  else
   s = (1_8 + int(dsqrt(1.d0+8.d0 * dble(q2)),8))/2_8
   q = q2 - s*(s-1_8)/2_8
  endif
  i(1) = int(p,4)
  j(1) = int(q,4)
  k(1) = int(r,4)
  l(1) = int(s,4)


  i(2) = i(1) !ilkj
  j(2) = l(1)
  k(2) = k(1)
  l(2) = j(1)

  i(3) = k(1) !kjil
  j(3) = j(1)
  k(3) = i(1)
  l(3) = l(1)

  i(4) = k(1) !klij
  j(4) = l(1)
  k(4) = i(1)
  l(4) = j(1)

  integer :: ii, jj
  do ii=2,4
    do jj=1,ii-1
      if ( (i(ii) == i(jj)).and. &
           (j(ii) == j(jj)).and. &
           (k(ii) == k(jj)).and. &
           (l(ii) == l(jj)) ) then
         i(ii) = 0
         exit
      endif
    enddo
  enddo
  do ii=1,4
    if (i(ii) /= 0) then
      call index_two_e_no_sym(i(ii),j(ii),k(ii),l(ii),n,i2)
      if (i1 /= i2) then
        print *,  i1, i2,ii
        print *,  i(ii), j(ii), k(ii), l(ii)
        stop 'index_reverse_two_e_no_sym failed'
      endif
    endif
  enddo

end

