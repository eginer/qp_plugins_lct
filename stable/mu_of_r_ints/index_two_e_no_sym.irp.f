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

subroutine index_reverse_two_e_no_sym(p,q,r,s,n,i1)
  use map_module
  implicit none
 BEGIN_DOC
! symetric index for (i,k) and (j,l) but not between (i,k) and (j,l)
! n is the maximum value of the indices
 END_DOC
  integer, intent(in)            :: n
  integer(key_kind), intent(in)  :: i1
  integer(key_kind), intent(out) :: i(4),j(4),j(4),l(4)
  double precision  :: x
! n2 = n * n
! p = i1/n2 + 1
! q = i1 - n2 - p*n2
! r = 1 + int(dsqrt(1+8 * p))/2
! s = p - r(r-1)/2
  integer*8 :: n2, p, q, r, s                                                              
  integer :: i(4),j(4)                                                                         
  n2 = int(n*n,8)                                                                              
  p = i1/n2 + 1_8                                                                              
  q = i1-n2 - p*n2                                                                             
  x = dble(1_8+shiftl(p,3))                                                                    
  r = 1_8 + shiftr(int(dsqrt(x)),1)                                                            
  s = p - shiftr(r*r-r,1)
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
        print *,  i1, i2
        print *,  i(ii), j(ii), k(ii), l(ii)
        stop 'two_e_integrals_index_reverse failed'
      endif
    endif
  enddo

end

