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

