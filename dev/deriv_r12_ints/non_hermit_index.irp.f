integer function two_int_d_dr12_index(i,j,k,l,imax,i1)
 use map_module
 implicit none
 BEGIN_DOC
! Gives a unique index for i,j,k,l using permtuation symmetry.
 ! <kl|ij> = <lk|ji>
 END_DOC
 integer, intent(in)            :: i,j,k,l,imax
 integer(key_kind), intent(out) :: i1
 integer(key_kind)              :: m_ik,n_kl,r,s,i2
! ! unique index for (i,k) using the usual matrix form
! ! m_ik = (i - 1) * imax + k
! m_ik = (i-1) * imax + k
! ! unique index for (j,l) using the usual matrix form
! ! n_kl = (j - 1) * imax + l
! n_jl = (j-1) * imax + l
! ! unique convention for (i,k) <-> (j,l)
! r = max(m_ik,n_jl)
! s = min(m_ij,n_jl) 
 i1 = 0
end
