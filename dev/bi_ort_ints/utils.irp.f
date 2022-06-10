
! ---

subroutine rotate_mo_coef(mo_coef_old, mo_mo_coef, mo_mo_overlap, mo_coef_new)

  BEGIN_DOC
  ! You have mo_coef_new which is based on a MO->MO transformation through mo_mo_coef
  END_DOC

  implicit none
  double precision, intent(in)  :: mo_coef_old(ao_num, mo_num)
  double precision, intent(in)  :: mo_mo_coef(mo_num, mo_num), mo_mo_overlap(mo_num, mo_num)
  double precision, intent(out) :: mo_coef_new(ao_num, mo_num)

! call dgemm('N','N',ao_num,mo_num,mo_num,1.d0,mo_coef_old,size(mo_coef_old,1),& 
!      mo_mo_coef,size(mo_mo_coef,1),&
!      0.d0,mo_coef_new,size(mo_coef_new,1))

  integer :: i, j, q, k
  mo_coef_new = 0.d0
  do i = 1, mo_num
    do j = 1, mo_num
      do k = 1, mo_num
        do q = 1, ao_num
          mo_coef_new(q,i) += mo_coef_old(q,j) * mo_mo_overlap(k,j) * mo_mo_coef(j,i) 
        enddo
      enddo
    enddo
  enddo

end subroutine rotate_mo_coef

! ---

subroutine ao_to_mo_bi_ortho(A_ao, LDA_ao, A_mo, LDA_mo)

  BEGIN_DOC
  ! Transform A from the |AO| basis to the BI ORTHONORMAL MOS 
  !
  ! $C_L^\dagger.A_{ao}.C_R$ where C_L and C_R are the LEFT and RIGHT MO coefs
  END_DOC

  implicit none
  integer, intent(in)           :: LDA_ao,LDA_mo
  double precision, intent(in)  :: A_ao(LDA_ao,ao_num)
  double precision, intent(out) :: A_mo(LDA_mo,mo_num)
  double precision, allocatable :: T(:,:)

  allocate ( T(ao_num,mo_num) )
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: T

  call dgemm('N', 'N', ao_num, mo_num, ao_num,  &
      1.d0, A_ao, LDA_ao,                       &
      mo_r_coef, size(mo_r_coef, 1),            &
      0.d0, T, size(T, 1))

  call dgemm('T', 'N', mo_num, mo_num, ao_num, &
      1.d0, mo_l_coef, size(mo_l_coef, 1),     &
      T, ao_num,                               &
      0.d0, A_mo, size(A_mo, 1))

!  call restore_symmetry(mo_num,mo_num,A_mo,size(A_mo,1),1.d-12)
  deallocate(T)

end subroutine ao_to_mo_bi_ortho

! ---

!subroutine two_e_integrals_index_reverse_tc(i,j,k,l,i1)
!  use map_module
!  implicit none
!  BEGIN_DOC
!! Computes the 4 indices $i,j,k,l$ from a unique index $i_1$.
!! For 2 indices $i,j$ and $i \le j$, we have
!! $p = i(i-1)/2 + j$.
!! The key point is that because $j < i$,
!! $i(i-1)/2 < p \le i(i+1)/2$. So $i$ can be found by solving
!! $i^2 - i - 2p=0$. One obtains $i=1 + \sqrt{1+8p}/2$
!! and $j = p - i(i-1)/2$.
!! This rule is applied 3 times. First for the symmetry of the
!! pairs (i,k) and (j,l), and then for the symmetry within each pair.
!  END_DOC
!  integer, intent(out)           :: i(2),j(2),k(2),l(2)
!  integer(key_kind), intent(in)  :: i1
!  integer(key_kind)              :: i2,i3
!
!  i = 0
!  i2   = ceiling(0.5d0*(dsqrt(dble(shiftl(i1,3)+1))-1.d0))
!  l(1) = ceiling(0.5d0*(dsqrt(dble(shiftl(i2,3)+1))-1.d0))
!  i3   = i1 - shiftr(i2*i2-i2,1)
!  k(1) = ceiling(0.5d0*(dsqrt(dble(shiftl(i3,3)+1))-1.d0))
!  j(1) = int(i2 - shiftr(l(1)*l(1)-l(1),1),4)
!  i(1) = int(i3 - shiftr(k(1)*k(1)-k(1),1),4)
!
!              !ijkl
!  i(2) = j(1) !jilk
!  j(2) = i(1)
!  k(2) = l(1)
!  l(2) = k(1)
!
!  integer :: ii, jj
!  ii = 2
!  jj = 1
!  if( (i(ii) == i(jj)).and. &
!      (j(ii) == j(jj)).and. &
!      (k(ii) == k(jj)).and. &
!      (l(ii) == l(jj)) ) then
!     i(ii) = 0
!     exit
!  endif
!! This has been tested with up to 1000 AOs, and all the reverse indices are
!! correct ! We can remove the test
!!    do ii=1,8
!!      if (i(ii) /= 0) then
!!        call two_e_integrals_index(i(ii),j(ii),k(ii),l(ii),i2)
!!        if (i1 /= i2) then
!!          print *,  i1, i2
!!          print *,  i(ii), j(ii), k(ii), l(ii)
!!          stop 'two_e_integrals_index_reverse failed'
!!        endif
!!      endif
!!    enddo
!
!
!end

! ---

