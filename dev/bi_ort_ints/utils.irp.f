
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
  integer :: i,j,p,q

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

