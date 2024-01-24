
BEGIN_PROVIDER [double precision, mo_j1b_gauss_hermI, (mo_num,mo_num)]

  BEGIN_DOC
  !
  ! Hermitian part of 1-body Jastrow factow in the |MO| basis set.
  !
  !  :math:`\langle \chi_A | -0.5 \Delta \tau_{1b} | \chi_B \rangle` 
  !
  END_DOC

  implicit none

  mo_j1b_gauss_hermI = 0.d0
  call ao_to_mo(    j1b_gauss_hermI, size(   j1b_gauss_hermI, 1) &
               , mo_j1b_gauss_hermI, size(mo_j1b_gauss_hermI, 1) )

END_PROVIDER



BEGIN_PROVIDER [double precision, mo_j1b_gauss_hermII, (mo_num,mo_num)]

  BEGIN_DOC
  !
  ! Hermitian part of 1-body Jastrow factow in the |MO| basis set.
  !
  !  :math:`\langle \chi_A | -0.5 \grad \tau_{1b} \cdot \grad \tau_{1b} | \chi_B \rangle` 
  !
  END_DOC

  implicit none

  mo_j1b_gauss_hermII = 0.d0
  call ao_to_mo(    j1b_gauss_hermII, size(   j1b_gauss_hermII, 1) &
               , mo_j1b_gauss_hermII, size(mo_j1b_gauss_hermII, 1) )

END_PROVIDER



BEGIN_PROVIDER [ double precision, mo_j1b_gauss_nonherm, (mo_num,mo_num)]

  BEGIN_DOC
  !
  ! Hermitian part of 1-body Jastrow factow in the |MO| basis set.
  !
  !              \langle \chi_i | - grad \tau_{1b} \cdot grad | \chi_j \rangle  = 
  !  2 \sum_A aA \langle \chi_i | exp[-aA riA^2] (ri-rA) \cdot grad | \chi_j \rangle
  !
  END_DOC

  implicit none
  integer                       :: i, j
  double precision, allocatable :: j1b_t(:,:)

  allocate(j1b_t(ao_num,ao_num))
  do i = 1, ao_num
    do j = 1, ao_num
      !j1b_t(j,i) = j1b_gauss_nonherm(i,j)
      j1b_t(j,i) = j1b_gauss_nonherm(j,i)
    enddo
  enddo

  mo_j1b_gauss_nonherm = 0.d0
  !call ao_to_mo(    j1b_gauss_nonherm, size(   j1b_gauss_nonherm, 1) &
  call ao_to_mo( j1b_t, size(j1b_t, 1) &
               , mo_j1b_gauss_nonherm, size(mo_j1b_gauss_nonherm, 1) )
  
  deallocate( j1b_t )
  
END_PROVIDER
