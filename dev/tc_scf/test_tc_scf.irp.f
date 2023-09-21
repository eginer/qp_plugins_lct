program tc_scf

  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC

  implicit none

  print *, 'starting ...'

  my_grid_becke  = .True.
  my_n_pt_r_grid = 30
  my_n_pt_a_grid = 50
  touch my_grid_becke my_n_pt_r_grid my_n_pt_a_grid

  bi_ortho = .True.
  touch bi_ortho

  !call print_fock
  !call test_completeness()
  call test_overlap()

end

subroutine print_fock
 implicit none
 integer :: i,j
 print*,'Diagonal of Fock matrix'
 do i = 1, mo_num
  print*,i,Fock_matrix_tc_mo_tot(i,i)
 enddo
 print*,'Fock matrix '
 do i = 1, mo_num
  write(*,'(100(F16.10,X))')Fock_matrix_tc_mo_tot(:,i)
 enddo
end

! ---

subroutine test_completeness()

  implicit none
  integer                       :: i, j, k
  double precision              :: diff
  double precision, allocatable :: M1(:,:), M2(:,:), R(:,:), L(:,:), tmp(:,:)

  allocate(R(ao_num,1), L(1,ao_num), tmp(ao_num,ao_num))
  allocate(M1(ao_num,ao_num), M2(ao_num,ao_num))

  M1 = 0.d0
  M2 = 0.d0

  !do i = 1, mo_num
  do i = 3, 3

    !do j = 1, ao_num
    !  R(j,1) = mo_r_coef(j,i)
    !  L(1,j) = mo_l_coef(j,i)
    !enddo
    !call dgemm( 'N', 'N', ao_num, ao_num, 1, 1.d0 &
    !          , R, size(R, 1), L, size(L, 1)      &
    !          , 0.d0, tmp, size(tmp, 1) )
    !do j = 1, ao_num
    !  do k = 1, ao_num
    !    M1(k,j) += tmp(k,j)
    !  enddo
    !enddo

    ! ---

    do j = 1, ao_num
      do k = 1, ao_num
        M1(k,j) += mo_r_coef(k,i) * mo_l_coef(j,i)
      enddo
    enddo

    do j = 1, ao_num
      do k = 1, ao_num
        M2(k,j) += mo_coef(k,i) * mo_coef(j,i)
      enddo
    enddo

  enddo

  print*, " comleteness HTC "
  do i = 1, ao_num
    write(*, '(100(F16.10,X))') M1(i,:)
  enddo

  print*, " comleteness H "
  do i = 1, ao_num
    write(*, '(100(F16.10,X))') M2(i,:)
  enddo

  diff = 0.d0
  do j = 1, ao_num
    do k = 1, ao_num
      diff += dabs(M1(k,j) - M2(k,j))
    enddo
  enddo
  print*, ' diff = ', diff

  deallocate(R, L, tmp)
  deallocate(M1)
  deallocate(M2)

end subroutine test_completeness

! ---

subroutine test_overlap()

  implicit none
  integer                       :: i
  double precision, allocatable :: C(:,:)
  double precision, allocatable :: L(:,:), R(:,:), S(:,:), M(:,:)

  allocate(L(ao_num,mo_num), R(ao_num,mo_num), S(ao_num,ao_num), M(ao_num,mo_num))
  L = mo_l_coef
  R = mo_r_coef
  M = mo_coef
  S = ao_overlap

  allocate(C(mo_num,mo_num))

  print*, " C.T x C"
  call LTxSxR(ao_num, mo_num, M, S, M, C)
  do i = 1, mo_num
    write(*, '(100(F16.10,X))') C(i,:)
  enddo

  print*, " L.T x R"
  call LTxSxR(ao_num, mo_num, L, S, R, C)
  do i = 1, mo_num
    write(*, '(100(F16.10,X))') C(i,:)
  enddo

  print*, " L.T x C"
  call LTxSxR(ao_num, mo_num, L, S, M, C)
  do i = 1, mo_num
    write(*, '(100(F16.10,X))') C(i,i)
  enddo

  print*, " C.T x R"
  call LTxSxR(ao_num, mo_num, M, S, R, C)
  do i = 1, mo_num
    write(*, '(100(F16.10,X))') C(i,i)
  enddo

  print*, " L.T x L"
  call LTxSxR(ao_num, mo_num, L, S, L, C)
  do i = 1, mo_num
    write(*, '(100(F16.10,X))') C(i,i)
  enddo

  print*, " R.T x R"
  call LTxSxR(ao_num, mo_num, R, S, R, C)
  do i = 1, mo_num
    write(*, '(100(F16.10,X))') C(i,i)
  enddo

  deallocate(C)
  deallocate(L, R, S, M)

end subroutine test_overlap

! ---
