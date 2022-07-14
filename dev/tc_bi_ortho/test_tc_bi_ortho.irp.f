program tc_bi_ortho
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  print *, 'Hello world'
  my_grid_becke = .True.
  my_n_pt_r_grid = 30
  my_n_pt_a_grid = 50
  read_wf = .True.
  touch read_wf
  touch  my_grid_becke my_n_pt_r_grid my_n_pt_a_grid

  call routine_2
end

subroutine routine_2
 implicit none
 integer :: i,j
 double precision :: bi_ortho_mo_ints
 print*,'H matrix'
 do i = 1, N_det
  write(*,'(1000(F16.5,X))')htilde_matrix_elmt_bi_ortho(:,i)
 enddo
 i = 1
 j = 2
 double precision :: phase
 integer :: degree,h1, p1, h2, p2, s1, s2, exc(0:2,2,2)
 call get_excitation_degree(psi_det(1,1,j), psi_det(1,1,i), degree, N_int)
 if(degree==2)then
  call get_double_excitation(psi_det(1,1,j), psi_det(1,1,i), exc, phase, N_int)
  call decode_exc(exc, 2, h1, p1, h2, p2, s1, s2)
  print*,'h1,h2,p1,p2'
  print*, h1,h2,p1,p2 
  print*,mo_bi_ortho_tc_two_e(p1,p2,h1,h2),mo_bi_ortho_tc_two_e(h1,h2,p1,p2)
 endif

 
 print*,'coef'
 do i = 1, ao_num
  print*,i,mo_l_coef(i,8),mo_r_coef(i,8)
 enddo
! print*,'mdlqfmlqgmqglj'
! print*,'mo_bi_ortho_tc_two_e()',mo_bi_ortho_tc_two_e(2,2,3,3)
! print*,'bi_ortho_mo_ints      ',bi_ortho_mo_ints(2,2,3,3)
 print*,'Overlap'
 do i = 1, mo_num
  write(*,'(100(F16.10,X))')overlap_bi_ortho(:,i)
 enddo

end

subroutine routine
 implicit none
 double precision :: hmono,htwoe,hthree,htot
 integer(bit_kind), allocatable  :: key1(:,:)
 integer(bit_kind), allocatable  :: key2(:,:)
 allocate(key1(N_int,2),key2(N_int,2))
 use bitmasks
 key1 = ref_bitmask
 call htilde_mu_mat_bi_ortho(key1,key1, N_int, hmono,htwoe,hthree,htot)
 key2 = key1
 integer :: h,p,i_ok
 h = 1
 p = 8
 call do_single_excitation(key2,h,p,1,i_ok) 
 call debug_det(key2,N_int)
 call htilde_mu_mat_bi_ortho(key2,key1, N_int, hmono,htwoe,hthree,htot)
! print*,'fock_matrix_tc_mo_alpha(p,h) =  ',fock_matrix_tc_mo_alpha(p,h)
 print*,'htot                         =  ',htot
 print*,'hmono                        =  ',hmono 
 print*,'htwoe                        =  ',htwoe
 double precision :: bi_ortho_mo_ints
 print*,'bi_ortho_mo_ints(1,p,1,h)',bi_ortho_mo_ints(1,p,1,h)

end
