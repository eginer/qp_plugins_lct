program tc_scf
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  my_grid_becke = .True.
  my_n_pt_r_grid = 30
  my_n_pt_a_grid = 50
  touch  my_grid_becke my_n_pt_r_grid my_n_pt_a_grid
!  call save_fock_mos
 call routine_scf
!  call routine_mo
end

subroutine routine_scf
 implicit none
 integer :: i
 i = 0
  print*,'iteration = ',i
  print*,'grad_good_hermit_tc_fock_mat = ',grad_good_hermit_tc_fock_mat
  print*,'***'
   print*,'TC HF total energy = ',TC_right_HF_energy
   print*,'TC HF 1 e   energy = ',TC_right_HF_one_electron_energy
   print*,'TC HF 2 e hermit   = ',TC_right_HF_two_e_hermit_energy
   print*,'TC HF 2 non hermit = ',TC_right_HF_two_e_n_hermit_energy
   print*,'TC HF 3 body       = ',diag_three_elem_hf
  print*,'***'
 do while(grad_good_hermit_tc_fock_mat.gt.thresh_scf)
  i += 1
! do i = 1, 10
  print*,'iteration = ',i
  print*,'grad_good_hermit_tc_fock_mat = ',grad_good_hermit_tc_fock_mat
  print*,'***'
   print*,'TC HF total energy = ',TC_right_HF_energy
   print*,'TC HF 1 e   energy = ',TC_right_HF_one_electron_energy
   print*,'TC HF 2 e hermit   = ',TC_right_HF_two_e_hermit_energy
   print*,'TC HF 2 non hermit = ',TC_right_HF_two_e_n_hermit_energy
   print*,'TC HF 3 body       = ',diag_three_elem_hf
  print*,'***'
  call save_good_hermit_tc_eigvectors
  touch mo_coef 
  call save_mos
 enddo
end

subroutine routine
 implicit none
 integer :: i,j
 double precision :: accu_alpha,accu_beta
 accu_alpha = 0.D0
 accu_beta = 0.D0
 do i = 1, ao_num
  do j = 1, ao_num
   accu_alpha += dabs(two_e_tc_hermit_integral_alpha(j,i) - ao_two_e_integral_alpha(j,i))
   accu_beta  += dabs(two_e_tc_hermit_integral_beta(j,i) - ao_two_e_integral_beta(j,i))
  enddo
 enddo
 print*,'accu_beta  = ',accu_beta
 print*,'accu_alpha = ',accu_alpha

end

subroutine routine_mo
 implicit none
 integer :: i,a,i_ok
 double precision :: f_tc
 double precision :: accu_alpha,accu_beta,hmono,heff,hderiv,hthree,htot
 integer(bit_kind), allocatable :: det_i(:,:)
 print*,'e_tilde_00 = ',e_tilde_00
 print*,'grad_good_hermit_tc_fock_mat = ',grad_good_hermit_tc_fock_mat
 allocate(det_i(N_int,2))
 do i = 1, elec_alpha_num
  do a = elec_alpha_num+1, mo_num
   det_i(:,1) = ref_bitmask(:,1)
   det_i(:,2) = ref_bitmask(:,2)
   call do_single_excitation(det_i,i,a,1,i_ok)
!   f_tc = Fock_matrix_tc_mo_alpha(i,a) ! <HF|H a^dagger_a a_i |HF > = F(i,a)
!   call htilde_mu_mat(ref_bitmask,det_i,hmono,heff,hderiv,hthree,htot)
   f_tc = Fock_matrix_tc_mo_alpha(a,i) ! <HF|H a^dagger_a a_i |HF > = F(i,a)
   call htilde_mu_mat(det_i,ref_bitmask,hmono,heff,hderiv,hthree,htot)
   if(dabs(htot).gt.1.d-10)then
    print*,'i,a',i,a
    print*,'ref,new,dabs'
    print*,htot,f_tc, dabs(f_tc - htot)
   endif
   accu_alpha += (f_tc - htot)
  enddo
 enddo
 print*,'accu_alpha = ',accu_alpha
 print*,''
 print*,'MO Fock matrix '
 do i = 1, mo_num
  write(*,'(100(F10.5,X))')good_hermit_tc_fock_mat(i,:)
 enddo


end

subroutine save_fock_mos
 implicit none
 character*(64) :: label
 integer :: sign,i,j,k
 logical       :: output
 output = .True.
 label = "Canonical"
 sign = 1
! call mo_as_eigvectors_of_mo_matrix(Fock_matrix_tc_mo_tot,mo_num,mo_num,label,sign,output)
 double precision, allocatable :: reigvec_tc_tmp(:,:),leigvec_tc_tmp(:,:),eigval_right_tmp(:)
 allocate(reigvec_tc_tmp(mo_num,mo_num),leigvec_tc_tmp(mo_num,mo_num),eigval_right_tmp(mo_num))
 integer :: n_real_tc_eigval_right

! do j = 1, mo_num
!  do i = 1, j-1
!   Fock_matrix_tc_mo_tot(i,j) = Fock_matrix_tc_mo_tot(j,i) 
!  enddo
! enddo
 print*,'MO Fock matrix '
 do i = 1, mo_num
  write(*,'(100(F10.5,X))')good_hermit_tc_fock_mat(i,:)
 enddo
 call non_hrmt_real_diag(mo_num,Fock_matrix_tc_mo_tot,reigvec_tc_tmp,leigvec_tc_tmp,m,eigval_right_tmp)

! do i = 1, elec_alpha_num
!  print*,'eigenvalues',eigval_right_tmp(i)
!  print*,'right eigenvectors '
!  write(*,'(100(F10.5,X))')reigvec_tc_tmp(:,i)
!  print*,'left  eigenvectors '
!  write(*,'(100(F10.5,X))')leigvec_tc_tmp(:,i)
! enddo
 double precision :: accu, s_mat(mo_num, mo_num)
 s_mat = 0.d0
 do i = 1, mo_num
  do k = 1, mo_num
   do j = 1, mo_num
    s_mat(k,i) += reigvec_tc_tmp(j,i) * reigvec_tc_tmp(j,k)
   enddo
  enddo
 enddo
 print*,'Overlap '
 do i = 1, elec_alpha_num
  write(*,'(100(F16.12,x))')s_mat(i,:)
 enddo
 double precision, allocatable  :: mo_coef_old(:,:)
 allocate(mo_coef_old(ao_num, mo_num))
 mo_coef_old = mo_coef
 integer :: m
 m = mo_num
 call dgemm('N','N',ao_num,m,m,1.d0,mo_coef_old,size(mo_coef_old,1),reigvec_tc_tmp,size(reigvec_tc_tmp,1),0.d0,mo_coef,size(mo_coef,1))
! call dgemm('N','N',ao_num,m,m,1.d0,mo_coef_old,size(mo_coef_old,1),R,size(R,1),0.d0,mo_coef,size(mo_coef,1))
 

end
