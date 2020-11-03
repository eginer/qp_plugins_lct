program write_rot_mat
 implicit none
 read_wf = .True.
 touch read_wf
 call write_rotation_matrix
end

subroutine write_rotation_matrix
 implicit none
 integer :: i,j,k,l
 integer :: iorb,jorb
 double precision :: eigvalues(mo_num),rot_mat_tmp(mo_num,mo_num), mo_array(mo_num, mo_num)
 double precision :: rot_mat_act(n_act_orb, n_act_orb)

 mo_array = -one_e_dm_mo
 print*,'Write rotation matrix from current orbitals to natural orbitals '
 call lapack_diagd(eigvalues,rot_mat_tmp,mo_array,mo_num,mo_num) 
 do i = 1, n_act_orb
  iorb = list_act(i)
  do j = 1, n_act_orb
   jorb = list_act(j)
   rot_mat_act(j,i) = rot_mat_tmp(jorb,iorb)
  enddo
 enddo
 open(1, file = 'rotation_matrix') 
 do i = 1, mo_num
  write(1,*)rot_mat_tmp(i,1:mo_num)
 enddo
 close(1) 
 double precision :: new_mo_coef_tmp(ao_num, mo_num)
 ! <AO_k| new_mo_j> = \sum_i U_ij <AO_k| old_mo_i>
 new_mo_coef_tmp = 0.d0
 do j = 1, mo_num ! 
  do i = 1, mo_num
   do k = 1, ao_num
     new_mo_coef_tmp(k,j) += mo_coef(k,i) * rot_mat_tmp(i,j) 
   enddo
  enddo
 enddo
 open(1, file = 'new_mo_coef') 
 do i = 1, mo_num
  write(1,'(100(F16.10,X))')new_mo_coef_tmp(i,1:mo_num)
 enddo
 close(1) 
 double precision :: mo_overlap_tmp(mo_num, mo_num)
 double precision :: one_e_rdm_alpha(n_act_orb, n_act_orb), one_e_rdm_beta(n_act_orb, n_act_orb)
 print*,'Computing the new overlap '
 call new_one_e_mat(ao_overlap, new_mo_coef_tmp, mo_overlap_tmp)
 open(1, file = 'mo_overlap_guess') 
 do i = 1, mo_num
  write(1,'(100(F16.13,X))')mo_overlap_tmp(i,:)
 enddo
 close(1) 
 print*,'computing the new one-rdm alpha'
! call new_one_e_mat_act(one_e_dm_mo_alpha_average, rot_mat_act, one_e_rdm_alpha)
 call new_one_e_mat_act(one_e_dm_mo_alpha, rot_mat_act, one_e_rdm_alpha)
 open(1, file = 'one_rdm_a') 
 do i = 1, n_act_orb
  write(1,'(100(F16.10,X))')one_e_rdm_alpha(i,1:n_act_orb)
 enddo
 close(1) 

 print*,'computing the new one-rdm beta'
 call new_one_e_mat_act(one_e_dm_mo_beta, rot_mat_act, one_e_rdm_beta)
 open(1, file = 'one_rdm_b') 
 do i = 1, n_act_orb
  write(1,'(100(F16.10,X))')one_e_rdm_beta(i,1:n_act_orb)
 enddo
 close(1) 

 call routine_active_only_test(act_2_rdm_ab_mo)
 double precision, allocatable :: two_rdm(:, :, :, :), two_e_ints(:,:,:,:)
 allocate( two_rdm(n_act_orb, n_act_orb, n_act_orb, n_act_orb) )
 allocate( two_e_ints(n_act_orb, n_act_orb, n_act_orb, n_act_orb) )

 ! transforming the two RDM with the new MOs
 call new_two_e_mat(rot_mat_act,act_2_rdm_ab_mo,two_rdm)
 ! transforming the two e integrals with the new MOs
 call new_two_e_mat(rot_mat_act,vee_big_array,two_e_ints)
 ! testing the alpha-beta two e energy
 call routine_active_only_test_bis(two_rdm, two_e_ints)
 ! writing in plain text the two RDM
 call write_two_rdm(two_rdm)
 deallocate(two_rdm)
end


subroutine new_one_e_mat(one_e_mat_ao, new_mos, one_e_mat_mo)
 implicit none
 double precision, intent(in) :: one_e_mat_ao(ao_num,ao_num), new_mos(ao_num, mo_num)
 double precision, intent(out):: one_e_mat_mo(mo_num, mo_num)
 integer :: i,j,k,l
 one_e_mat_mo = 0.d0
 do i = 1, mo_num
  do j = 1, mo_num
   do k = 1, ao_num
    do l = 1, ao_num
     one_e_mat_mo(j,i) += new_mos(l,j) * new_mos(k,i) * one_e_mat_ao(l,k)
    enddo
   enddo
  enddo
 enddo
end

subroutine new_one_e_mat_bis(one_e_mat_mo, rot_mat, one_e_mat_mo_2)
 implicit none
 double precision, intent(in) :: one_e_mat_mo(mo_num,mo_num), rot_mat(mo_num, mo_num)
 double precision, intent(out):: one_e_mat_mo_2(mo_num, mo_num)
 integer :: i,j,k,l
 one_e_mat_mo_2 = 0.d0
 do i = 1, mo_num
  do j = 1, mo_num
   do k = 1, mo_num
    do l = 1, mo_num
     one_e_mat_mo_2(j,i) += rot_mat(l,j) * rot_mat(k,i) * one_e_mat_mo(l,k)
    enddo
   enddo
  enddo
 enddo
end

subroutine new_one_e_mat_act(one_e_mat_mo, rot_mat, one_e_mat_mo_2)
 implicit none
 double precision, intent(in) :: one_e_mat_mo(mo_num,mo_num), rot_mat(n_act_orb, n_act_orb)
 double precision, intent(out):: one_e_mat_mo_2(n_act_orb, n_act_orb)
 integer :: i,j,k,l
 double precision :: mo_act_array(n_act_orb,n_act_orb)

 do i = 1, n_act_orb
  j = list_act(i)
  do k = 1, n_act_orb
   l = list_act(k)
   mo_act_array(i,k) = one_e_mat_mo(j,l)
  enddo
 enddo
! do i =1, n_act_orb
!  print*,'',mo_act_array(i,:)
! enddo
 one_e_mat_mo_2 = 0.d0
 do i = 1, n_act_orb
  do j = 1, n_act_orb
   do k = 1, n_act_orb
    do l = 1, n_act_orb
     one_e_mat_mo_2(j,i) += rot_mat(l,j) * rot_mat(k,i) * mo_act_array(l,k)
    enddo
   enddo
  enddo
 enddo
end


subroutine new_two_e_mat(rot_mat,two_e_array_in,two_e_array)
 implicit none
 double precision, intent(in) :: rot_mat(n_act_orb, n_act_orb),two_e_array_in(n_act_orb, n_act_orb, n_act_orb, n_act_orb)
 double precision, intent(out) :: two_e_array(n_act_orb, n_act_orb, n_act_orb, n_act_orb)
 double precision, allocatable :: two_e_array_tmp_1(:,:,:,:),two_e_array_tmp_2(:,:,:,:)
 integer :: i,j,k,l
 integer :: m,n,p,q
 allocate( two_e_array_tmp_1(n_act_orb, n_act_orb, n_act_orb, n_act_orb) )
 allocate( two_e_array_tmp_2(n_act_orb, n_act_orb, n_act_orb, n_act_orb) )
 two_e_array_tmp_1 = 0.d0
 do l = 1, n_act_orb
  do k = 1, n_act_orb
   do j = 1, n_act_orb
    do m = 1, n_act_orb
     do i = 1, n_act_orb
      ! <mj|kl> = \sum_i <i|m> * <ij|kl>
!     two_e_array_tmp_1(m,j,k,l) += rot_mat(i,m) * act_2_rdm_ab_mo(i,j,k,l,1)
      two_e_array_tmp_1(m,j,k,l) += rot_mat(i,m) * two_e_array_in(i,j,k,l)
     enddo
    enddo
   enddo
  enddo
 enddo
 two_e_array_tmp_2 = 0.d0
 do l = 1, n_act_orb
  do k = 1, n_act_orb
   do n = 1, n_act_orb
    do j = 1, n_act_orb
     do m = 1, n_act_orb
     ! <mn|kl> = \sum_j <j|n> <mj|kl>
      two_e_array_tmp_2(m,n,k,l) += rot_mat(j,n) * two_e_array_tmp_1(m,j,k,l)
     enddo
    enddo
   enddo
  enddo
 enddo

 two_e_array_tmp_1 = 0.d0
 do l = 1, n_act_orb
  do p = 1, n_act_orb
   do k = 1, n_act_orb
    do n = 1, n_act_orb
     do m = 1, n_act_orb
      ! <mn|pl> = \sum_k <k|p> <mn|kl>
      two_e_array_tmp_1(m,n,p,l) += rot_mat(k,p) * two_e_array_tmp_2(m,n,k,l) 
     enddo
    enddo
   enddo
  enddo
 enddo
 deallocate(two_e_array_tmp_2)

 two_e_array = 0.d0
 do q = 1, n_act_orb
  do l = 1, n_act_orb
   do p = 1, n_act_orb
    do n = 1, n_act_orb
     do m = 1, n_act_orb
      ! <mn|pq> = \sum_l <l|q> <mn|pl>
      two_e_array(m,n,p,q) += rot_mat(l,q) * two_e_array_tmp_1(m,n,p,l) 
     enddo
    enddo
   enddo
  enddo
 enddo
end

BEGIN_PROVIDER [double precision, vee_big_array, (n_act_orb, n_act_orb, n_act_orb, n_act_orb)]
 implicit none
 integer :: i,j,k,l,iorb,jorb,korb,lorb,istate
 double precision :: vijkl,get_two_e_integral

   do i = 1, n_act_orb
    iorb = list_act(i)
    do j = 1, n_act_orb
     jorb = list_act(j)
     do k = 1, n_act_orb
      korb = list_act(k)
      do l = 1, n_act_orb
       lorb = list_act(l)

       vijkl = get_two_e_integral(lorb,korb,jorb,iorb,mo_integrals_map)                                 
       vee_big_array(l,k,j,i) = vijkl
      enddo
     enddo
    enddo
   enddo
END_PROVIDER 

subroutine routine_active_only_test_bis(two_rdm, two_e_ints)
 implicit none
 double precision, intent(in) :: two_rdm(n_act_orb, n_act_orb, n_act_orb, n_act_orb)
 double precision, intent(in) :: two_e_ints(n_act_orb, n_act_orb, n_act_orb, n_act_orb)
 BEGIN_DOC
! This routine computes the two electron repulsion within the active space using various providers 
! 
 END_DOC
 double precision :: vijkl
 double precision :: wee_ab,rdmab
 integer :: iorb,jorb,korb,lorb
 wee_ab  = 0.d0
 print*,'**************************'
 print*,'**************************'
   !! PURE ACTIVE PART 
   !! 
   do iorb = 1, n_act_orb
    do jorb = 1, n_act_orb
     do korb = 1, n_act_orb
      do lorb = 1, n_act_orb
       vijkl = two_e_ints(lorb,korb,jorb,iorb)
       rdmab =  two_rdm(lorb,korb,jorb,iorb)
   !   print*,iorb,jorb,korb,lorb
   !   print*,rdmab
   !   print*,vijkl
       wee_ab             += vijkl * rdmab
      enddo
     enddo
    enddo
   enddo
   print*,''
   print*,''
   print*,'wee_ab(istate)          = ',wee_ab
   print*,''
end

 subroutine write_two_rdm(two_rdm)
 implicit none
 integer :: i,j,k,l
 double precision, intent(in) :: two_rdm(n_act_orb, n_act_orb, n_act_orb, n_act_orb)
 double precision :: value_rdm
 character*(1) :: coma
 coma = ","
 open(1, file = 'two_rdm') 
 do i = 1, n_act_orb
  do j = 1, n_act_orb
   do k = 1, n_act_orb
    do l = 1, n_act_orb
     value_rdm = two_rdm(l,k,j,i)
     write(1,'(4(I3,A1),F16.13)')l,coma, k, coma, j, coma, i, coma, value_rdm
    enddo
   enddo
  enddo
 enddo
 close(1)

 end
