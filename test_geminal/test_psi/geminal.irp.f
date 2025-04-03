
subroutine geminal_at_r1_r2(r1,r2,geminal)
 implicit none
 BEGIN_DOC
 ! f(r1,r2) = \sum_{i,j} 1/(phi_(r1)phi_j(r1))\sum_{a,b}t_ij^ab phi_a(r1)phi_b(r2)
 ! 
 ! where i/j is in the set of alpha/beta OCCUPIED orbital,  
 ! 
 ! and t_ij^{ab} is the MP1 wave function
 END_DOC
 double precision, intent(in) :: r1(3), r2(3)
 double precision, intent(out) :: geminal
 double precision, allocatable :: mos_array_r2(:),mos_array_r1(:)
 double precision :: mo_i, mo_j, mo_i_occ, mo_j_occ, inv_hf
 integer :: ii,jj,i_occ,j_occ 
 integer :: a,b,k,aa,bb
 allocate(mos_array_r1(mo_num), mos_array_r2(mo_num))
 call give_all_mos_at_r(r1,mos_array_r1)
 call give_all_mos_at_r(r2,mos_array_r2)
 double precision :: geminal_ij
 geminal = 0.d0
 do i_occ = 1, n_occ_orb_cc
  do j_occ = 1, n_occ_orb_cc
   ii = list_occ_orb_cc(i_occ) ! index in the total MOs set 
   jj = list_occ_orb_cc(j_occ) ! index in the total MOs set 
   mo_i_occ = mos_array_r1(ii) ! values of MO_i(r1)
   mo_i_occ = 1.d0/mo_i_occ    ! and its inverse
   mo_j_occ = mos_array_r2(jj) ! values of MO_j(r2)
   mo_j_occ = 1.d0/mo_j_occ    ! and its inverse
   geminal_ij=0.D0
   do a = 1, n_virt_orb_cc
    aa = list_virt_orb_cc(a)
    do b = 1, n_virt_orb_cc
     bb  = list_virt_orb_cc(b)
     geminal_ij +=  amplitudes_prov(b,a,j_occ,i_occ) * mos_array_r1(aa) * mos_array_r2(bb)
    enddo
   enddo
   geminal += geminal_ij*mo_i_occ* mo_j_occ 
  enddo
 enddo
end

subroutine geminal_ij_at_r1_r2(r1,r2,i_occ,j_occ,geminal)
 implicit none
 BEGIN_DOC
 ! f_ij(r1,r2) = 1/(phi_(r1)phi_j(r1))\sum_{a,b}t_ij^ab phi_a(r1)phi_b(r2)
 ! 
 ! where i_occ/j_occ is in the set of alpha/beta OCCUPIED orbital,  
 ! 
 ! and t_ij^{ab} is the MP1 wave function
 END_DOC
 integer, intent(in) :: i_occ,j_occ ! alpha/beta OCCUPIED orbital pair in the list_virt_orb_cc
 double precision, intent(in) :: r1(3), r2(3)
 double precision, intent(out) :: geminal
 double precision, allocatable :: mos_array_r2(:),mos_array_r1(:)
 double precision :: mo_i, mo_j, mo_i_occ, mo_j_occ, inv_hf
 integer :: ii,jj
 integer :: a,b,k,aa,bb
 allocate(mos_array_r1(mo_num), mos_array_r2(mo_num))
 call give_all_mos_at_r(r1,mos_array_r1)
 call give_all_mos_at_r(r2,mos_array_r2)
 ii = list_occ_orb_cc(i_occ) ! index in the total MOs set 
 jj = list_occ_orb_cc(j_occ) ! index in the total MOs set 
 mo_i_occ = mos_array_r1(ii) ! values of MO_i(r1)
 mo_i_occ = 1.d0/mo_i_occ    ! and its inverse
 mo_j_occ = mos_array_r2(jj) ! values of MO_j(r2)
 mo_j_occ = 1.d0/mo_j_occ    ! and its inverse
 geminal=0.D0
 do a = 1, n_virt_orb_cc
  aa = list_virt_orb_cc(a)
  do b = 1, n_virt_orb_cc
   bb  = list_virt_orb_cc(b)
   geminal +=  amplitudes_prov(b,a,j_occ,i_occ) * mos_array_r1(aa) * mos_array_r2(bb)
  enddo
 enddo
 geminal *= mo_i_occ* mo_j_occ 
end

BEGIN_PROVIDER [integer, n_virt_orb_cc]
 implicit none
 n_virt_orb_cc = mo_num - elec_alpha_num
END_PROVIDER 

BEGIN_PROVIDER [ integer, list_virt_orb_cc, (n_virt_orb_cc)]
 implicit none
 integer :: i,j
 j=1
 do i = elec_alpha_num+1,mo_num
  list_virt_orb_cc(j) = i
  j+=1
 enddo
END_PROVIDER 

BEGIN_PROVIDER [ integer, n_occ_orb_cc]
 implicit none
 n_occ_orb_cc = elec_alpha_num - n_core_orb
END_PROVIDER 

BEGIN_PROVIDER [ integer, list_occ_orb_cc, (n_occ_orb_cc)]
 implicit none
 integer :: i,j
 j = 1
 do i = n_core_orb+1,elec_alpha_num
  list_occ_orb_cc(j) = i
  j+=1
 enddo
END_PROVIDER 

BEGIN_PROVIDER [ double precision, amplitudes_prov, (n_virt_orb_cc, n_virt_orb_cc,n_occ_orb_cc, n_occ_orb_cc)]
 implicit none
 integer :: i,j,a,b,ii,jj,aa,bb
 double precision :: hij,eps_i,eps_j,eps_a,eps_b,amplitude, get_two_e_integral
 do i = 1, n_occ_orb_cc
  ii = list_occ_orb_cc(i)
  eps_i = fock_matrix_diag_mo(ii)
  do j = 1, n_occ_orb_cc
   jj = list_occ_orb_cc(j)
   eps_j = fock_matrix_diag_mo(jj)
   do a = 1, n_virt_orb_cc
    aa = list_virt_orb_cc(a)
    eps_a = fock_matrix_diag_mo(aa)
    do b = 1, n_virt_orb_cc
     bb  = list_virt_orb_cc(b)
     eps_b = fock_matrix_diag_mo(bb)
     hij = get_two_e_integral(ii,jj,aa,bb,mo_integrals_map)
     amplitude = hij/(eps_i+eps_j-eps_a-eps_b)
     amplitudes_prov(b,a,j,i) = amplitude
    enddo
   enddo
  enddo
 enddo
END_PROVIDER 
