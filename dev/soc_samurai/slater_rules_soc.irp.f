subroutine i_soc_j(conf_i,conf_j,det_i, det_j,sz_i,sz_j,soc_ij)
 implicit none
 use bitmasks ! you need to include the bitmasks_module.f90 features
 integer(bit_kind), intent(in) :: det_i(N_int,2), det_j(N_int,2)
 integer(bit_kind), intent(in) :: conf_i(N_int,2), conf_j(N_int,2)
 integer          , intent(in) :: sz_i, sz_j
 double precision, intent(out) :: soc_ij
 BEGIN_DOC
! det_i and det_j are two determinants of possibly different ms
!
! i_soc_j = <det_i | SOC | det_j>
 END_DOC
 integer :: i,j, degree_doc, degree_somo,n_doc_i,n_doc_j,n_somo_i,n_somo_j
 integer :: degree
 integer(bit_kind) :: det_i_tmp(N_int), det_j_tmp(N_int)

 soc_ij = 0.d0
 if(dabs(dble(sz_i-sz_j)).gt.1.d0)return
 ! degree of excitation of doubly occupied part 
 degree_doc = 0
 n_doc_i = 0
 n_doc_j = 0
 do i = 1, N_int
  degree_doc += xor(conf_i(i,1),conf_j(i,1))
  n_doc_i += popcnt(conf_i(i,1))
  n_doc_j += popcnt(conf_j(i,1))
 enddo
 degree_doc = degree_doc / 2
 ! degree of excitation of single occupied part 
 degree_somo = 0
 do i = 1, N_int
  degree_somo += xor(conf_i(i,2),conf_j(i,2))
  n_somo_i += popcnt(conf_i(i,2))
  n_somo_j += popcnt(conf_j(i,2))
 enddo
 degree_somo = degree_somo / 2

 if(degree_somo .gt. 1)then
  return
 endif
 if(degree_doc .gt. 2)then
  return
 endif
 
 if(degree_doc == 0 .and. degree_somo == 1)then
  ! single from singly occupied to empty
 else if(degree_somo == 2 .and. degree_doc == 1)then
  ! single from somo to somo 
 else if(degree_doc == 0 .and. degree_somo == 0)then
  ! either diagonal or spin-flip
  if(sz_j == sz_i)then
   call i_sz_j(det_i,det_j,sz_i,sz_j,soc_ij)
  endif
 endif
! if(degree_doc == 0)then
!  if(degree_somo == 0)then
!   ! same conf 
!  else if(degree_somo == 1)then
!   ! differ by a single exc of somo
!  else 
!   return
!  endif
! else if(degree_doc == 1)then
!  if(n_doc_i == n_doc_j .and. degree_somo == 0)then
!  ! double excitation from a closed shell to a closed shell
!  return
!  else if(n_doc_i == n_doc_j .and. degree_somo == 1)then
!  ! single excitation from a somo to another somo
!  else if(abs(n_doc_i - n_doc_j) == 1 .and. abs(n_somo_i - n_somo_j)==2)then
!  ! single excitation from a doubly occupied to an empty orbital 
!  else if(abs(n_doc_i - n_doc_j) == 0 .and. abs(n_somo_i - n_somo_j)==0 .and. degree_somo == 1)then
!  ! single excitation from a doubly occupied to a somo ==> transfer of doubly occupied orbitals 
!  endif
! endif

end

subroutine i_sz_j(det_i,det_j,sz_i,sz_j,soc_ij)
 implicit none
 BEGIN_DOC
! i_sz_j = <det_i| V^z | det_j> where det_i and det_j can have different S_z values and where 
! 
! V^z = \sum_{p,q} V_pq^z [ a^{dagger}_{p,alpha} a_{q,\alpha}  -  a^{dagger}_{p,beta} a_{q,beta} ]
 END_DOC
 use bitmasks ! you need to include the bitmasks_module.f90 features
 integer(bit_kind), intent(in) :: det_i(N_int,2), det_j(N_int,2)
 integer          , intent(in) :: sz_i, sz_j
 double precision, intent(out) :: soc_ij

 integer, allocatable           :: occ(:,:)
 integer                        :: n_occ_ab(2),iorb
 integer          :: exc(0:2,2,2),h,p,i
 integer          :: degree,ispin
 double precision :: phase,sign_spin(2)

 soc_ij = 0.d0
 if(sz_i.ne.sz_j)return ! S_z does not couple two determinants with different S_z
 sign_spin(1) =+1.d0
 sign_spin(2) =-1.d0
 call get_excitation_degree(det_i,det_j,degree,N_int)

 if (degree == 0)then! diagonal case 
  allocate(occ(N_int*bit_kind_size,2))
  call bitstring_to_list_ab(det_i, occ, n_occ_ab, N_int)        
  do ispin = 1,2 ! loop over the spins of electrons : ispin = 1 :: alpha , ispin = 2 :: beta
   do i = 1, n_occ_ab(ispin) ! loop over electrons of spin "ispin"
    iorb = occ(i,ispin) ! ith ispin-occupied orbital in det_i
    soc_ij += v_soc_tot_mo(iorb,iorb,3) * sign_spin(ispin) ! contribution from orbital iorb with the good sign
   enddo
  enddo
 else if(degree == 1)then
  ! check wether det_i = a^{dagger}_p a_h det_j or the opposite
  call get_single_excitation(det_i,det_j,exc,phase,N_int)
  if (exc(0,1,1) == 1) then ! single alpha : goes with +1
    h = exc(1,1,1)
    p = exc(1,2,1)
    ispin = 1
  else                      ! single beta : goes with -1
    h = exc(1,1,2)
    p = exc(1,2,2)
    ispin = 2
  endif
  soc_ij += v_soc_tot_mo(p,h,3) * sign_spin(ispin) * phase ! contribution of a^{dagger}_{p,ispin} a_{h,ispin} 
 endif

end
