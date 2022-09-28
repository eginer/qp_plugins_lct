subroutine i_soc_j(conf_i,conf_j,det_i, det_j,sz_i,sz_j,soc_ij)
 implicit none
 BEGIN_DOC
! i_soc_j = <det_i, sz_i| SOC | det_j, sz_j>
 END_DOC
 use bitmasks ! you need to include the bitmasks_module.f90 features
 integer(bit_kind), intent(in) :: det_i(N_int,2), det_j(N_int,2)
 integer(bit_kind), intent(in) :: conf_i(N_int,2), conf_j(N_int,2)
 double precision, intent(in) :: sz_i, sz_j
 complex*8, intent(out) :: soc_ij
 BEGIN_DOC
! det_i and det_j are two determinants of possibly different ms
!
! i_soc_j = <det_i | SOC | det_j>
!
! |det_i> is also written with conf_i:: conf_i(1,:) = doubly occupied, conf_i(2,:) = somo
 END_DOC
 integer :: i,j, degree_doc, degree_somo,n_doc_i,n_doc_j,n_somo_i,n_somo_j,h,p
 integer :: degree,i_ok
 integer(bit_kind) :: somo_i(N_int,2),somo_j(N_int,2),somo_dif_ab(N_int,2),conf_dif(N_int,2)
 integer(bit_kind) :: det_tmp(N_int,2),det_tmp_2(N_int,2),domo_dif(N_int,2)
 integer :: somo_dif_ab_occ(N_int*bit_kind_size,2), n_somo_dif_ab(2)
 double precision :: phase,phase_2
 integer          :: exc(0:2,2,2),pm,spin,orb_j

 soc_ij = complex(0.d0,0.d0)
 if(dabs(sz_i-sz_j).gt.1.d0)return
 call get_good_things_dets(conf_i,conf_j,det_i, det_j,n_doc_i,n_doc_j,degree_doc,n_somo_i,n_somo_j,degree_somo)
 if(abs(n_doc_j - n_doc_i).gt.1)return
 if(abs(n_somo_j - n_somo_i).gt.2)return
 call bitstring_to_list_ab(somo_dif_ab, somo_dif_ab_occ, n_somo_dif_ab, N_int)
 call get_diff_dets_somo_doc(det_i,det_j,conf_i,conf_j, somo_i, somo_j, somo_dif_ab, conf_dif)
 if(dabs(sz_i- sz_j).gt.1.d-10)then    
  ! S^+ V^- or S^-V^+
  if(degree_doc == 0 .and. degree_somo.lt.2)then ! doubly occupied untouched 
   if(sz_i .gt. sz_j)then  ! <det_i|S^+ V^-|det_j> == <det_i| a^{dagger}_p_alpha a_q_beta | det_j>
    p = somo_dif_ab_occ(1,1) ! particle  a^{dagger}_p_alpha | det_i> 
    h = somo_dif_ab_occ(1,2) ! hole      a_q_beta           | det_i> 
    spin = 2
    pm = +1
    soc_ij = mo_v_soc_tot(p,h,1) ! V^-
   else                    ! <det_i|S^- V^+|det_j> == <det_i| a^{dagger}_p_beta a_q_alpha | det_j>
    p = somo_dif_ab_occ(1,2) ! particle  a^{dagger}_p_beta | det_i> 
    h = somo_dif_ab_occ(1,1) ! hole      a_q_alpha         | det_i> 
    spin = 1
    pm = -1
    soc_ij = mo_v_soc_tot(p,h,2) ! V^+
   endif
   if(h==p)then ! spin flip without single excitation 
    orb_j = h
    call s_plus_or_s_minus_det_orb_j(det_j,orb_j,pm,det_tmp,phase) ! S^+|det_i> == phase * |det_j>
    call get_excitation_degree(det_i,det_tmp,degree,N_int)
    if(degree .ne. 0 )then
     print*,'Probleeeemmm '
     stop
    endif
   else ! firt do excitation from h_beta --> p_beta and then spin-flip
    det_tmp = det_j
    call do_single_excitation(det_tmp,h,p,spin,i_ok) ! case S^+ :: |det_tmp> = a^\dagger_p_beta a_h_beta | det_j>
    call get_single_excitation(det_j,det_tmp,exc,phase,N_int) ! phase associated with singl exc
    call s_plus_or_s_minus_det_orb_j(det_tmp,p,det_tmp_2,phase_2) ! case S^+ :: S^+|det_j> == phase_2 * phase_1 |det_j>
    phase *= phase_2
    if(i_ok .ne. 1)then
     print*,'PROBLEMMMM'
     stop
    endif
    call get_excitation_degree(det_i,det_tmp_2,degree,N_int)
    if(degree .ne. 0 )then
     print*,'Probleeeemmm bis '
     stop
    endif
   endif
   soc_ij *= complex(phase,0.d0)
  else ! doubly occupied have changed 
   do i = 1, N_int
    det_tmp(i,1) = conf_dif(i,1) ! new DOMO created/destroyed 
    det_tmp(i,2) = xor(conf_dif(i,1),conf_dif(i,2)) ! only the new SOMO created/destroyed
   enddo
   call bitstring_to_list_ab(det_tmp, somo_dif_ab_occ, n_somo_dif_ab, N_int)
   if(n_doc_i .gt. n_doc_j)then ! particle is on domo, hole is on somo
    p = somo_dif_ab_occ(1,1)
    h = somo_dif_ab_occ(1,2) 
   else                         ! particle is on somo, hole is on domo
    p = somo_dif_ab_occ(1,2)
    h = somo_dif_ab_occ(1,1) 
   endif
   if(sz_i .gt. sz_j)then  ! <det_i|S^+ V^-|det_j> == <det_i| a^{dagger}_p_alpha a_q_beta | det_j>
    pm = +1
    spin = +1 
    soc_ij = mo_v_soc_tot(p,h,1) ! V^-
   else
    pm = -1
    spin = -1
    soc_ij = mo_v_soc_tot(p,h,2) ! V^+
   endif
   det_tmp = det_i
   call do_single_excitation(det_tmp,p,h,spin,i_ok) ! case S^- :: |det_tmp> = a^\dagger_h_beta a_p_beta | det_i>
   if(i_ok .ne. 1)then
    print*,'something went wrong in do_single_excitation'
   endif
   call get_single_excitation(det_i,det_tmp,exc,phase,N_int) ! phase associated with singl exc
   call s_plus_or_s_minus_det_orb_j(det_j,pm,det_tmp_2,phase_2) ! case S^- :: S^-|det_j> == phase_2 * phase_1 |det_j>
   phase *= phase_2
   call get_excitation_degree(det_tmp,det_tmp_2,degree,N_int)
   if(degree .ne. 0 )then
    print*,'Probleeeemmm tris '
    stop
   endif
   soc_ij *= complex(phase,0.d0)
  endif
 else if(sz_j == sz_i)then ! <det_i|S^z V^z|det_j>
  call i_sz_j(det_i,det_j,sz_i,sz_j,soc_ij)
 endif
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
 double precision , intent(in) :: sz_i, sz_j
 complex*8, intent(out) :: soc_ij

 integer, allocatable           :: occ(:,:)
 integer                        :: n_occ_ab(2),iorb
 integer          :: exc(0:2,2,2),h,p,i
 integer          :: degree,ispin
 double precision :: phase,sign_spin(2)

 soc_ij = complex(0.d0,0.d0)
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
    soc_ij += mo_v_soc_tot(iorb,iorb,3) * complex(sign_spin(ispin),0.d0) ! contribution from orbital iorb with the good sign
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
  soc_ij += mo_v_soc_tot(p,h,3) * complex(sign_spin(ispin) * phase,0.d0) ! contribution of a^{dagger}_{p,ispin} a_{h,ispin} 
 endif

end
