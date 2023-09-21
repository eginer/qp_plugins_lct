! For a derivation of the transcorrelated Hamiltonian matrix elements, see Luo-JCP-10 
! http://dx.doi.org/10.1063/1.3505037
subroutine diag_htilde_mu_mat(key_i,hmono,herf,heff,hderiv,hthree,htot)
  use bitmasks
  BEGIN_DOC
!  diagonal element of htilde 
  END_DOC
  implicit none
  integer(bit_kind), intent(in)  :: key_i(N_int,2)
  double precision, intent(out)  :: hmono,herf,heff,hderiv,htot,hthree
  integer                        :: occ(N_int*bit_kind_size,2)
  integer                        :: Ne(2),i,j,ii,jj,ispin,jspin,k,kk
  double precision :: get_mo_two_e_integral_erf,mo_two_e_integral_eff_pot

  call bitstring_to_list_ab(key_i,occ,Ne,N_int)

  hmono = 0.d0
  herf  = 0.d0
  heff  = 0.d0
  hderiv= 0.d0
  hthree = 0.d0
  htot = 0.d0

  do ispin = 1, 2 
   do i = 1, Ne(ispin) ! 
    ii = occ(i,ispin) 
    hmono += mo_one_e_integrals(ii,ii)
   enddo
  enddo


  ! alpha/beta two-body
  ispin = 1
  jspin = 2 
  do i = 1, Ne(ispin) 
   ii = occ(i,ispin) 
   do j = 1, Ne(jspin)
    jj = occ(j,jspin) 
    herf += get_mo_two_e_integral_erf(ii,jj,ii,jj,mo_integrals_erf_map)
    heff += mo_two_e_integral_eff_pot(ii,jj,ii,jj) 
    hderiv += mo_non_hermit_term(ii,jj,ii,jj) 
   enddo
  enddo

  ! alpha/alpha two-body
  do i = 1, Ne(ispin)
   ii = occ(i,ispin) 
   do j = i+1, Ne(ispin)
    jj = occ(j,ispin) 
    herf += get_mo_two_e_integral_erf(ii,jj,ii,jj,mo_integrals_erf_map) - get_mo_two_e_integral_erf(ii,jj,jj,ii,mo_integrals_erf_map)
    heff += mo_two_e_integral_eff_pot(ii,jj,ii,jj) - mo_two_e_integral_eff_pot(ii,jj,jj,ii)
    hderiv += mo_non_hermit_term(ii,jj,ii,jj) & 
    - 0.5d0 * mo_non_hermit_term(ii,jj,jj,ii) - 0.5d0 * mo_non_hermit_term(jj,ii,ii,jj) 
   enddo
  enddo

  ! beta/beta two-body
  do i = 1, Ne(jspin)
   ii = occ(i,jspin) 
   do j = i+1, Ne(jspin)
    jj = occ(j,jspin) 
    herf += get_mo_two_e_integral_erf(ii,jj,ii,jj,mo_integrals_erf_map) - get_mo_two_e_integral_erf(ii,jj,jj,ii,mo_integrals_erf_map)
    heff += mo_two_e_integral_eff_pot(ii,jj,ii,jj) - mo_two_e_integral_eff_pot(ii,jj,jj,ii)
    hderiv += mo_non_hermit_term(ii,jj,ii,jj) & 
    - 0.5d0 * mo_non_hermit_term(ii,jj,jj,ii) & 
    - 0.5d0 * mo_non_hermit_term(jj,ii,ii,jj) 
   enddo
  enddo
 if(three_body_h_tc)then
   if(Ne(1)+Ne(2).ge.3)then
!!   ! alpha/alpha/beta three-body
!   print*,'In diag H mu for three body'
    do i = 1, Ne(1)
     ii = occ(i,1) 
     do j = i+1, Ne(1)
      jj = occ(j,1) 
      do k = 1, Ne(2)
       kk = occ(k,2) 
!      print*,kk,jj,ii,kk,jj,ii
!      print*,three_body_ints(kk,jj,ii,kk,jj,ii),three_body_ints(kk,jj,ii,kk,ii,jj)
       hthree += three_body_ints(kk,jj,ii,kk,jj,ii) - three_body_ints(kk,jj,ii,kk,ii,jj) 
      enddo
     enddo
    enddo
 
    ! beta/beta/alpha three-body
    do i = 1, Ne(2)
     ii = occ(i,2) 
     do j = i+1, Ne(2)
      jj = occ(j,2) 
      do k = 1, Ne(1)
       kk = occ(k,1) 
       hthree += three_body_ints(kk,jj,ii,kk,jj,ii) - three_body_ints(kk,ii,jj,kk,jj,ii) 
      enddo
     enddo
    enddo
 
!!   ! alpha/alpha/alpha three-body
!!   do i = 1, Ne(1)
!!    ii = occ(i,1) 
!!    do j = i+1, Ne(1)
!!     jj = occ(j,1) 
!!     do k = j+1, Ne(1)
!!      kk = occ(k,1) 
!!      ! < i j k | i j k > - < i j k | j i k > - < i j k | i k j > - < i j k | k j i >
!!      hthree += three_body_ints(ii,jj,kk,ii,jj,kk) - three_body_ints(ii,jj,kk,jj,ii,kk) & 
!!                                                   - three_body_ints(ii,jj,kk,ii,kk,jj) & 
!!                                                   - three_body_ints(ii,jj,kk,kk,jj,ii)
!!     enddo
!!    enddo
!!   enddo
 
!!   ! beta/beta/beta three-body
!!   do i = 2, Ne(2)
!!    ii = occ(i,2) 
!!    do j = i+2, Ne(2)
!!     jj = occ(j,2) 
!!     do k = j+2, Ne(2)
!!      kk = occ(k,2) 
!!      ! < i j k | i j k > - < i j k | j i k > - < i j k | i k j > - < i j k | k j i >
!!      hthree += three_body_ints(ii,jj,kk,ii,jj,kk) - three_body_ints(ii,jj,kk,jj,ii,kk) & 
!!                                                   - three_body_ints(ii,jj,kk,ii,kk,jj) & 
!!                                                   - three_body_ints(ii,jj,kk,kk,jj,ii)
!!     enddo
!!    enddo
!!   enddo
 
   endif
  endif
  htot = hmono + herf + heff + hderiv + hthree

end

subroutine single_htilde_mu_mat(key_j,key_i,hmono,herf,heff,hderiv,hthree,htot)
  use bitmasks
  BEGIN_DOC
! <key_j | H_tilde | key_i> for single excitation  
!!
!! WARNING !!
! 
! Non hermitian !!
  END_DOC
  implicit none
  integer(bit_kind), intent(in)  :: key_j(N_int,2),key_i(N_int,2)
  double precision, intent(out)  :: hmono,herf,heff,hderiv,hthree, htot
  integer                        :: occ(N_int*bit_kind_size,2)
  integer                        :: Ne(2),i,j,ii,jj,ispin,jspin,k,kk
  integer                        :: degree,exc(0:2,2,2)
  integer                        :: h1, p1, h2, p2, s1, s2
  double precision :: get_mo_two_e_integral_erf,mo_two_e_integral_eff_pot,phase
  integer :: other_spin(2)
  other_spin(1) = 2
  other_spin(2) = 1

  call get_excitation_degree(key_i,key_j,degree,N_int)

  hmono = 0.d0
  herf  = 0.d0
  heff  = 0.d0
  hderiv= 0.d0
  hthree = 0.d0
  htot = 0.d0

  if(degree.ne.1)then
   return
  endif
  call bitstring_to_list_ab(key_i,occ,Ne,N_int)
  call get_single_excitation(key_i,key_j,exc,phase,N_int)
  call decode_exc(exc,1,h1,p1,h2,p2,s1,s2)

  hmono = mo_one_e_integrals(h1,p1) * phase
  
  ! alpha/beta two-body 
  ispin = other_spin(s1)
  do i = 1, Ne(ispin)
   ii = occ(i,ispin) 
   herf   += get_mo_two_e_integral_erf(ii,p1,ii,h1,mo_integrals_erf_map)
   heff   += mo_two_e_integral_eff_pot(ii,p1,ii,h1) 
   hderiv += mo_non_hermit_term(ii,p1,ii,h1) 
  enddo
  ! same spin two-body 
  do i = 1, Ne(s1)
   ii = occ(i,s1) 
   ! (h1p1|ii ii) - (h1 ii| p1 ii)
   herf   += get_mo_two_e_integral_erf(ii,p1,ii,h1,mo_integrals_erf_map) - get_mo_two_e_integral_erf(ii,p1,h1,ii,mo_integrals_erf_map) 
   heff   += mo_two_e_integral_eff_pot(ii,p1,ii,h1) - mo_two_e_integral_eff_pot(ii,p1,h1,ii)
   hderiv += mo_non_hermit_term(ii,p1,ii,h1) & 
          -0.5d0 * mo_non_hermit_term(ii,p1,h1,ii) & 
          -0.5d0 * mo_non_hermit_term(p1,ii,ii,h1)
  enddo

 if(three_body_h_tc)then
  ! alpha/alpha/beta three-body
  if(Ne(1)+Ne(2).ge.3)then
   if(s1 == 2)then ! single beta 
    ! alpha-alpha + hole/particle beta 
    do i = 1, Ne(1)
     ii = occ(i,1) 
     do j = i+1, Ne(1)
      jj = occ(j,1) 
      !                         b  a a     b a a       b  a a   b a a
      !                       < h1 j  i | p1 j i > - < h1 j i | p1 i j >
      hthree += three_body_ints(h1,jj,ii,p1,jj,ii) - three_body_ints(h1,jj,ii,p1,ii,jj)
     enddo
    enddo

    ! alpha-beta + hole/particle beta
    do i = 1, Ne(1)
     ii = occ(i,1) 
     do j = 1, Ne(2)
      jj = occ(j,2) 
      !                         b  b a   b b a         b  b a   b b a
      !                       < h1 j  i | p1 j i > - < h1 j i | j p1 i >
      hthree += three_body_ints(h1,jj,ii,p1,jj,ii) - three_body_ints(h1,jj,ii,jj,p1,ii)
     enddo
    enddo

   else ! single alpha 
    ! beta-beta + hole/particle alpha 
    do i = 1, Ne(2)
     ii = occ(i,2) 
     do j = i+1, Ne(2)
      jj = occ(j,2)
      !                         a  b b   a  b b       a  b b   a  b b
      !                       < h1 j i | p1 j i > - < h1 j i | p1 i j >
      hthree += three_body_ints(h1,ii,jj,p1,ii,jj) - three_body_ints(h1,ii,jj,p1,jj,ii)
     enddo
    enddo
    ! alpha-beta + hole/particle alpha 
    do i = 1, Ne(2)
     ii = occ(i,2) 
     do j = 1, Ne(1)
      jj = occ(j,1)
      !                         a  a b   a  a b                       a  a b   a  a b   
      !                       < h1 j i | p1 j i >  -                < h1 j i | j p1 i >  
      hthree += three_body_ints(h1,jj,ii,p1,jj,ii) - three_body_ints(h1,jj,ii,jj,p1,ii)
     enddo
    enddo

   endif
  endif
 endif

  herf    *= phase
  heff    *= phase
  hderiv  *= phase
  hthree  *= phase
  htot = hmono + herf + heff + hderiv + hthree
end

subroutine double_htilde_mu_mat(key_j,key_i,hmono,herf,heff,hderiv,hthree,htot)
  use bitmasks
  BEGIN_DOC
! <key_j | H_tilde | key_i> for double excitation  
!!
!! WARNING !!
! 
! Non hermitian !!
  END_DOC
  implicit none
  integer(bit_kind), intent(in)  :: key_j(N_int,2),key_i(N_int,2)
  double precision, intent(out)  :: hmono,herf,heff,hderiv,hthree,htot
  integer                        :: occ(N_int*bit_kind_size,2)
  integer                        :: Ne(2),i,j,ii,jj,ispin,jspin,k,kk
  integer                        :: degree,exc(0:2,2,2)
  integer                        :: h1, p1, h2, p2, s1, s2
  double precision :: get_mo_two_e_integral_erf,mo_two_e_integral_eff_pot,phase
  integer :: other_spin(2)
  other_spin(1) = 2
  other_spin(2) = 1

  call bitstring_to_list_ab(key_i,occ,Ne,N_int)
  call get_excitation_degree(key_i,key_j,degree,N_int)

  hmono = 0.d0
  herf  = 0.d0
  heff  = 0.d0
  hderiv= 0.d0
  hthree = 0.d0
  htot = 0.d0

  if(degree.ne.2)then
   return
  endif
  call get_double_excitation(key_i,key_j,exc,phase,N_int)
  call decode_exc(exc,2,h1,p1,h2,p2,s1,s2)
  ! opposite spin two-body 
  herf    = get_mo_two_e_integral_erf(p1,p2,h1,h2,mo_integrals_erf_map)   
  heff    = mo_two_e_integral_eff_pot(p1,p2,h1,h2) 
  hderiv  = mo_non_hermit_term(p1,p2,h1,h2) 
  ! same spin two-body 
  if(s1.eq.s2)then
   herf   -= get_mo_two_e_integral_erf(p1,p2,h2,h1,mo_integrals_erf_map)   
   heff   -= mo_two_e_integral_eff_pot(p1,p2,h2,h1) 
   hderiv -= 0.5d0 * mo_non_hermit_term(p1,p2,h2,h1) & 
            +0.5d0 * mo_non_hermit_term(p2,p1,h1,h2)
  endif
  if(three_body_h_tc)then
   ! alpha/alpha/beta threee-body 
   if(Ne(1)+Ne(2).ge.3)then
    if(s1.eq.s2.and.s2.eq.1)then ! double alpha 
     do k = 1, Ne(2)
      kk = occ(k,2)
      hthree += three_body_ints(h1,h2,kk,p1,p2,kk) 
     enddo
    else if(s1.eq.s2.and.s2.eq.2)then ! double beta 
     do k = 1, Ne(1)
      kk = occ(k,1)
      hthree +=   three_body_ints(h1,h2,kk,p1,p2,kk) 
     enddo
    else ! double alpha/beta 
     if(s1.eq.1.and.s2.eq.2)then ! s1 == alpha , s2 == beta 
      do k = 1, Ne(1)
       kk = occ(k,1) ! direct - exchange in alpha 
       hthree +=  three_body_ints(h1,h2,kk,p1,p2,kk) - three_body_ints(h1,h2,kk,kk,p2,p1) 
      enddo
      do k = 1, Ne(2)
       kk = occ(k,2)! direct - exchange in beta 
       hthree +=  three_body_ints(h1,h2,kk,p1,p2,kk) - three_body_ints(h1,h2,kk,p1,kk,p2)
      enddo
     else if(s1.eq.2.and.s2.eq.1)then  ! s1 == beta, s2 == alpha 
      do k = 1, Ne(2)
       kk = occ(k,2) ! direct - exchange in beta 
       hthree +=  three_body_ints(h1,h2,kk,p1,p2,kk) - three_body_ints(h1,h2,kk,kk,p2,p1) 
      enddo
      do k = 1, Ne(1)
       kk = occ(k,1)! direct - exchange in alpha 
       hthree +=  three_body_ints(h1,h2,kk,p1,p2,kk) - three_body_ints(h1,h2,kk,p1,kk,p2)
      enddo
     endif 
    endif
   endif
  endif
  herf   *= phase
  heff   *= phase
  hderiv *= phase
  hthree  *= phase
  htot = herf + heff + hderiv + hthree
 end


subroutine htilde_mu_mat(key_j,key_i,hmono,herf,heff,hderiv,hthree,htot)
  use bitmasks
  BEGIN_DOC
! <key_j | H_tilde | key_i> 
!!
!! WARNING !!
! 
! Non hermitian !!
  END_DOC
  implicit none
  integer(bit_kind), intent(in)  :: key_j(N_int,2),key_i(N_int,2)
  double precision, intent(out)  :: hmono,herf,heff,hderiv,hthree,htot
  integer                        :: degree
   call get_excitation_degree(key_j,key_i,degree,N_int)
   hmono = 0.d0
   herf = 0.d0
   heff = 0.d0
   hderiv = 0.d0
   hthree = 0.d0
   htot = 0.d0
   if(degree.gt.3)then
    return
   else if(degree == 3.and.three_body_h_tc)then
    if(pure_three_body_h_tc)then
     call triple_htilde_mu_mat(key_j,key_i,hmono,herf,heff,hderiv,hthree,htot)
    endif
   else if(degree == 2)then
    call double_htilde_mu_mat(key_j,key_i,hmono,herf,heff,hderiv,hthree,htot)
   else if(degree == 1)then
    call single_htilde_mu_mat(key_j,key_i,hmono,herf,heff,hderiv,hthree,htot)
   else if(degree == 0)then
    call diag_htilde_mu_mat(key_i,hmono,herf,heff,hderiv,hthree,htot)
   endif
   
end

subroutine htilde_mat(key_j,key_i,hmono,herf,heff,hderiv,hthree,htot)
  use bitmasks
  BEGIN_DOC
! <key_j | H_tilde | key_i> 
!!
!! WARNING !!
! 
! Non hermitian !!
  END_DOC
  implicit none
  integer(bit_kind), intent(in)  :: key_j(N_int,2),key_i(N_int,2)
  double precision, intent(out)  :: hmono,herf,heff,hderiv,hthree,htot
  integer                        :: degree
  if(.not.ten_no_jastrow)then
   call htilde_mu_mat(key_j,key_i,hmono,herf,heff,hderiv,hthree,htot)
  else 
   call htilde_ten_no_mat(key_j,key_i,hmono,herf,heff,hderiv,hthree,htot)
  endif
end


subroutine triple_htilde_mu_mat(key_j,key_i,hmono,herf,heff,hderiv,hthree,htot)
  use bitmasks
  BEGIN_DOC
! <key_j | H_tilde | key_i> for triple excitation  
!!
!! WARNING !!
! 
! Genuine triple excitations of the same spin are not yet implemented
  END_DOC
  implicit none
  integer(bit_kind), intent(in)  :: key_j(N_int,2),key_i(N_int,2)
  double precision, intent(out)  :: hmono,herf,heff,hderiv,hthree,htot
  integer                        :: occ(N_int*bit_kind_size,2)
  integer                        :: Ne(2),i,j,ii,jj,ispin,jspin,k,kk
  integer                        :: degree,exc_double(0:2,2,2),exc_single(0:2,2,2)
  integer                        :: degree_alpha,degree_beta
  integer                        :: h1, p1, h2, p2, s1, s2, h3, p3, s3, h4, p4, s4
  double precision :: phase_double, phase_single
  double precision :: get_mo_two_e_integral_erf,mo_two_e_integral_eff_pot
  integer(bit_kind)  :: key_j_alpha(N_int,2),key_i_alpha(N_int,2)
  integer(bit_kind)  :: key_j_beta(N_int,2),key_i_beta(N_int,2)
  integer :: other_spin(2)
  other_spin(1) = 2
  other_spin(2) = 1

  call bitstring_to_list_ab(key_i,occ,Ne,N_int)
  call get_excitation_degree(key_i,key_j,degree,N_int)
  hmono = 0.d0
  herf = 0.d0
  heff = 0.d0
  hderiv = 0.d0
  hthree = 0.d0
  htot = 0.d0
  if(degree.ne.3)then
   return
  endif
  do i = 1, N_int
   key_j_alpha(i,1) = key_j(i,1)
   key_j_alpha(i,2) = 0_bit_kind
   key_i_alpha(i,1) = key_i(i,1)
   key_i_alpha(i,2) = 0_bit_kind

   key_j_beta(i,2) = key_j(i,2)
   key_j_beta(i,1) = 0_bit_kind
   key_i_beta(i,2) = key_i(i,2)
   key_i_beta(i,1) = 0_bit_kind
  enddo
  ! check whether it is a triple excitation of the same spin
  
  call get_excitation_degree(key_i_alpha,key_j_alpha,degree_alpha,N_int)
  call get_excitation_degree(key_i_beta,key_j_beta,degree_beta,N_int)
  if(degree_alpha==3.or.degree_beta==3)then
   return
  endif
  if(degree_alpha == 2.and.degree_beta == 1)then ! double alpha + single beta
   call get_double_excitation(key_i_alpha,key_j_alpha,exc_double,phase_double,N_int)
   call decode_exc(exc_double,2,h1,p1,h2,p2,s1,s2)
   call get_single_excitation(key_i_beta,key_j_beta,exc_single,phase_single,N_int)
   call decode_exc(exc_single,1,h3,p3,h4,p4,s3,s4)
  else if(degree_beta == 2 .and. degree_alpha == 1)then ! double beta + single alpha 
   call get_double_excitation(key_i_beta,key_j_beta,exc_double,phase_double,N_int)
   call decode_exc(exc_double,2,h1,p1,h2,p2,s1,s2)
   call get_single_excitation(key_i_alpha,key_j_alpha,exc_single,phase_single,N_int)
   call decode_exc(exc_single,1,h3,p3,h4,p4,s3,s4)
  else 
   print*,'PB !!'
   print*,'degree_beta, degree_alpha',degree_beta, degree_alpha
   print*,'degree',degree
   stop
  endif
  hthree = three_body_ints(h1,h2,h3,p1,p2,p3) -  three_body_ints(h1,h2,h3,p2,p1,p3)

  hthree  *= phase_single * phase_double
  htot = herf + heff + hderiv + hthree
 end

