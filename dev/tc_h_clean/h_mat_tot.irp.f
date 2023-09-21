
subroutine htilde_mu_mat_tot(key_j, key_i, Nint, htot)

  BEGIN_DOC
  ! <key_j | H_tilde | key_i> 
  !!
  !! WARNING !!
  ! 
  ! Non hermitian !!
  END_DOC

  use bitmasks

  implicit none
  integer, intent(in)           :: Nint
  integer(bit_kind), intent(in) :: key_j(Nint,2),key_i(Nint,2)
  double precision, intent(out) :: htot
  double precision              :: hmono, heff, hderiv, hthree

  call htilde_mu_mat(key_j, key_i, Nint, hmono, heff, hderiv, hthree, htot)
  htot = hmono + heff + hderiv + hthree

end subroutine htilde_mu_mat_tot



subroutine htilde_mu_mat(key_j, key_i, Nint, hmono, heff, hderiv, hthree, htot)

  BEGIN_DOC
  ! <key_j | H_tilde | key_i> 
  !!
  !! WARNING !!
  ! 
  ! Non hermitian !!
  END_DOC

  use bitmasks
  implicit none
  integer,           intent(in) :: Nint
  integer(bit_kind), intent(in) :: key_j(Nint,2), key_i(Nint,2)
  double precision, intent(out) :: hmono, heff, hderiv, hthree, htot
  integer                       :: degree

  call get_excitation_degree(key_j, key_i, degree, Nint)
  hmono  = 0.d0
  heff   = 0.d0
  hderiv = 0.d0
  hthree = 0.d0
  htot   = 0.d0


  if(degree.gt.3)then
    return
  else if(degree == 2) then
    call double_htilde_mu_mat_scal_map(Nint, key_j, key_i, hmono, heff, hderiv, htot)
  else if(degree == 1) then
    call single_htilde_mu_mat_scal_map(Nint, key_j, key_i, hmono, heff, hderiv, htot)
  else if(degree == 0) then
    call diag_htilde_mu_mat_scal_map(Nint, key_i, hmono, heff, hderiv, htot)
  endif

  if(three_body_h_tc) then
    if(degree == 2) then
      if(.not.double_normal_ord) then
        call double_htilde_mu_mat_three_body(Nint, key_j, key_i, hthree)
      endif
    else if(degree == 1)then
      call single_htilde_mu_mat_three_body(Nint, key_j, key_i, hthree)
    else if(degree == 0)then
      call diag_htilde_mu_mat_three_body(Nint, key_i, hthree)
    endif
  endif

  htot += hthree
   
end subroutine htilde_mu_mat




! -------------------------------------------------------------------------------------------------

! ---

subroutine htildedag_mu_mat_tot(key_j, key_i, Nint, htot)

  BEGIN_DOC
  ! <key_j | H_tilde_dag | key_i> 
  !!
  !! WARNING !!
  ! 
  ! Non hermitian !!
  END_DOC

  use bitmasks

  implicit none
  integer,           intent(in) :: Nint
  integer(bit_kind), intent(in) :: key_j(Nint,2), key_i(Nint,2)
  double precision, intent(out) :: htot
  double precision              :: hmono, heff, hderiv, hthree

  call htildedag_mu_mat(key_j, key_i, Nint, hmono, heff, hderiv, hthree, htot)
  htot = hmono + heff + hderiv + hthree

end subroutine htildedag_mu_mat_tot

! ---

subroutine htildedag_mu_mat(key_j, key_i, Nint, hmono, heff, hderiv, hthree, htot)

  BEGIN_DOC
  ! <key_j | H_tilde_dag | key_i> 
  !!
  !! WARNING !!
  ! 
  ! Non hermitian !!
  END_DOC

  use bitmasks

  implicit none
  integer,           intent(in) :: Nint
  integer(bit_kind), intent(in) :: key_j(Nint,2), key_i(Nint,2)
  double precision, intent(out) :: hmono, heff, hderiv, hthree, htot
  integer                       :: degree

  hmono  = 0.d0
  heff   = 0.d0
  hderiv = 0.d0
  hthree = 0.d0
  htot   = 0.d0

  call get_excitation_degree(key_j, key_i, degree, Nint)

  if(degree.gt.3) then
    return
  else if(degree == 2) then
    call double_htildedag_mu_mat_scal_map(Nint, key_j, key_i, hmono, heff, hderiv, htot)
  else if(degree == 1) then
    call single_htildedag_mu_mat_scal_map(Nint, key_j, key_i, hmono, heff, hderiv, htot)
  else if(degree == 0) then
    call diag_htildedag_mu_mat_scal_map(Nint, key_i, hmono, heff, hderiv, htot)
  endif

  if(three_body_h_tc) then
    if(degree == 2) then
      if(.not.double_normal_ord) then
        call double_htilde_mu_mat_three_body(Nint, key_j, key_i, hthree)
        hthree = -hthree !dag
      endif
    else if(degree == 1) then
      call single_htilde_mu_mat_three_body(Nint, key_j, key_i, hthree)
    else if(degree == 0) then
      call diag_htilde_mu_mat_three_body(Nint, key_i, hthree)
    endif
  endif
  
  htot += hthree
   
end subroutine htildedag_mu_mat

! ---

! -------------------------------------------------------------------------------------------------

subroutine hji_hij_mu_mat_tot(key_j, key_i, Nint, htot_ji,htot_ij)

  BEGIN_DOC
  ! htot_ji = <key_j | H_tilde | key_i> 
  !!
  ! htot_ij = <key_i | H_tilde | key_j> 
  !!
  !! WARNING !!
  ! 
  ! Non hermitian !!
  END_DOC

  use bitmasks

  implicit none
  integer, intent(in)           :: Nint
  integer(bit_kind), intent(in) :: key_j(Nint,2),key_i(Nint,2)
  double precision, intent(out) :: htot_ji,htot_ij
  double precision              :: hmono_ij, hmono_ji,heff, hderiv_ij, hderiv_ji, hthree

  call hji_hij_mu_mat(key_j, key_i, Nint, hmono_ji, hmono_ij, heff, hderiv_ji, hderiv_ij, hthree, htot_ji, htot_ij)
  htot_ji = hmono_ji + heff + hderiv_ji + hthree
  htot_ij = hmono_ij + heff + hderiv_ij + hthree

end subroutine htilde_mu_mat_tot



subroutine hji_hij_mu_mat(key_j, key_i, Nint, hmono_ji, hmono_ij, heff, hderiv_ji, hderiv_ij, hthree, htot_ji, htot_ij)
                
  BEGIN_DOC
  ! htot_ji = <key_j | H_tilde | key_i> 
  !!
  ! htot_ij = <key_i | H_tilde | key_j> 
  !!
  !! WARNING !!
  ! 
  ! Non hermitian !!
  END_DOC

  use bitmasks
  implicit none
  integer,           intent(in) :: Nint
  integer(bit_kind), intent(in) :: key_j(Nint,2), key_i(Nint,2)
  double precision, intent(out) :: hmono_ij, hmono_ji,heff, hderiv_ij, hderiv_ji, htot_ij, htot_ji, hthree
  integer                       :: degree

  call get_excitation_degree(key_j, key_i, degree, Nint)

  hmono_ji  = 0.d0
  hmono_ij  = 0.d0
  heff      = 0.d0
  hderiv_ji = 0.d0
  hderiv_ij = 0.d0
  hthree    = 0.d0
  htot_ji   = 0.d0
  htot_ij   = 0.d0


  if(degree.gt.3)then
    return
  else if(degree == 2) then
    call double_hij_ji_mu_mat_scal_map(Nint, key_j, key_i, hmono_ji, heff, hderiv_ji, hderiv_ij, htot_ij, htot_ji)
    hmono_ij = hmono_ji
  else if(degree == 1) then
    call single_hij_ji_mu_mat_scal_map(Nint, key_j, key_i, hmono_ji,hmono_ij, heff, hderiv_ji, hderiv_ij, htot_ij, htot_ji)
  else if(degree == 0) then
    call diag_htilde_mu_mat_scal_map(Nint, key_i, hmono_ji, heff, hderiv_ji, htot_ji)
    hmono_ij  = hmono_ji
    hderiv_ij = hderiv_ji
    htot_ij   = htot_ji
  endif

  if(three_body_h_tc) then
    if(degree == 2) then
      if(.not.double_normal_ord) then
        call double_htilde_mu_mat_three_body(Nint, key_j, key_i, hthree)
      endif
    else if(degree == 1)then
      call single_htilde_mu_mat_three_body(Nint, key_j, key_i, hthree)
    else if(degree == 0)then
      call diag_htilde_mu_mat_three_body(Nint, key_i, hthree)
    endif
  endif

  htot_ij += hthree
  htot_ji += hthree
   
end subroutine hji_hij_mu_mat
