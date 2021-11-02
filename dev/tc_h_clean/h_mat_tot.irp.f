subroutine htilde_mu_mat_tot(key_j,key_i,Nint,htot)
  use bitmasks
  BEGIN_DOC
! <key_j | H_tilde | key_i> 
!!
!! WARNING !!
! 
! Non hermitian !!
  END_DOC
  implicit none
  integer(bit_kind), intent(in)  :: key_j(Nint,2),key_i(Nint,2)
  integer, intent(in)            :: Nint
  double precision, intent(out)  :: htot
  double precision  :: hmono,heff,hderiv,hthree
  call htilde_mu_mat(key_j,key_i,hmono,heff,hderiv,hthree,htot)
  htot = hmono+heff+hderiv+hthree
end

subroutine htilde_mu_mat(key_j,key_i,hmono,heff,hderiv,hthree,htot)
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
  double precision, intent(out)  :: hmono,heff,hderiv,hthree,htot

  integer                        :: degree
   call get_excitation_degree(key_j,key_i,degree,N_int)
   hmono = 0.d0
   heff = 0.d0
   hderiv = 0.d0
   hthree = 0.d0
   htot = 0.d0
   if(degree.gt.3)then
    return
   else if(degree == 2)then
    call double_htilde_mu_mat_scal_map(key_j,key_i,hmono,heff,hderiv,htot)
   else if(degree == 1)then
    call single_htilde_mu_mat_scal_map(key_j,key_i,hmono,heff,hderiv,htot)
   else if(degree == 0)then
    call diag_htilde_mu_mat_scal_map(key_i,hmono,heff,hderiv,htot)
   endif
   if(three_body_h_tc)then
    if(degree == 2)then
     if(.not.double_normal_ord.and.double_3_body_tc)then
       call double_htilde_mu_mat_three_body(key_j,key_i,hthree)
     endif
    else if(degree == 1)then
     call single_htilde_mu_mat_three_body(key_j,key_i,hthree)
    else if(degree == 0)then
     call diag_htilde_mu_mat_three_body(key_i,hthree)
    endif
   endif
   htot += hthree
   
end
