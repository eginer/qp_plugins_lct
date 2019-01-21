BEGIN_TEMPLATE

subroutine pt2_hcore ($arguments)
  use bitmasks
  implicit none
  $declarations
  
  BEGIN_DOC
  ! compute the standard Epstein-Nesbet perturbative first order coefficient and second order energetic contribution
  !
  ! for the various N_st states.
  !
  ! c_pert(i) = <psi(i)|H|det_pert>/( E(i) - <det_pert|H|det_pert> )
  !
  ! e_2_pert(i) = <psi(i)|H|det_pert>^2/( E(i) - <det_pert|H|det_pert> )
  !
  END_DOC
  
  integer                        :: i,j
  double precision               :: diag_H_mat_elem_fock, h
  double precision               :: i_H_psi_array(N_st)
  PROVIDE  selection_criterion

  ASSERT (Nint == N_int)
  ASSERT (Nint > 0)
  call i_h_core_psi_minilist(det_pert,minilist,idx_minilist,N_minilist,psi_selectors_coef,Nint,N_minilist,psi_selectors_size,N_st,i_H_psi_array)
  
  
  h = diag_H_mat_elem_fock(det_ref,det_pert,fock_diag_tmp,Nint)
  do i =1,N_st
    if(electronic_energy(i)>h.and.electronic_energy(i).ne.0.d0)then
      c_pert(i) = -1.d0
      e_2_pert(i) = selection_criterion*selection_criterion_factor*2.d0
    else if  (dabs(electronic_energy(i) - h) > 1.d-6) then
        c_pert(i) = i_H_psi_array(i) / (electronic_energy(i) - h)
        H_pert_diag(i) = h*c_pert(i)*c_pert(i)
        e_2_pert(i) = c_pert(i) * i_H_psi_array(i)
    else
      c_pert(i) = -1.d0
      e_2_pert(i) = -dabs(i_H_psi_array(i))
      H_pert_diag(i) = h
    endif
  enddo
  
end

SUBST [ arguments, declarations ]

electronic_energy,det_ref,det_pert,fock_diag_tmp,c_pert,e_2_pert,H_pert_diag,Nint,ndet,N_st,minilist,idx_minilist,N_minilist ;

    integer, intent(in)             :: Nint
    integer, intent(in)             :: ndet
    integer, intent(in)             :: N_st
    integer, intent(in)             :: N_minilist
    integer(bit_kind), intent(in)   :: det_ref (Nint,2)
    integer(bit_kind), intent(in)   :: det_pert(Nint,2)
    double precision , intent(in)   :: fock_diag_tmp(2,mo_num+1)
    double precision , intent(in)    :: electronic_energy(N_st)
    double precision , intent(out)  :: c_pert(N_st)
    double precision , intent(out)  :: e_2_pert(N_st)
    double precision, intent(out)   :: H_pert_diag(N_st)
    integer, intent(in)             :: idx_minilist(0:N_det_selectors)
    integer(bit_kind), intent(in)   :: minilist(Nint,2,N_det_selectors)
;;


END_TEMPLATE

! Note : If the arguments are changed here, they should also be changed accordingly in
! the perturbation.template.f file.

