
 BEGIN_PROVIDER [double precision, htilde_matrix_elmt, (N_det,N_det)]
&BEGIN_PROVIDER [double precision, htilde_matrix_elmt_tranp, (N_det,N_det)]
&BEGIN_PROVIDER [double precision, htilde_matrix_elmt_eff, (N_det,N_det)]
&BEGIN_PROVIDER [double precision, htilde_matrix_elmt_deriv, (N_det,N_det)]
&BEGIN_PROVIDER [double precision, htilde_matrix_elmt_hcore, (N_det,N_det)]
&BEGIN_PROVIDER [double precision, htilde_matrix_elmt_hthree, (N_det,N_det)]

  BEGIN_DOC
  ! htilde_matrix_elmt(j,i) = <J| H^tilde |I> 
  !
  ! WARNING !!!!!!!!! IT IS NOT HERMITIAN !!!!!!!!!
  END_DOC
 
  implicit none
  integer          :: i, j
  double precision :: hmono,heff,hderiv,hthree,htot

  PROVIDE N_int
 
  if(zero_tc_eff_map) then

    htilde_matrix_elmt = H_matrix_all_dets
    htilde_matrix_elmt_eff = 0.d0
    htilde_matrix_elmt_deriv = 0.d0
    htilde_matrix_elmt_hcore = 0.d0
    htilde_matrix_elmt_hthree = 0.d0
    do i = 1, N_det
      do j = 1, N_det
        ! < J | Htilde | I >
        call htilde_mu_mat(psi_det(1,1,j), psi_det(1,1,i), N_int, hmono, heff, hderiv, hthree, htot)
        htilde_matrix_elmt(j,i) += hderiv + hthree
      enddo
    enddo

  else

    do i = 1, N_det
      do j = 1, N_det
        ! < J | Htilde | I >
        call htilde_mu_mat(psi_det(1,1,j), psi_det(1,1,i), N_int, hmono, heff, hderiv, hthree, htot)
        htilde_matrix_elmt(j,i)        = htot
        htilde_matrix_elmt_eff(j,i)    = heff
        htilde_matrix_elmt_deriv(j,i)  = hderiv
        htilde_matrix_elmt_hcore(j,i)  = hmono
        htilde_matrix_elmt_hthree(j,i) = hthree
      enddo
    enddo

  endif

  do i = 1, N_det
    do j = 1, N_det
      htilde_matrix_elmt_tranp(j,i) = htilde_matrix_elmt(i,j)
    enddo
  enddo

END_PROVIDER 




BEGIN_PROVIDER [ double precision, diag_htilde, (N_det)]

  implicit none
  integer          :: i
  double precision :: hmono, heff, hderiv, hthree, htot
 
  PROVIDE N_int

  do i = 1, N_det
    call htilde_mu_mat(psi_det(1,1,i), psi_det(1,1,i), N_int, hmono, heff, hderiv, hthree, htot)
    diag_htilde(i) = htot
  enddo

END_PROVIDER 

