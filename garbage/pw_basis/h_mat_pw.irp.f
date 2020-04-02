BEGIN_PROVIDER  [double precision, pot_mat_pw, (n_pw_basis,n_pw_basis) ]
 implicit none
 include 'constants.include.F'
 integer :: i,j,m,n
 pot_mat_pw = 0.d0
 double precision  :: v_diag_pw,diag,k_i(3), k_j(3),k_ij
 integer :: nx
 nx = 1000000
 diag = v_diag_pw(pw_box_x,pw_box_y,pw_box_z,nx)
 diag = diag 
 do i = 1, n_pw_basis
  pot_mat_pw(i,i) = diag
 enddo
 do i = 1, n_pw_basis
  do m = 1, 3
   k_i(m) = d_list_pw(m,i) * delta_k_pw(m)
  enddo
  do j = i+1, n_pw_basis
   do m = 1, 3
    k_j(m) = d_list_pw(m,j) * delta_k_pw(m)
   enddo
   k_ij = 0.d0
   do m = 1, 3
    k_ij += (k_i(m) - k_j(m))**2.d0
   enddo
   pot_mat_pw(i,j) = dfour_pi / k_ij
   pot_mat_pw(j,i) = dfour_pi / k_ij
  enddo
 enddo
 pot_mat_pw *= -nucl_charge(1) / vol_box_pw

END_PROVIDER 

BEGIN_PROVIDER [double precision, kin_mat_pw, (n_pw_basis)]
 implicit none
 integer :: i,m
 double precision :: k_i
 do i = 1, n_pw_basis
  k_i = 0.d0
  do m = 1, 3
   k_i += ( d_list_pw(m,i) * delta_k_pw(m) )**2.d0
  enddo
  kin_mat_pw(i) = 0.5d0 * k_i 
 enddo

END_PROVIDER

BEGIN_PROVIDER [double precision, h_mat_pw, (n_pw_basis,n_pw_basis) ]
 implicit none
 integer :: i
 h_mat_pw = pot_mat_pw * lambda_v_pw
 do i = 1, n_pw_basis
  h_mat_pw(i,i) += kin_mat_pw(i)
 enddo

END_PROVIDER 

 BEGIN_PROVIDER [double precision, pw_eigval, (n_pw_basis) ]
&BEGIN_PROVIDER [double precision, pw_eigvec, (n_pw_basis,n_pw_basis) ]
 implicit none
 integer :: i
 call lapack_diagd(pw_eigval,pw_eigvec,h_mat_pw,n_pw_basis,n_pw_basis)
 do i = 1, n_pw_basis
  write(36,*)dble(i)/dble(n_pw_basis),pw_eigval(i)
 enddo
END_PROVIDER 
