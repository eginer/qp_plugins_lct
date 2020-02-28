
BEGIN_PROVIDER [double precision, a_resp_hf_mat, (n_singles_resp,n_singles_resp)
 implicit none
 integer :: i_single,j_single,i_orb,a_orb,j_orb,b_orb
 double precision :: eps_i,eps_a,get_two_e_integral,int_hartree,int_exchange
 PROVIDE mo_two_e_integrals_in_map mo_integrals_map big_array_exchange_integrals
 a_resp_hf_mat = 0.d0
 do i_single = 1, n_singles_resp
  i_orb = list_singles(i_single,1) 
  a_orb = list_singles(i_single,2) 
  eps_i = fock_matrix_diag_mo(i_orb)
  eps_a = fock_matrix_diag_mo(a_orb)
  do j_single = i_single, n_singles_resp
   j_orb = list_singles(j_single,1) 
   b_orb = list_singles(j_single,2) 
   !                                  <aj|ib>
   int_hartree  = get_two_e_integral(a_orb,j_orb,i_orb,b_orb,mo_integrals_map)
   !                                  <aj|bi>
   int_exchange = get_two_e_integral(a_orb,j_orb,b_orb,i_orb,mo_integrals_map)
   a_resp_hf_mat(j_single,i_single) = 2.d0 * int_hartree - int_exchange 
   a_resp_hf_mat(i_single,j_single) = a_resp_hf_mat(j_single,i_single)
  enddo
  ! special case for the diagonal 
  a_resp_hf_mat(i_single,i_single) += eps_a - eps_i
 enddo

END_PROVIDER 

BEGIN_PROVIDER [double precision, eigval_a_resp_h_mat, (n_singles_resp)]
 implicit none
 double precision :: eigvectors(n_singles_resp,n_singles_resp)
 call lapack_diagd(eigval_a_resp_h_mat,eigvectors,a_resp_hf_mat,n_singles_resp,n_singles_resp)
END_PROVIDER 


BEGIN_PROVIDER [double precision, b_resp_hf_mat, (n_singles_resp,n_singles_resp)
 implicit none
 integer :: i_single,j_single,i_orb,a_orb,j_orb,b_orb
 double precision :: eps_i,eps_a,get_two_e_integral,int_hartree,int_exchange
 PROVIDE mo_two_e_integrals_in_map mo_integrals_map big_array_exchange_integrals
 b_resp_hf_mat = 0.d0
 do i_single = 1, n_singles_resp
  i_orb = list_singles(i_single,1) 
  a_orb = list_singles(i_single,2) 
  eps_i = fock_matrix_diag_mo(i_orb)
  eps_a = fock_matrix_diag_mo(a_orb)
  do j_single = i_single, n_singles_resp
   j_orb = list_singles(j_single,1) 
   b_orb = list_singles(j_single,2) 
   !                                  <aj|ib>
   int_hartree  = get_two_e_integral(a_orb,j_orb,i_orb,b_orb,mo_integrals_map)
   !                                  <ab|ji>
   int_exchange = get_two_e_integral(a_orb,b_orb,j_orb,i_orb,mo_integrals_map)
   b_resp_hf_mat(j_single,i_single) = 2.d0 * int_hartree - int_exchange 
   b_resp_hf_mat(i_single,j_single) = b_resp_hf_mat(j_single,i_single)
  enddo
 enddo
END_PROVIDER 

BEGIN_PROVIDER [double precision, ab_resp_hf_mat, (2*n_singles_resp,2*n_singles_resp)]
 implicit none
 ab_resp_hf_mat = 0.d0 
 ab_resp_hf_mat( 1:n_singles_resp , 1:n_singles_resp ) = a_resp_hf_mat(1:n_singles_resp,1:n_singles_resp)
 ab_resp_hf_mat( 1:n_singles_resp, n_singles_resp+1:2*n_singles_resp ) = b_resp_hf_mat(1:n_singles_resp,1:n_singles_resp)
 ab_resp_hf_mat( n_singles_resp+1:2*n_singles_resp, 1:n_singles_resp ) = - b_resp_hf_mat(1:n_singles_resp,1:n_singles_resp)
 ab_resp_hf_mat( n_singles_resp+1:2*n_singles_resp, n_singles_resp+1:2*n_singles_resp) = - a_resp_hf_mat(1:n_singles_resp,1:n_singles_resp)
 print*,ab_resp_hf_mat(1:n_singles_resp,1:n_singles_resp)
END_PROVIDER 

 BEGIN_PROVIDER [double precision, eigval_ab_resp_hf_mat, (2*n_singles_resp) ]
&BEGIN_PROVIDER [double precision, eigvect_r_ab_resp_hf_mat, (2*n_singles_resp,2*n_singles_resp) ]
 implicit none
 double precision :: eigvalues_r(2*n_singles_resp),eigvalues_i(2*n_singles_resp)
 double precision :: eigvectors_r(2*n_singles_resp,2*n_singles_resp)
 integer :: i

 double precision :: mat(2*n_singles_resp,2*n_singles_resp)
 mat = 0.d0
 mat( 1:n_singles_resp , 1:n_singles_resp ) = a_resp_hf_mat(1:n_singles_resp,1:n_singles_resp)
 mat( 1:n_singles_resp, n_singles_resp+1:2*n_singles_resp ) = b_resp_hf_mat(1:n_singles_resp,1:n_singles_resp)
 mat( n_singles_resp+1:2*n_singles_resp, 1:n_singles_resp ) = - b_resp_hf_mat(1:n_singles_resp,1:n_singles_resp)
 mat( n_singles_resp+1:2*n_singles_resp, n_singles_resp+1:2*n_singles_resp) = - a_resp_hf_mat(1:n_singles_resp,1:n_singles_resp) 
! call non_sym_lapack_diagd(eigvalues_r,eigvalues_i,eigvectors_r,ab_resp_hf_mat,2*n_singles_resp,2*n_singles_resp)
 call non_sym_lapack_diagd(eigvalues_r,eigvalues_i,eigvectors_r,mat,2*n_singles_resp,2*n_singles_resp)
 do i = 1, 2 *n_singles_resp
  print*,'r/i', eigvalues_r(i),eigvalues_i(i)
 enddo

END_PROVIDER 

