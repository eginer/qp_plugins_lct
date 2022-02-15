subroutine i_soc_j(conf_i,conf_j,det_i, det_j,soc_ij)
 implicit none
 use bitmasks ! you need to include the bitmasks_module.f90 features
 integer(bit_kind), intent(in) :: det_i(N_int,2), det_j(N_int,2)
 integer(bit_kind), intent(in) :: conf_i(N_int,2), conf_j(N_int,2)
 double precision, intent(out) :: soc_ij
 BEGIN_DOC
! det_i and det_j are two determinants of possibly different ms
 END_DOC
 integer :: i,j, degree_doc, degree_somo,n_doc_i,n_doc_j,n_somo_i,n_somo_j
 integer(bit_kind) :: det_i_tmp(N_int), det_j_tmp(N_int)

 soc_ij = 0.d0
 
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
 
 if(degree_doc == 0)then
  if(degree_somo == 0)then
   ! same conf 
  else if(degree_somo == 1)then
   ! differ by a single exc of somo
  else 
   return
  endif
 else if(degree_doc == 1)then
  if(n_doc_i == n_doc_j .and. degree_somo == 0)then
  ! double excitation from a closed shell to a closed shell
  return
  else if(n_doc_i == n_doc_j .and. degree_somo == 1)then
  ! single excitation from a somo to another somo
  else if(abs(n_doc_i - n_doc_j) == 1 .and. abs(n_somo_i - n_somo_j)==2)then
  ! single excitation from a doubly occupied to an empty orbital 
  else if(abs(n_doc_i - n_doc_j) == 0 .and. abs(n_somo_i - n_somo_j)==0 .and. degree_somo == 1)then
  ! single excitation from a doubly occupied to a somo ==> transfer of doubly occupied orbitals 
  endif
 endif

end
