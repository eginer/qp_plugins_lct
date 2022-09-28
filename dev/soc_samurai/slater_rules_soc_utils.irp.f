subroutine get_good_things_dets(conf_i,conf_j,det_i, det_j,n_doc_i,n_doc_j,degree_doc,n_somo_i,n_somo_j,degree_somo)
 implicit none
 use bitmasks ! you need to include the bitmasks_module.f90 features
 integer(bit_kind), intent(in) :: det_i(N_int,2), det_j(N_int,2)
 integer(bit_kind), intent(in) :: conf_i(N_int,2), conf_j(N_int,2)
 integer          , intent(out):: n_doc_i,n_doc_j,degree_doc,n_somo_i,n_somo_j,degree_somo
 ! degree of excitation of doubly occupied part 
 integer :: i
 degree_doc = 0
 n_doc_i = 0
 n_doc_j = 0
 do i = 1, N_int
  degree_doc += popcnt(xor(conf_i(i,1),conf_j(i,1)))
  n_doc_i += popcnt(conf_i(i,1))
  n_doc_j += popcnt(conf_j(i,1))
 enddo
 degree_doc = degree_doc / 2
 ! degree of excitation of single occupied part 
 degree_somo = 0
 do i = 1, N_int
  degree_somo += popcnt(xor(conf_i(i,2),conf_j(i,2)))
  n_somo_i += popcnt(conf_i(i,2))
  n_somo_j += popcnt(conf_j(i,2))
 enddo
 degree_somo = degree_somo / 2

end

subroutine get_diff_dets_somo_doc(det_i,det_j,conf_i,conf_j, somo_i, somo_j, somo_dif_ab, conf_dif)
 implicit none
 use bitmasks ! you need to include the bitmasks_module.f90 features
 integer(bit_kind), intent(in) :: det_i(N_int,2), det_j(N_int,2)
 integer(bit_kind), intent(in) :: conf_i(N_int,2), conf_j(N_int,2)
 integer(bit_kind), intent(out):: somo_i(N_int,2), somo_j(N_int,2)
 integer(bit_kind), intent(out):: somo_dif_ab(N_int,2), conf_dif(N_int,2)
 integer :: i
 do i = 1, N_int
  somo_i(i,1) = iand(det_i(i,1),conf_i(i,2)) ! extract alpha part of somo for det_i
  somo_i(i,2) = iand(det_i(i,2),conf_i(i,2)) ! extract beta  part of somo for det_i
  somo_j(i,1) = iand(det_j(i,1),conf_j(i,2)) ! extract alpha part of somo for det_j
  somo_j(i,2) = iand(det_j(i,2),conf_j(i,2)) ! extract beta  part of somo for det_j
 enddo
 do i = 1, N_int
  somo_dif_ab(i,1) = xor(somo_i(i,1),somo_j(i,1)) ! extract alpha differences in somo
  somo_dif_ab(i,2) = xor(somo_i(i,2),somo_j(i,2)) ! extract beta  differences in somo
  conf_dif(i,1) = xor(conf_i(i,1),conf_j(i,1))    ! difference in doubly occupied 
  conf_dif(i,2) = xor(conf_i(i,2),conf_j(i,2))    ! difference in singly occupied 
 enddo
end
