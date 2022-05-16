program print_emu

  implicit none

  my_grid_becke = .True.
  my_n_pt_r_grid = 30
  my_n_pt_a_grid = 50
  touch  my_grid_becke my_n_pt_r_grid my_n_pt_a_grid

  read_wf = .True.
  touch read_wf

  PROVIDE N_det N_int
  call main()

end


subroutine main()

  use bitmasks

  implicit none
  integer                       :: k

  print *, 'mu_erf = ', mu_erf

  print *, ''
  do k = 1, nucl_num
    print *, 'k = ', k
    print *, j1b_gauss_pen(k)
    print *, nucl_coord(k,1:3)
  enddo

  return
end subroutine main()




