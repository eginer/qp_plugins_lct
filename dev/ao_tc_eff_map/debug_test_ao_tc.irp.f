program debug_test_ao_tc

  implicit none
  call main()

end


subroutine main()

  implicit none 

  integer          :: j, num_A, power_A(3) 
  double precision :: A_center(3), C_center(3), gama

  do j = 1, ao_num
    num_A         = ao_nucl(j)
    power_A(1:3)  = ao_power(j,1:3)
    A_center(1:3) = nucl_coord(num_A,1:3)
    print *, j, num_A
  enddo

  do j = 1, nucl_num
    num_A         = ao_nucl(j)
    gama          = j1b_gauss_pen(j)
    !C_center(1:3) = nucl_coord(j,1:3)
    print *, j, num_A
  enddo

  return
end subroutine main
