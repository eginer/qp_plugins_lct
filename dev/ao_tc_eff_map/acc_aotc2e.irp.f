program acc_aotc2e

  implicit none

  !print *, ' coul integrals:'
  !call test_coul()

  !print *, ' erf integrals:'
  !call test_erf()

  call test_coulerf()

end


subroutine test_coul()

  implicit none 

  integer          :: iao, jao, kao, lao
  double precision :: integral_1, integral_2, delta, accu_abs
  double precision :: t1, t2, t_old, t_new

  double precision :: j1b_gauss_coul, j1b_gauss_coul_acc

  accu_abs = 0.d0
  t_old    = 0.d0
  t_new    = 0.d0

  do iao = 1, ao_num ! r1
    do jao = 1, ao_num ! r2
      do kao = 1, ao_num ! r1
        do lao = 1, ao_num ! r2

          CALL CPU_TIME(t1)
          integral_1 = j1b_gauss_coul(iao, kao, jao, lao)
          CALL CPU_TIME(t2)
          t_old += t2 - t1

          CALL CPU_TIME(t1)
          integral_2 = j1b_gauss_coul_acc(iao, kao, jao, lao)
          CALL CPU_TIME(t2)
          t_new += t2 - t1

          delta = dabs( integral_1 - integral_2 )
          if( delta .gt. 1.d-10 ) then
            print *, ' WARNING !!!!!!!!!!!' 
            print *, iao, jao, kao, lao
            print *, integral_1, integral_2, delta
          endif

          accu_abs += delta
        enddo
      enddo
    enddo
  enddo

  print *, ' accu_abs = ', accu_abs
  print *, ' t_old = ', t_old
  print *, ' t_new = ', t_new

  return
end subroutine test_coul


subroutine test_erf()

  implicit none 

  integer          :: iao, jao, kao, lao
  double precision :: integral_1, integral_2, delta, accu_abs
  double precision :: t1, t2, t_old, t_new

  double precision :: j1b_gauss_erf, j1b_gauss_erf_acc

  accu_abs = 0.d0
  t_old    = 0.d0
  t_new    = 0.d0

  do iao = 1, ao_num ! r1
    do jao = 1, ao_num ! r2
      do kao = 1, ao_num ! r1
        do lao = 1, ao_num ! r2

          CALL CPU_TIME(t1)
          integral_1 = j1b_gauss_erf(iao, kao, jao, lao)
          CALL CPU_TIME(t2)
          t_old += t2 - t1

          CALL CPU_TIME(t1)
          integral_2 = j1b_gauss_erf_acc(iao, kao, jao, lao)
          CALL CPU_TIME(t2)
          t_new += t2 - t1

          delta = dabs( integral_1 - integral_2 )
          if( delta .gt. 1.d-10 ) then
            print *, ' WARNING !!!!!!!!!!!' 
            print *, iao, jao, kao, lao
            print *, integral_1, integral_2, delta
          endif

          accu_abs += delta
        enddo
      enddo
    enddo
  enddo

  print *, ' accu_abs = ', accu_abs
  print *, ' t_old = ', t_old
  print *, ' t_new = ', t_new

  return
end subroutine test_erf



subroutine test_coulerf()

  implicit none 

  integer          :: iao, jao, kao, lao
  double precision :: integral_1, integral_2, integral_3, integral_4
  double precision :: accu_abs1, accu_abs2, accu_abs3
  double precision :: t1, t2, tot1, tot2, tot3, tot4

  double precision :: j1b_gauss_coul, j1b_gauss_coul_acc
  double precision :: j1b_gauss_erf, j1b_gauss_erf_acc
  double precision :: j1b_gauss_coulerf
  double precision :: j1b_gauss_coulerf_schwartz

  accu_abs1 = 0.d0
  accu_abs2 = 0.d0
  accu_abs3 = 0.d0

  tot1      = 0.d0
  tot2      = 0.d0
  tot3      = 0.d0
  tot4      = 0.d0

  do iao = 1, ao_num ! r1
    do jao = 1, ao_num ! r2
      do kao = 1, ao_num ! r1
        do lao = 1, ao_num ! r2

          integral_1 = 0.d0
          CALL CPU_TIME(t1)
          integral_1 += j1b_gauss_coul(iao, kao, jao, lao)
          integral_1 += j1b_gauss_erf (iao, kao, jao, lao)
          CALL CPU_TIME(t2)
          tot1 += t2 - t1

          integral_2 = 0.d0
          CALL CPU_TIME(t1)
          integral_2 += j1b_gauss_coul_acc(iao, kao, jao, lao)
          integral_2 += j1b_gauss_erf_acc (iao, kao, jao, lao)
          accu_abs1  += dabs( integral_2 - integral_1 )
          CALL CPU_TIME(t2)
          tot2 += t2 - t1

          CALL CPU_TIME(t1)
          integral_3 = j1b_gauss_coulerf(iao, kao, jao, lao)
          accu_abs2 = accu_abs2 + dabs( integral_3 - integral_1 )
          CALL CPU_TIME(t2)
          tot3 += t2 - t1

          CALL CPU_TIME(t1)
          integral_4 = j1b_gauss_coulerf_schwartz(iao, kao, jao, lao)
          accu_abs3 += dabs( integral_4 - integral_1 )
          CALL CPU_TIME(t2)
          tot4 += t2 - t1

        enddo
      enddo
    enddo
  enddo

  print *, ' accu_abs1 = ', accu_abs1
  print *, ' accu_abs2 = ', accu_abs2
  print *, ' accu_abs3 = ', accu_abs3

  print *, ' tot1 = ', tot1
  print *, ' tot2 = ', tot2
  print *, ' tot3 = ', tot3
  print *, ' tot4 = ', tot4

  return
end subroutine test_coulerf
