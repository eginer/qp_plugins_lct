program test_energy_integral_lda
 implicit none
 BEGIN_DOC
 !
 END_DOC

 read_wf = .true.
 touch read_wf
 basis_cor_func = "su_pbe_ot"
 touch basis_cor_func

 ! total one-e integrals 
 io_mo_one_e_integrals = "None"
 touch io_mo_one_e_integrals  
 ! Vne integrals on the MO basis 
 io_mo_integrals_n_e = "None"
 touch io_mo_integrals_n_e
 ! kinetic integrals on the MO basis 
 io_mo_integrals_kinetic = "None"
 touch io_mo_integrals_kinetic 
 ! Vne integrals on the AO basis 
 io_ao_integrals_n_e = "None"
 touch io_ao_integrals_n_e 
 ! kinetic integrals on the AO basis 
 io_ao_integrals_kinetic = "None"
 touch io_ao_integrals_kinetic 

   print*, '___________LDA___________'
   print*, "ex_lda [n]                        =", ex_lda_at_n
   print*, "ex_lda [n + delta_n]              =", ex_lda_at_n_plus_delta_n
   print*, "ex_lda [n + delta_n] - ex_lda [n] =", ex_lda_at_n_plus_delta_n - ex_lda_at_n
   print*, "n variation integral              :", int_vx_lda_at_n 

   print*, '___________PBE___________'
   print*, "ex_pbe [n] =", ex_pbe_at_n
   print*, "ex_pbe [n + delta_n]              =", ex_pbe_at_n_plus_delta_n
   print*, "ex_pbe [n + delta_n] - ex_pbe [n] =", ex_pbe_at_n_plus_delta_n - ex_pbe_at_n
   print*, "n variation integral              :", int_vx_pbe_at_n 
   print*, "ec_pbe [n]                        =", ec_pbe_at_n
   print*, "ec_pbe [n + delta_n]              =", ec_pbe_at_n_plus_delta_n
   print*, "ec_pbe [n + delta_n] - ex_pbe [n] =", ec_pbe_at_n_plus_delta_n - ec_pbe_at_n
   print*, "n variation integral              :", int_vc_pbe_at_n 

   print*, '___________PBEUEG___________'
   print*, "ex_pbeUEG [n, grad_n]                                                  =", ex_pbeUEG_at_n
   print*, "ex_pbeUEG [n + delta_n, grad_n + delta_grad_n]                         =", ex_pbeUEG_at_n_plus_delta_n
   print*, "ex_pbeUEG [n + delta_n, grad_n + delta_grad_n] - ex_pbeUEG [n, grad_n] =", ex_pbeUEG_at_n_plus_delta_n - ex_pbeUEG_at_n
   print*, "n and grad_n variations integral                                       :", int_vx_pbeUEG_at_n 
   print*, "ec_pbeUEG [n, grad_n]                                                  =", ec_pbeUEG_at_n
   print*, "ec_pbeUEG [n + delta_n, grad_n + delta_grad_n]                         =", ec_pbeUEG_at_n_plus_delta_n
   print*, "ec_pbeUEG [n + delta_n, grad_n + delta_grad_n] - ex_pbeUEG [n, grad_n] =", ec_pbeUEG_at_n_plus_delta_n - ec_pbeUEG_at_n
   print*, "n and grad_n variations integral                                       :", int_vc_pbeUEG_at_n 

   print*, '___________PBEn2___________'
   print*, "ec_pben2 [n,n2]                                         =", ec_pben2_at_n_n2
   print*, "ec_pben2 [n + delta_n, n2 + delta_n2]                   =", ec_pben2_at_n_plus_delta_n_n2_plus_delta_n2
   print*, "ec_pben2 [n + delta_n, n2]                              =", ec_pben2_at_n_plus_delta_n
   print*, "ec_pben2 [n, n2 + delta_n2]                             =", ec_pben2_at_n2_plus_delta_n2
   print*, "ec_pben2 [n + delta_n, n2 + delta_n2] - ec_pben2 [n,n2] =", ec_pben2_at_n_plus_delta_n_n2_plus_delta_n2 - ec_pben2_at_n_n2
   print*, "n and n2 variation integral                             :", int_vc_pben2_one_e_at_n_n2 + int_vc_pben2_two_e_at_n_n2
   print*, "ec_pben2 [n + delta_n, n2] - ec_pben2 [n,n2]            =", ec_pben2_at_n_plus_delta_n - ec_pben2_at_n_n2
   print*, " only n variation integral                              :", int_vc_pben2_one_e_at_n_n2
   print*, "ec_pben2 [n, n2 + delta_n2] - ec_pben2 [n,n2]           =", ec_pben2_at_n2_plus_delta_n2 - ec_pben2_at_n_n2
   print*, " only n2 variation integral                             :", int_vc_pben2_two_e_at_n_n2 
end program

