===============
mu_of_r_tc_ints
===============

gaussian integrals       : ao_one_ints_slater/gauss_ints.irp.f
erf(mu r12)/r12          : deriv_r12_ints/deriv_one_e_ints.irp.f
                         : deriv_r12_ints/semi_num_ints_ao.irp.f provider: v_ij_erf_rk
erf(mu r12)/r12 * x/y/z  : deriv_r12_ints/semi_num_ints_ao.irp.f function: NAI_pol_x_mult_erf_ao

needed     : x/y/z * e^{-alpha * |r-R|^2} * erf(mu r12)/r12 
