#!/bin/bash                                                                                                                                  
QP_ROOT=~/programs/qp2/
source ${QP_ROOT}/quantum_package.rc

#### Assumes that an "xyz" file exists for H2 
file_xyz=He.xyz
basis=cc-pvdz
ezfio=He

rm -rf $ezfio

#start with a RHF calculation
qp create_ezfio -b $basis -o $ezfio $file_xyz
qp run scf | tee ${ezfio}.scf.out
qp run fci | tee ${ezfio}.fci.out
qp set mu_of_r mu_of_r_potential cas_ful
qp run basis_correction | tee ${ezfio}.DFT.out

#define the main options : mu, functional, number of determinants, max pt2 and reading or not two-elec integrals (to save time)

qp run write_eff_mu_of_r_ints | tee ${ezfio}.eff_ints_1_out
qp run diagonalize_h | tee ${ezfio}.diag_h_1_out

qp run write_eff_mu_of_r_ints | tee ${ezfio}.eff_ints_2_out
qp run diagonalize_h | tee ${ezfio}.diag_h_2_out

qp run write_eff_mu_of_r_ints | tee ${ezfio}.eff_ints_3_out
qp run diagonalize_h | tee ${ezfio}.diag_h_3_out                      
