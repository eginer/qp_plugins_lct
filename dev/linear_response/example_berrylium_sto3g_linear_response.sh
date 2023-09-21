#/bin/bash                                                                                              
QP_ROOT=~/programs/qp2/
source ${QP_ROOT}/quantum_package.rc

system=Be
INPUT=${system}.xyz

multiplicity=1
basis="sto-3g"

 mkdir ${system}
 cd ${system}
 cp ../${system}.xyz .

  for b in ${basis}
  do
   mkdir $b
   cd ${b}
   cp ../${system}.xyz .

   qp create_ezfio -b ${b} -m $multiplicity -o ${system}_${b}.ezfio ${INPUT}

   qp run scf | tee scf_${system}_${b}.out
   qp run cisd | tee cisd_${system}_${b}.out
   qp run save_natorb | tee save_natorb_${system}_${b}.out
   qp run scf | tee scf2_${system}_${b}.out

   qp_set_frozen_core -l ${system}_${b}.ezfio

   qp_sc_basis_corr -m hf -f lda_ueg -p 0.0001 -n 10000000 ${system}_${b}.ezfio -i 3 | tee sc_bc_${system}_${b}.out 
   
   # Save the one electron density matrix from the scf basis-set correction
   qp run save_one_e_dm | tee save_one_e_dm.out

   qp reset -d

   # Generate the needed set of determinants for the linear response calculation
   qp set determinants n_states 16 
   # This line permits to consider all possible determinants:
   qp set determinants save_threshold 0.0 
   qp run cisd | tee cisd_post_bsc.out
   qp run fci | tee fci_post_bsc.out

   # Save only the ground state coefficient
   qp edit -s 1
   qp run diagonalize_h | tee diagonalize_h.out

   qp set density_for_dft density_for_dft input_density_ao
   qp run linear_response_K_equal_0 | tee linear_response_K0.out
   qp run linear_response_Tamm_Dancoff | tee linear_response_TD.out
   qp run linear_response | tee linear_response.out

   cd ..
  done
  cd ../

