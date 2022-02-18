# source /home_lct/eginer/qp2/quantum_package.rc
 atom=$1 # ATOM 
 basis=$2 # BASE 
 echo $atom > ${atom}.xyz
 rm -rf $atom.$basis
 qp create_ezfio -b $basis ${atom}.xyz -m 1 -o $atom.$basis
 qp run scf
 # prints out typical values of mu, mu_lda and mu_rsc_lda seem good 
 qp run print_mu_av_tc | tee ${EZFIO_FILE}/mu_av.out
 mu_rsc_lda=` grep "average_mu_rs_c_lda" ${EZFIO_FILE}/mu_av.out | cut -d "=" -f 2`
 qp set ao_two_e_erf_ints mu_erf $mu_rsc_lda 

 qp set tc_h_clean three_body_h_tc true       # THREE BODY TERMS 
 qp set tc_h_clean pure_three_body_h_tc False # NO PURE TRIPLE EXCITATIONS 
 qp set tc_h_clean full_tc_h_solver False     # DO NOT DIAGONALIZE THE FULL TC MATRIX, USE DAVIDSON 
 qp set tc_h_clean comp_left_eigv True        # COMPUTE ALSO THE LEFT EIGENVECTOR 

 qp set perturbation pt2_max 0.00001          # WHERE TO STOP FOR PT2 

 qp set tc_h_clean double_normal_ord False    # NORMAL ORDER  

############### SELECTION WITH TC H 
 qp set fci_tc cipsi_tc h_tc                  # NON HERMITIAN SELECTION
 qp run fci_tc | tee ${EZFIO_FILE}.fci_tc_tc_h.out

 qp reset -d # RESET THE WAVE FUNCTION 
 qp set determinants weight_selection 6       # TO SELECT WITH THE COEFFICIENTS  
 qp run fci_tc | tee ${EZFIO_FILE}.fci_tc_tc_h_select_coef.out


############### SELECTION WITH USUAL H 
 qp set fci_tc cipsi_tc reg_h #### USUAL H SELECTION
 qp reset -d
 qp run fci_tc | tee ${EZFIO_FILE}.fci_tc_reg_h.out

############### SELECTION WITH USUAL H_TC + H_TC^dagger 
 qp set fci_tc cipsi_tc sym_h_tc #### USUAL H SELECTION
 qp reset -d
 qp run fci_tc | tee ${EZFIO_FILE}.fci_tc_sym_h_tc.out

# rm -rf  ${EZFIO_FILE}_tc_scf
# cp -r ${EZFIO_FILE} ${EZFIO_FILE}_tc_scf
# qp set_file ${EZFIO_FILE}_tc_scf
# qp run tc_scf | tee ${EZFIO_FILE}.out 
# 
# qp reset -d 
# qp run fci_tc | tee ${EZFIO_FILE}.fci_tc_non_h.out 

