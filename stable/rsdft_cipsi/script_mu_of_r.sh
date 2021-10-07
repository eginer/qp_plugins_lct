source ~/programs/new_qp2/qp2/quantum_package.rc 

basis=aug-cc-pvdz 
mu="rsc"
ezfio=He_${basis}_${mu}
rm -rf $ezfio
rm data_$mu
qp create_ezfio -b $basis He.xyz -o $ezfio
qp set_file $ezfio 
qp run scf 
qp run fci | tee ${ezfio}.fci.out 
qp set dft_keywords mu_dft_type $mu 
qp run save_mu_of_r_ints 
cd $ezfio/work
for i in *erf*   
do  
 ln -s $i ${i/_erf/} 
done
cd ../../

for j in `seq 1 10`
#for j in `seq 1 1`
do 
qp run write_effective_rsdft_hamiltonian_mu_of_r | tee ${ezfio}.$j.${mu}.out 
e=`grep "TOTAL ENERGY        =" He.ezfio.$j.out | cut -d "=" -f 2`
echo $j $e >> data_$mu
qp run diagonalize_h | tee ${ezfio}.diag.$j.${mu}.out 
done
qp run print_ecmd_components | tee ${ezfio}.ecmd.${mu}.out 
