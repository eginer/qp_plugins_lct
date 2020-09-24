source ~/programs/qp2/quantum_package.rc

#basis="d"
basis="d t q"
#mu="0.005 0.05 0.1 0.2 0.3 0.325 0.35 0.375 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.2 1.4 1.6 1.8 2.0 2.5 3.0 3.5 4.0 5.0 6.0 7.0 8.0 9.0 10.0"
mu="0.005 0.05 0.1 0.2 0.3 0.325 0.35 0.375 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.2 1.4 1.6 1.8 2.0 2.5 3.0 3.5 4.0"
#mu=" 0.5"
for i in $basis
do
 rm data_$i
 prefix=aug-cc-pv
 dir=${prefix}${i}z
 mkdir $dir
 cd $dir
 echo "He" > He.xyz 
 ezfio=He_$dir
 rm -rf $ezfio
 qp create_ezfio -b ${dir} He.xyz -o $ezfio 
 qp run scf 
 qp run fci | tee ${ezfio}.fci.out 
 qp set becke_numerical_grid my_grid_becke true
 qp set becke_numerical_grid my_n_pt_r_grid 30
 qp set becke_numerical_grid my_n_pt_a_grid 74
 for j in $mu
 do
  ezfio_mu=${ezfio}_$j
  cp -r $ezfio $ezfio_mu
  qp set_file $ezfio_mu 
  qp set ao_two_e_erf_ints mu_erf $j 
  qp run transcorr_h | tee ${ezfio_mu}.htilde.out 
  file=${ezfio_mu}.htilde.out
  E=`grep "E0_tilde =" $file | cut -d "=" -f 2`
  echo $j $E >> data_$i
 done
 cd ../
done
