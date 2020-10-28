source ~/programs/qp2/quantum_package.rc

basis="d t q"
basis="5"
declare -i b=2
rm data_mu_lda 
for i in $basis
 do 
 prefix=aug-cc-pv
 dir=${prefix}${i}z
 cd $dir
 ezfio=H_$dir
 qp set_file $ezfio
 qp run mu_lda_average | tee ${ezfio}.mu_av.out 
 mu_av=`grep "average_mu_lda ="  ${ezfio}.mu_av.out | cut -d "=" -f 2`
 mu_av=${mu_av:0:10}
 mu_av="$(echo -e "${mu_av}" | tr -d '[:space:]')"
 echo ${mu_av}
 mu_av2=`grep "average_mu_lda " ${ezfio}.mu_av.out  | grep "1/2"  | cut -d "=" -f 2`
 mu_av2=${mu_av2:0:10}
# mu_av2=${mu_av2##*( )}
 mu_av2="$(echo -e "${mu_av2}" | tr -d '[:space:]')"
 echo ${mu_av2}
 ezfio2=${ezfio}_mu_lda_${mu_av}
 echo $ezfio2
 cp -r $ezfio $ezfio2
 qp set_file $ezfio2 
 qp set ao_two_e_erf_ints mu_erf $mu_av
 qp run transcorr_general | tee ${ezfio2}.htilde.out
 E_mu_av=`grep "E0_tilde =" ${ezfio2}.htilde.out | cut -d "=" -f 2`

 ezfio2=${ezfio}_mu_lda_${mu_av2}
 cp -r $ezfio $ezfio2
 qp set_file $ezfio2 
 qp set ao_two_e_erf_ints mu_erf $mu_av2
 qp run transcorr_general | tee ${ezfio2}.htilde.out
 E_mu_av2=`grep "E0_tilde =" ${ezfio2}.htilde.out | cut -d "=" -f 2`
 echo  $b $mu_av $mu_av2 $E_mu_av $E_mu_av2 >> ../data_mu_lda 
 b=${b}+1
 cd ../
done 
