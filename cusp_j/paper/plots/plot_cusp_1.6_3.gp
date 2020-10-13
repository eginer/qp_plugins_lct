
set terminal pdf enhanced
set output "He_mu_1.6_cusp_avdz_avtz_3.pdf"
set ylabel "Wave functions, {/ Symbol m} = 1.6"
set xlabel "{/ Symbol q}_{12}"
set key bottom right
set key spacing 1.3

plot "aug-cc-pvdz/He_aug-cc-pvdz_1.6.cusp_wf_3" u 1:3 w l dt 3 lw 3 title "{/ Symbol Y}_0^{ex}",      "aug-cc-pvdz/He_aug-cc-pvdz_1.6.cusp_wf_3" u 1:7 w l title "AVDZ: {/ Symbol Y}_0^{B,{/ Symbol m}}",      "aug-cc-pvtz/He_aug-cc-pvtz_1.6.cusp_wf_3" u 1:7 w l title "AVTZ: {/ Symbol Y}_0^{B,{/ Symbol m}}",      "aug-cc-pvdz/He_aug-cc-pvdz_1.6.cusp_wf_3" u 1:4 w l title "AVDZ: FCI",       "aug-cc-pvtz/He_aug-cc-pvtz_1.6.cusp_wf_3" u 1:4 w l title "AVTZ: FCI",       "aug-cc-pvdz/He_aug-cc-pvdz_1.6.cusp_wf_3" u 1:5 w l title "AVDZ: {/ Symbol F}_{/ Symbol m}^{B}",      "aug-cc-pvtz/He_aug-cc-pvtz_1.6.cusp_wf_3" u 1:5 w l title "AVTZ: {/ Symbol F}_{/ Symbol m}^{B}"


set terminal pdf enhanced
set output "He_mu_1.6_cusp_avqz_av5z_3.pdf"
set ylabel "Wave functions, {/ Symbol m} = 1.6"
set xlabel "{/ Symbol q}_{12}"
set key bottom right
set key spacing 1.3

plot "aug-cc-pvqz/He_aug-cc-pvqz_1.6.cusp_wf_3" u 1:3 w l dt 3 lw 3 title "{/ Symbol Y}_0^{ex}",      "aug-cc-pvqz/He_aug-cc-pvqz_1.6.cusp_wf_3" u 1:7 w l title "AVQZ: {/ Symbol Y}_0^{B,{/ Symbol m}}",      "aug-cc-pv5z/He_aug-cc-pv5z_1.6.cusp_wf_3" u 1:7 w l title "AV5Z: {/ Symbol Y}_0^{B,{/ Symbol m}}",      "aug-cc-pvqz/He_aug-cc-pvqz_1.6.cusp_wf_3" u 1:4 w l title "AVQZ: FCI",       "aug-cc-pv5z/He_aug-cc-pv5z_1.6.cusp_wf_3" u 1:4 w l title "AV5Z: FCI",       "aug-cc-pvqz/He_aug-cc-pvqz_1.6.cusp_wf_3" u 1:5 w l title "AVQZ: {/ Symbol F}_{/ Symbol m}^{B}",      "aug-cc-pv5z/He_aug-cc-pv5z_1.6.cusp_wf_3" u 1:5 w l title "AV5Z: {/ Symbol F}_{/ Symbol m}^{B}"
 
