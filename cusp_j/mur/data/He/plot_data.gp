set terminal pdf enhanced
set output "conv_He_lda.pdf"
set ylabel "Tota energy" font ",12"
#set xlabel "Cardinal number" font ",12"
set key font "10,12"
set key spacing 1.5
set key right top
set xtics ("aug-cc-pVDZ" 2, "aug-cc-pVTZ" 3, "aug-cc-pVQZ" 4, "aug-cc-pV5Z" 5 )

plot "data_fci" w lp title "FCI","data_lda_0.7" w lp title "LDA-0.7", "data_lda_0.6" w lp title "LDA-0.6", "data_lda_0.5" w lp title "LDA-0.5", "data_lda_0.4" w lp title "LDA-0.4", "data_lda_0.3" w lp title "LDA-0.3", "data_lda_0.2" w lp title "LDA-0.2", -2.903724 lw 2 title "Exact energy"

set terminal pdf enhanced
set output "conv_He_rsc.pdf"
set ylabel "Tota energy" font ",12"
#set xlabel "Cardinal number" font ",12"
set key font "10,12"
set key spacing 1.5
set key right top
set xtics ("aug-cc-pVDZ" 2, "aug-cc-pVTZ" 3, "aug-cc-pVQZ" 4, "aug-cc-pV5Z" 5 )

plot  "data_fci" w lp title "FCI","data_rsc_0.7" w lp title "rsc-0.7", "data_rsc_0.6" w lp title "rsc-0.6", "data_rsc_0.5" w lp title "rsc-0.5", "data_rsc_0.4" w lp title "rsc-0.4", "data_rsc_0.3" w lp title "rsc-0.3", "data_rsc_0.2" w lp title "rsc-0.2", -2.903724 lw 2 title "Exact energy"

set terminal pdf enhanced
set output "conv_He_rsc_lda.pdf"
set ylabel "Tota energy" font ",12"
#set xlabel "Cardinal number" font ",12"
set key font "10,12"
set key spacing 1.5
set key right top
set xtics ("aug-cc-pVDZ" 2, "aug-cc-pVTZ" 3, "aug-cc-pVQZ" 4, "aug-cc-pV5Z" 5 )

plot "data_fci" w lp title "FCI","data_rsc_lda_0.7" w lp title "rsc-lda-0.7", "data_rsc_lda_0.6" w lp title "rsc-lda-0.6", "data_rsc_lda_0.5" w lp title "rsc-lda-0.5", "data_rsc_lda_0.4" w lp title "rsc-lda-0.4", "data_rsc_lda_0.3" w lp title "rsc-lda-0.3", "data_rsc_lda_0.2" w lp title "rsc-lda-0.2", -2.903724 lw 2 title "Exact energy"

set terminal pdf enhanced
set output "conv_He_grad_n.pdf"
set ylabel "Tota energy" font ",12"
#set xlabel "Cardinal number" font ",12"
set key font "10,12"
set key spacing 1.5
set key right bottom
set xtics ("aug-cc-pVDZ" 2, "aug-cc-pVTZ" 3, "aug-cc-pVQZ" 4, "aug-cc-pV5Z" 5 )

plot "data_fci" w lp title "FCI","data_grad_n_0.7" w lp title "grad n-0.7", "data_grad_n_0.6" w lp title "grad n-0.6", "data_grad_n_0.5" w lp title "grad n-0.5", "data_grad_n_0.4" w lp title "grad n-0.4", "data_grad_n_0.3" w lp title "grad n-0.3", "data_grad_n_0.2" w lp title "grad n-0.2", -2.903724 lw 2 title "Exact energy"
