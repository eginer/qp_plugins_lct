
set terminal pdf enhanced
set output "He_E_conv_basis_small_mu.pdf"
set ylabel "Energy (a.u.)" font ",15" 
set key font "10,15"
set key spacing 1.5
set xlabel "1/X"
set key left top
set xtics ("AV2Z" 0.5, "AV3Z" 0.333, "AV4Z" 0.25, "AV5Z" 0.20 , "AV6Z" 0.16666)



plot  "data_mu_0.2" u (1./$1):2 w lp lw 2 title  "{/Symbol m} = 0.2", "data_mu_0.3" u (1./$1):2 w lp lw 2  title "{/Symbol m} = 0.3", "data_mu_0.5" u (1./$1):2 w lp lw 2  title "{/Symbol m} = 0.5",  "data_mu_0.7" u (1./$1):2 w lp lw 2  title "{/Symbol m} = 0.7", "data_ten_no" u (1./$1):3 w lp lw 2 title "FROGG",  "data_fci" u (1./$1):2 w lp lw 2 title "FCI",  -2.903724 title "Exact" dt 3 lw 4

set output "He_E_conv_basis_large_mu.pdf"
set ylabel "Energy (a.u.)" font ",15" 
set xlabel "1/X"
set key left top
plot   "data_mu_0.7" u (1./$1):2 w lp lt 4 lw 2  title "{/Symbol m} = 0.7",  "data_ten_no" u (1./$1):3 w lp lt 5 lw 2 title "FROGG", "data_mu_2.0" u (1./$1):2 w lp lw 2  title "{/Symbol m} = 2.0", "data_mu_3.0" u (1./$1):2 w lp lw 2 lt 8 title "{/Symbol m} = 3.0", "data_fci" u (1./$1):2 w lp lw 2 lt 6 title "FCI",  -2.903724 title "Exact" dt 3 lw 4 lt 7

#set output "He_conv_basis_mu_lda.pdf"
#set ylabel "Total energy"
#set xlabel "1/X"
#set key left top
#plot "data_mu_lda" u (1./$1):4 w lp lw 2  title "{/Symbol m}_{LDA}",  "" u (1./$1):5 w lp lw 2  title "{/Symbol m}_{LDA} * 0.5", "data_fci" u (1./$1):2 w lp lw 2 title "FCI",  -2.903724 title "Exact" dt 3
