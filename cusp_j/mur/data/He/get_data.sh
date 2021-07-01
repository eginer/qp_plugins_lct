method=grad_n 
for mu in 0.2 0.3 0.4 0.5 0.6 0.7 
 do 
 grep "$mu -" data_$method > data_${method}_$mu
done
