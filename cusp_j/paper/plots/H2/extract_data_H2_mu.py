#!/usr/bin/env python3

import sys
import os
def float_round_list(list_1):
  list_2 = []
  list_tmp = [float(element) for element in list_1]
  list_2   = [round(element*1000.,2) for element in list_tmp ]
  return list_2

def print_in_file(list_name,filename,eref):
  ns = 0
  if (os.path.isfile(filename)): 
   os.remove(filename)
  with open(filename, 'a') as the_file:
   new_list = []
   for sys in system: 
    new_list = float_round_list(list_name[ns])
    the_file.write(sys+'     '+str(eref[ns])+ ' ' + str(new_list)+ '\n')
    ns +=1

def print_e_in_file(list_name,filename,eref):
  ns = 0
  if (os.path.isfile(filename)): 
   os.remove(filename)
  with open(filename, 'a') as the_file:
   for sys in system: 
    list_tmp = [float(element) for element in list_name[ns]]
    the_file.write(sys+'     '+str(eref[ns])+ ' ' + str(list_tmp)+ '\n')
    ns +=1


def subst_in_file(filename):
  comand='sed -i "s/\[/  /g" '+ filename
  os.system(comand)
  comand='sed -i "s/\]/  /g" '+ filename
  os.system(comand)
  comand='sed -i "s/,/  /g" '+ filename
  os.system(comand)


mu=sys.argv[1]

basis = []

basis.append(2)
basis.append(3)
basis.append(4)
basis.append(5)

eref = []
system = []
fileH2="exact-H2_bis"
nb=0
with open(fileH2, "r") as fp:
 for line in fp:
   a=line.split()
   if(a[0].find("#") == 0):  
    continue
   else:
    eref.append(float(a[1]))
    system.append(a[0])

e_mu_full = []
de_mu_full = []
ns=0
for sys in system: 
 print(sys)
 eref_tmp=eref[ns]

 filemu='H2/'+sys+'/data_mu_'+mu
 
 nb=0
 e_mu = []
 de_mu = []
 with open(filemu, "r") as fp:
   for line in fp:
     a=line.split()
     if(a[0].find("#") == 0):  
      continue
     else:
      #MU LDA 
      emu=round(float(a[2]),8)
      print(emu)
      demu=round((emu - eref_tmp),8)
      e_mu.append(str(emu)) 
      de_mu.append(str(demu)) 
     nb+= 1
 de_mu_full.append(de_mu)
 e_mu_full.append(e_mu)
 ns += 1

###################### fci 
filename="data_H2_de_mu_"+mu
print_in_file(de_mu_full,filename,eref)
subst_in_file(filename)
filename="data_H2_e_mu_"+mu
print_e_in_file(e_mu_full,filename,eref)
subst_in_file(filename)
