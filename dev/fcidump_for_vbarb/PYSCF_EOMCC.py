#!/usr/bin/env python
#
# Author: Qiming Sun <osirpt.sun@gmail.com>
#

'''
Writing FCIDUMP file for given integrals or SCF orbitals
'''

import os
import sys
from functools import reduce
import numpy
import pyscf


from pyscf import gto, scf, ao2mo, cc
from pyscf.cc import ccsd_t
from pyscf import symm
from pyscf.tools import fcidump

#basis_sets = ['aug-cc-pvdz', 'aug-cc-pvtz', 'aug-cc-pvqz', 'aug-cc-pv5z']
#basis =str(sys.argv[1])
#### usage : 
## PYSCF_EOMCC.py xyz_file basis n_roots n_frozen 


xyz_file=str(sys.argv[1]) 
basis=str(sys.argv[2])
n_roots=int(sys.argv[3])
n_frozen=int(sys.argv[4])

mol = gto.Mole()
#mol.atom = 'C  0.00000000 0.00000000 -1.24942055; O  0.00000000 0.00000000 0.89266692'
mol.atom = xyz_file
mol.unit = 'B'
mol.basis = basis
mol.build()

mf = scf.RHF(mol)
mf.verbose = 7
mf.scf()

mf.mol.verbose = 4

mycc = cc.CCSD(mf,frozen=n_frozen).run()

#e_ee, c_ee = mycc.eeccsd(nroots=10)
eS, c_S = mycc.eomee_ccsd_singlet(nroots=n_roots)
iroot=0
for e in eS:
  iroot+= 1
  print('E_S '+str(iroot)+' =  '+str(e)+'  eV = '+str(e*27.2114))

eT, c_T = mycc.eomee_ccsd_triplet(nroots=n_roots)
for e in eT:
  iroot+= 1
  print('E_T '+str(iroot)+' =  '+str(e)+'  eV = '+str(e*27.2114))
