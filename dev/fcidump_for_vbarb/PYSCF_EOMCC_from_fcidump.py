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

#### usage : 
## PYSCF_EOMCC_from_fcidump.py xyz_file basis FCIDUMP n_roots 
xyz_file=str(sys.argv[1])
basis=str(sys.argv[2])
FCIDUMP=str(sys.argv[3])
n_roots=int(sys.argv[4])


mol = gto.Mole()
#mol.atom = 'C  0.00000000 0.00000000 -1.24942055; O  0.00000000 0.00000000 0.89266692'
mol.atom = xyz_file
mol.unit = 'B'
mol.basis = basis
mol.build()

print('reading the following FCIDUMP '+ FCIDUMP)

ctx = fcidump.read(FCIDUMP)
mf = fcidump.to_scf(FCIDUMP, molpro_orbsym=True)
mf.run()

mf.mol.verbose = 4

mycc = cc.CCSD(mf,frozen=0).run()

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
