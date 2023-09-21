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
basis =str(sys.argv[1])


mol = gto.Mole()
mol.atom = 'O  0.00000000 0.00000000 -0.13209669; H  0.00000000 1.43152878 0.97970006; H  0.00000000 -1.43152878 0.97970006' #Angstrom
mol.unit = 'B'
mol.basis = basis
mol.build()

FCIDUMP=str(sys.argv[2])
print('reading the following FCIDUMP '+ FCIDUMP)

ctx = fcidump.read(FCIDUMP)
mf = fcidump.to_scf(FCIDUMP, molpro_orbsym=True)
mf.run()

mf.mol.verbose = 4

mycc = cc.CCSD(mf,frozen=0).run()

e_ee, c_ee = mycc.eeccsd(nroots=10)
