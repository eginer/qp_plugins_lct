#!/usr/bin/env python
#
# Author: Qiming Sun <osirpt.sun@gmail.com>
#

'''
Writing FCIDUMP file for given integrals or SCF orbitals
'''

import os
from functools import reduce
import numpy


from pyscf import gto, scf, ao2mo, cc
from pyscf.cc import ccsd_t
from pyscf import symm
from pyscf.tools import fcidump

#mol = gto.M(atom='Li 0. 0. 0.', basis='cc-pcvtz')
mol = gto.Mole()
mol.atom = 'Li 0. 0. 0.'
mol.basis = 'cc-pcvtz'
mol.charge = 0
mol.spin = 1  
mol.build()
FCIDUMP='B.ezfio.FCIDUMP'


#
# Hamiltonians of FCIDUMP file can be load
#
ctx = fcidump.read(FCIDUMP)

#
# Construct an SCF object using the quantities defined in FCIDUMP
# (pyscf-1.7.4 or newer)
#
mf = fcidump.to_scf(FCIDUMP, molpro_orbsym=True)
mf.mol.verbose = 4
mf.run()
#mf.MP2().run()

mycc = cc.CCSD(mf).run()
et=mycc.ccsd_t()
