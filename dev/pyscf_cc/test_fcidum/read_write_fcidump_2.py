#!/usr/bin/env python
#
# Author: Qiming Sun <osirpt.sun@gmail.com>
#

'''
Writing FCIDUMP file for given integrals or SCF orbitals
'''

from functools import reduce
import numpy
from pyscf import gto, scf, ao2mo, cc
from pyscf.cc import ccsd_t
from pyscf import symm
from pyscf.tools import fcidump

mol = gto.M(atom='H 0 0 0; H 0 0 0.7', basis='cc-pvtz')

#
# Hamiltonians of FCIDUMP file can be load
#
ctx = fcidump.read('H2_0.7.ezfio.FCIDUMP')

#
# Construct an SCF object using the quantities defined in FCIDUMP
# (pyscf-1.7.4 or newer)
#
mf = fcidump.to_scf('H2_0.7.ezfio.FCIDUMP', molpro_orbsym=True)
mf.mol.verbose = 4
mf.run()
mf.MP2().run()

mycc = cc.CCSD(mf)
mycc.kernel()
e_ee, c_ee = mycc.eeccsd(nroots=5)
