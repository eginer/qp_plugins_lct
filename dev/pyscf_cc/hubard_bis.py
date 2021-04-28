#!/usr/bin/env python

# 10-site Hubbard model at half-filling with U/t = 4
import numpy as np
from pyscf import gto, scf, ao2mo, cc
mol = gto.Mole(verbose=4)
mol.nelectron = n = 10

# t,u definition 
t, u = 1., 4.


mf = scf.RHF(mol)
h1 = np.zeros((n,n))

for i in range(n-1):
  h1[i,i+1] = h1[i+1,i] = t

mf.get_hcore = lambda *args: h1
mf.get_ovlp = lambda *args: np.eye(n)
mf._eri = np.zeros((n,n,n,n))

for i in range(n):
  mf._eri[i,i,i,i] = u

# 2e Hamiltonian in 4-fold symmetry
mf._eri = ao2mo.restore(4, mf._eri, n)
mf.run()
cc.CCSD(mf).run()
