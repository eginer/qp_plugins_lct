use bitmasks
BEGIN_SHELL [ /usr/bin/env python2 ]
from generate_h_apply import *

s = H_apply("just_mono",do_double_exc=False)
s.set_selection_pt2("epstein_nesbet_2x2")
print s

s = H_apply("just_mono_hcore",do_double_exc=False)
s.set_selection_pt2("hcore")
print s

H = H_apply("just_mono_no_selection",do_double_exc=False)
print H



END_SHELL

