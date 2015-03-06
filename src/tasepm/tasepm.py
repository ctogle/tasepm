# imports from modular_core
import modular_core.fundamental as lfu
import modular_core.simulationmodule as smd
import modular_core.parameterspaces as lpsp
import modular_core.cython.writer as cwr
import modular_core.io.liboutput as lo

# imports from standard python libraries
import pdb,os,sys,time
from cStringIO import StringIO

# this is true if this file is imported as a module
if __name__ == 'tasepm.tasepm':
    # provide access to modular_core.gui package as needed...
    lfu.check_gui_pack()
    lgm = lfu.gui_pack.lgm
    lgb = lfu.gui_pack.lgb
    lgd = lfu.gui_pack.lgd
# this is true if you run this file like a script
if __name__ == '__main__':print 'tasepm module'

###############################################################################

# modular_core requires this attribute
module_name = 'tasepm'






