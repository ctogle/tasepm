# imports from modular_core
import modular_core.fundamental as lfu
import modular_core.modules.simulationmodule as smd
import modular_core.parameterspace.parameterspaces as lpsp
import modular_core.cython.writer as cwr
import modular_core.io.liboutput as lo

# imports from standard python libraries
import pdb,os,sys,time,numpy
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

# adding this class allows this module to specify more information
# for the core to use
class simulation_module(smd.simulation_module):

    # the 'simulation' attribute of this class will be a function
    # which takes one argument and returns the data from one realization
    def __init__(self,*args,**kwargs):
        self.simulation = simulate
        smd.simulation_module.__init__(self,*args,**kwargs)

    # this is called before subprocesses are created which inherit the state
    # of this module - used to get info in multiprocessing runs...
    # no information set here should depend on parameters being scanned
    def _set_parameters_prepoolinit(self):
        tcount = len(self.parent.simulation_plan.plot_targets)
        ccount = self.parent.simulation_plan._capture_count()
        self.dshape = (tcount,ccount)
        self.captureincrement = self.parent.simulation_plan._capture_increment()

    # whenever the pspace location changes this function is called
    # this happens once even when not mapping parameter space
    # setting attribute 'sim_args';it can be passed to simulation
    def _set_parameters(self):
        plottargets = self.parent.run_params['plot_targets']
        ccrit = self.parent.run_params['capture_criteria'][0]
        ecrit = self.parent.run_params['end_criteria'][0]
        self.sim_args = [self.dshape,self.captureincrement,ccrit,ecrit,plottargets]

# this function which receives the attribute 'sim_args' must return
# the data of one realization of the simulation
def simulate(args):
    dshape,cincr,ccrit,ecrit,ptargs = args
    print 'ptargets',ptargs
    data = numpy.zeros(dshape,dtype = numpy.float)
    print 'resultant data has shape:',dshape
    


    # this is where the simulation happens
    #del_t = -1.0*log(<float>random.random())*tpinv



    return data



