# imports from modular_core
import modular_core.fundamental as lfu
import modular_core.modules.simulationmodule as smd
import modular_core.parameterspace.parameterspaces as lpsp
import modular_core.cython.writer as cwr
import modular_core.io.liboutput as lo

# imports from standard python libraries
import pdb,os,sys,time,numpy,random
from cStringIO import StringIO
from math import log

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
        #tcount = len(self.parent.simulation_plan.plot_targets)
        tcount = len(self.parent.run_params['plot_targets'])
        ccount = self.parent.simulation_plan._capture_count()
        self.dshape = (tcount,ccount)
        self.capture = self.parent.simulation_plan._capture_increment()
        self.endtime = self.parent.simulation_plan._max_time()
        self.targets = self.parent.run_params['plot_targets']

    def _reset_parameters(self):
        def_targeted = ['time']
        smd.simulation_module._reset_parameters(self,def_targeted)

    # whenever the pspace location changes this function is called
    # this happens once even when not mapping parameter space
    # setting attribute 'sim_args';it can be passed to simulation
    def _set_parameters(self):
        self.sim_args = [self.dshape,self.capture,self.endtime,self.targets]

    def _panel_templates(self,*args,**kwargs):
        window = args[0]
        ensem = args[1]
        self._gui_memory()
        plot_target_labels = ['time','initflux','gateflux','termflux']
        panel_template_lookup =\
            smd.simulation_module._panel_templates(
                self,window,ensem,plot_target_labels)
        return panel_template_lookup

# from a position in the lattice, are all sites up to distance away vacant?
def scan_ahead(latt,position,distance,count):
    #if position + distance > latt.shape[0]:
    #    distance = latt.shape[0] - position
    if count - position < distance:
        distance = count - position

    for x in range(1,distance):  
        try:
            if latt[position+x] > 0:
                return 0
        except:
            pdb.set_trace()
    return 1

# set the propensities for movement based on current state
def set_propensities(props,latt,initrate,speeds,count,awidth):
    # props[0] is for initializing a new agent at site 0
    if latt[0] == 0 and scan_ahead(latt,0,awidth,count):props[0] = initrate
    else:props[0] = 0.0
    totalprop = props[0]

    # props[k+1] corresponds to the propensity to move of an
    # agent located at lattice site k, that is, an agent whose
    # position spans sites k to k+awidth
    for x in range(count):
        if latt[x] > 0:
            if scan_ahead(latt,x,awidth,count):
                thismoveprop = 1.0 * speeds[x]
                thisfallprop = 1.0 - thismoveprop
            else:
                thismoveprop = 0.0
                thisfallprop = 1.0
            props[x+1] = totalprop + thismoveprop
            #props[x+count+1] = totalprop + thisfallprop
            props[x+count+1] = thisfallprop
            totalprop += thismoveprop
        else:
            props[x+1] = 0.0
            props[x+count+1] = 0.0

    for x in range(count):
        if latt[x] > 0:
            thisfallprop = props[x+count+1]
            props[x+count+1] += totalprop
            totalprop += thisfallprop

    return totalprop
    
# given the propensities for motion associated with each lattice site
# and a propensity for initiating a new agent at site 0, pick a
# site, which will by construction correspond to an agent
def pick_action(lattice,props,totalp,ccount):
    rand = random.random()
    for adx in range(ccount+1):
        if rand < props[adx]/totalp:
            return adx

    pdb.set_trace()
    return -1

pmin = 10**-30  # threshold to consider a propensity equal to zero
# this function which receives the attribute 'sim_args' must return
# the data of one realization of the simulation
def simulate(args):
    dshape,capt,end,ptargs = args
    data = numpy.zeros(dshape,dtype = numpy.float)

    alltargets = ['time','iteration','initflux','gateflux','termflux','stalled']
    state = numpy.zeros((len(alltargets)),dtype = numpy.float)
    ptdexes = [alltargets.index(pt) for pt in ptargs]
    ptcount = len(ptdexes)

    time = 0.0
    iteration = 0
    completed = 0
    stalled = 0

    capture = 0
    captures = dshape[1]
    tcount = dshape[0]

    realtime = 0.0
    lasttime = 0.0
    fixed_del_t = capt
    
    rwidth = 12
    ccount = 300
    ocount = 2
    irate = 10.0
    speeds = 1.0*numpy.ones((ccount),dtype = numpy.float)
    speeds[int(0.4*ccount):int(0.6*ccount)] = 0.75

    lattice = numpy.zeros((ccount),dtype = numpy.float)
    fluxes = numpy.zeros((ccount),dtype = numpy.float)
    props = numpy.zeros((ocount*ccount+1),dtype = numpy.float)

    # this is where the simulation happens
    while capture < captures:
        totalprop = set_propensities(props,lattice,irate,speeds,ccount,rwidth)
        agent = pick_action(lattice,props,totalprop,ocount*ccount)

        #if agent == -1:del_t = fixed_del_t
        if totalprop < pmin:del_t = fixed_del_t
        else:
            # calculate a timestep
            tpinv = 1.0/totalprop
            del_t = -1.0*log(random.random())*tpinv
            if agent == 0:
                lattice[0] = 1
                for x in range(agent+rwidth):fluxes[x] += 1
                #print "birth an agent at codon site zero"
            elif agent < ccount:
                if not agent > ccount - rwidth:
                    fluxes[agent-1+rwidth] += 1
                lattice[agent-1] = 0
                lattice[agent] = 1
                #print "move an agent along the lattice"
            elif agent == ccount:
                completed += 1
                lattice[ccount-1] = 0
                #print "agent reached the end of the lattice"
            else:
                stalled += 1
                lattice[agent-ccount-1] = 0
                #print "remove agent from lattice ... stalled"

        time += del_t
        iteration += 1

        # set the state of the system
        state[0] = time
        state[1] = iteration
        state[2] = fluxes[0]
        state[3] = fluxes[int(ccount/2)]
        state[4] = completed
        state[5] = stalled

        # perform captures of the state
        realtime = state[0]
        while lasttime < realtime and capture < captures:
            state[0] = lasttime
            lasttime += fixed_del_t
            for tdex in range(ptcount):
                data[tdex,capture] = state[ptdexes[tdex]]
            capture += 1
        state[0] = realtime

    return data



