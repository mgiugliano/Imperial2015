#!/usr/bin/env python

####
#  This script is used to simulate a stochastic HH neuron model incorporating channel noise.
#  It uses configuration files (see the included example) to determine the characteristics
#  of the neuron to be simulated.
#
#  Author: Daniele Linaro - daniele.linaro@unige.it
#  Date: September 2010
####


###
#  FUNCTIONS
###

def printHelp(scriptname):
    print('\nUsage: ' + scriptname + ' [options]\n')
    print('\twhere options are:\n');
    print('\t\t-f --conf-file: parse the specified configuration file.\n')
    print('\t\t-h --help: print this help message.\n')
    print('\nAuthor: Daniele Linaro -- daniele.linaro@unige.it\n')

def checkDuplicates(filename):
    import os
    s = filename.split('.')
    if len(s) > 1:
        s = [filename[0:-(len(s[-1])+1)], s[-1]]
    suffix = 0
    while os.path.isfile(filename):
        suffix = suffix + 1
        filename = s[0] + '_' + str(suffix)
        if len(s) > 1:
            filename = filename + '.' + s[1]
    return filename

###

from sys import argv, stdout
import os
import time

if len(argv) == 2 and (argv[1] == '-h' or argv[1] == '--help'):
    printHelp(argv[0])
    exit(1)

conffile = 'neuron.cfg'
if len(argv) > 2:
    if argv[1] == '-f' or argv[1] == '--config-file':
        conffile = argv[2]
    else:
        printHelp(argv[0])
        exit(1)

from ConfigParser import ConfigParser
if not os.path.isfile(conffile):
    printHelp(argv[0])
    print(conffile + ' is not a valid configuration file. Aborting...')
    exit(1)
fid = open(conffile,'r')
config = ConfigParser()
config.readfp(fid)
fid.close()

from neuron import h
from nrn import Section
import numpy

## stimulus properties
stimprop = {'delay': config.getfloat('Stimulus','delay'),
	    'dur': config.getfloat('Stimulus','duration'),
	    'amp': config.getfloat('Stimulus','amplitude')}

## soma properties
type = config.get('Neuron','type')
L = config.getfloat('Neuron','length')
gamma_na = config.getfloat('Neuron','gamma_na')
gamma_k = config.getfloat('Neuron','gamma_k')
diam = config.getfloat('Neuron','diameter')
Ra = 35.4

## stuff
savedata = config.get('Misc','savedata').lower()
savespikes = config.get('Misc','savespikes').lower()
showplot = config.get('Misc','showplot').lower();
saveplot = config.get('Misc','saveplot').lower()

## build the neuron
soma = Section()     # Create the soma of our single-compartment model
soma.L = L           # [um] length of the soma
soma.diam = diam     # [um] diameter of the soma
soma.nseg = 1        # number of segments: 1 means a single-compartment model
soma.cm = 1.0        # [uF] membrane capacitance
soma.Ra = Ra         # [Ohm cm] citoplasmic resistance
if type.lower() == 'effective':
    soma.insert('HHcn')  # insert effective HH mechanism
    mechanism = soma(0.5).HHcn
elif type.lower() == 'microscopic':
    soma.insert('HHmicro')  # insert microscopic HH mechanism
    mechanism = soma(0.5).HHmicro
else:
    print('Unknown neuron type: ' + type)
    exit(1)
mechanism.gamma_na = gamma_na
mechanism.gamma_k = gamma_k
mechanism.seed = round(numpy.random.rand() * 10000)

## stimulus
stimulus = h.IClamp(soma(0.5))
stimulus.delay = stimprop['delay']
stimulus.dur = stimprop['dur']
stimulus.amp = stimprop['amp']
dt = 0.001
tstop = stimprop['dur'] + 2*stimprop['delay']

## recorders for time and membrane voltage
if savedata == 'yes' or showplot == 'yes':
    t = numpy.arange(0,tstop+dt,dt)
    recorders = {}
    labels = ['v']
    for lbl in labels:
        recorders[lbl] = h.Vector(t)
    recorders['v'].record(soma(0.5)._ref_v)

## a spike counter
apc = h.APCount(soma(0.5))
apc.thresh = 0
spikes = h.Vector()
apc.record(spikes)

## run the simulation
h.load_file('stdrun.hoc')
h.tstop = tstop
h.dt = dt

stdout.write('\n\n||| Started on ' + time.strftime('%a, %d %b %Y %H:%M:%S', time.localtime()) + '. |||\n\n')
stdout.write('>>> Simulating the ' + type + ' HH model for ' + str(h.tstop) + ' ms. <<<\n')
stdout.flush()

start = time.time()
h.run()
stop = time.time()

stdout.write('>>> Elapsed time: ' + str(int(stop-start)) + ' seconds. <<<\n')
stdout.flush()

## save the results
id = '_' + type + '_amp='+('%.3f' % stimprop['amp']) + \
     '_gamma_na='+str(gamma_na) + '_gamma_k='+str(gamma_k) + \
     '_L='+str(L) + '_diam='+str(diam)

if savedata == 'yes':
    ds = 10
    suffix = 0
    datafile = checkDuplicates('voltage' + id + '.dat')
    data = numpy.zeros((len(recorders)+1,len(t[::ds])))
    data[0,:] = t[::ds]
    for k,lbl in enumerate(labels):
        data[k+1,:] = numpy.array(recorders[lbl])[::ds]
    numpy.savetxt(datafile,data.transpose(),'%.4f')

if savespikes == 'yes':
    suffix = 0
    spikesfile = checkDuplicates('spikes' + id + '.dat')
    numpy.savetxt(spikesfile, numpy.array(spikes), '%.6f')

## plot the results
if showplot == 'yes' or saveplot == 'yes':
    from pylab import figure, plot, show, xlabel, ylabel
    figure()
    plot(t,recorders['v'],'k')
    xlabel('t (ms)')
    ylabel('V (mV)')
    show()
    if saveplot == 'yes':
        from pylab import savefig
        suffix = 0
        plotfile = checkDuplicates('voltage' + id + '.eps')
        savefig(plotfile)

stdout.write('\n||| Ended on ' + time.strftime('%a, %d %b %Y %H:%M:%S', time.localtime()) + '. |||\n\n')
stdout.flush()
