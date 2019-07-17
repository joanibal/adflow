# dot product test for the heat flux derivatives
from __future__ import print_function
# DO NOT USE THIS SCRIPT AS A REFERENCE FOR HOW TO USE ADFLOW
# THIS SCRIPT USES PRIVATE INTERNAL FUNCTIONALITY THAT IS
# SUBJECT TO CHANGE!!
############################################################

from commonUtils import *
from mdo_regression_helper import *
from baseclasses import AeroProblem
from mpi4py import MPI
import copy
import os
import sys
# sys.path.append(os.path.abspath('../../../adflow_2'))
sys.path.append(os.path.abspath('../../'))
from python.pyADflow import ADFLOW
# from python.pyADflow_C import ADFLOW_C as ADFLOW

# ###################################################################
# DO NOT USE THIS IMPORT STRATEGY! THIS IS ONLY USED FOR REGRESSION
# SCRIPTS ONLY. Use 'from adflow import ADFLOW' for regular scripts.

# ###################################################################

# ****************************************************************************
print('Test Heat: MDO tutorial -- Rans -- Scalar JST')

# ****************************************************************************
# gridFile = '../inputFiles/fc_therm_000_vol.cgns'
# gridFile = '../inputFiles/floating_plate_hot.cgns'
gridFile = '../inputFiles/debug.cgns'
# gridFile = '../inputFiles/mdo_tutorial_viscous_scalar_jst.cgns'
options = {
    'printTiming': False,

    # Common Parameters
    'gridFile': gridFile,
    # 'restartFile': gridFile,
    'outputDirectory': './',

    # 'oversetupdatemode': 'full',
    'volumevariables': ['temp', ],
    'surfacevariables': ['yplus', 'vx', 'vy', 'vz', 'temp', ],
    'monitorVariables': ['resturb', 'yplus'],

    # Physics Parameters
    # 'equationType': 'euler',
    # 'equationType': 'laminar NS',
    'equationType': 'rans',
    # 'vis2':0.0,
    'liftIndex': 2,
    'CFL': 1.0,

    'useANKSolver': False,
    # 'ANKswitchtol': 1e0,
    # 'ankcfllimit': 1e4,
    'anksecondordswitchtol': 1e-3,
    'ankcoupledswitchtol': 1e-6,

    # NK parameters
    'useNKSolver': False,
    'nkswitchtol': 1.0e-7,
    # 'rkreset': True,
    # 'nrkreset': 20,
    'MGCycle': 'sg',
    # 'MGStart': -1,
    # Convergence Parameters
    'L2Convergence': 1e-14,
    'nCycles': 1000,
    'nCyclesCoarse': 250,

}


temp_air = 273  # kelvin
Pr = 0.72
mu = 1.81e-5  # kg/(m * s)
# Rex =
u_inf = 30  # m/s\
p_inf = 101e3
rho_inf = p_inf / (287 * temp_air)

ap = AeroProblem(name='fc_therm', V=u_inf, T=temp_air,
                 P=p_inf, areaRef=1.0, chordRef=1.0, evalFuncs=['heatflux'],
                 alpha=0.0, beta=0.00, xRef=0.0, yRef=0.0, zRef=0.0)


# Create the solver
CFDSolver = ADFLOW(options=options, debug=False)
# res = CFDSolver.getResidual(ap)
# totalR0 = CFDSolver.getFreeStreamResidual(ap)
# res /= totalR0
# print('Norm of residual')
# reducedSum = MPI.COMM_WORLD.reduce(numpy.sum(res**2))
# print(reducedSum)

res = CFDSolver.getResidual(ap)


print('# ---------------------------------------------------#')
print('#                 Dot product Tests                  #')
print('# ---------------------------------------------------#')

# Create a common set of seeds
wDot = CFDSolver.getStatePerturbation(314)
xVDot = CFDSolver.getSpatialPerturbation(314)
dwBar = CFDSolver.getStatePerturbation(314)
fBar = CFDSolver.getSurfacePerturbation(314)
hfBar = CFDSolver.getSurfacePerturbation(314)[:, 1]

print('Dot product test for w -> R')

dwDot = CFDSolver.computeJacobianVectorProductFwd(wDot=wDot, residualDeriv=True)
wBar = CFDSolver.computeJacobianVectorProductBwd(resBar=dwBar, wDeriv=True)

dotLocal1 = numpy.sum(dwDot * dwBar)
dotLocal2 = numpy.sum(wDot * wBar)

reg_par_write_sum(dotLocal1, 1e-10, 1e-10)
reg_par_write_sum(dotLocal2, 1e-10, 1e-10)

print('Dot product test for Xv -> R')
dwDot = CFDSolver.computeJacobianVectorProductFwd(xVDot=xVDot, residualDeriv=True)
xVBar = CFDSolver.computeJacobianVectorProductBwd(resBar=dwBar, xVDeriv=True)

dotLocal1 = numpy.sum(dwDot * dwBar)
dotLocal2 = numpy.sum(xVDot * xVBar)

# For laminar/Rans the DP test for nodes is generally appears a
# little less accurate. This is ok. It has to do with the very
# small offwall spacing on laminar/RANS meshes.
reg_par_write_sum(dotLocal1, 1e-9, 1e-9)
reg_par_write_sum(dotLocal2, 1e-9, 1e-9)


print('Dot product test for w -> HF')
hfDot = CFDSolver.computeJacobianVectorProductFwd(wDot=wDot, hfDeriv=True)
wBar = CFDSolver.computeJacobianVectorProductBwd(hfBar=hfBar, wDeriv=True)
# import ipdb
# ipdb.set_trace()
dotLocal1 = numpy.sum(hfDot.flatten() * hfBar.flatten())
dotLocal2 = numpy.sum(wDot * wBar)

reg_par_write_sum(dotLocal1, 1e-10, 1e-10)
reg_par_write_sum(dotLocal2, 1e-10, 1e-10)

print('Dot product test for xV -> HF')

hfDot = CFDSolver.computeJacobianVectorProductFwd(xVDot=xVDot, hfDeriv=True)
xVBar = CFDSolver.computeJacobianVectorProductBwd(hfBar=hfBar, xVDeriv=True)


dotLocal1 = numpy.sum(hfDot.flatten() * hfBar.flatten())
dotLocal2 = numpy.sum(xVDot * xVBar)

reg_par_write_sum(dotLocal1, 1e-10, 1e-10)
reg_par_write_sum(dotLocal2, 1e-10, 1e-10)


print('Dot product test for w -> F')
wBar = CFDSolver.computeJacobianVectorProductBwd(fBar=fBar, wDeriv=True)
fDot = CFDSolver.computeJacobianVectorProductFwd(wDot=wDot, fDeriv=True)

dotLocal1 = numpy.sum(fDot.flatten() * fBar.flatten())
dotLocal2 = numpy.sum(wDot * wBar)

reg_par_write_sum(dotLocal1, 1e-10, 1e-10)
reg_par_write_sum(dotLocal2, 1e-10, 1e-10)

print('Dot product test for xV -> F')

fDot = CFDSolver.computeJacobianVectorProductFwd(xVDot=xVDot, fDeriv=True)
xVBar = CFDSolver.computeJacobianVectorProductBwd(fBar=fBar, xVDeriv=True)

dotLocal1 = numpy.sum(fDot.flatten() * fBar.flatten())
dotLocal2 = numpy.sum(xVDot * xVBar)

reg_par_write_sum(dotLocal1, 1e-10, 1e-10)
reg_par_write_sum(dotLocal2, 1e-10, 1e-10)

print('Dot product test for (w, xV) -> (dw, F)')

dwDot, fDot = CFDSolver.computeJacobianVectorProductFwd(wDot=wDot, xVDot=xVDot,
                                                        residualDeriv=True, fDeriv=True)
wBar, xVBar = CFDSolver.computeJacobianVectorProductBwd(resBar=dwBar, fBar=fBar,
                                                        wDeriv=True, xVDeriv=True)

dotLocal1 = numpy.sum(dwDot * dwBar) + numpy.sum(fDot.flatten() * fBar.flatten())
dotLocal2 = numpy.sum(wDot * wBar) + numpy.sum(xVDot * xVBar)

reg_par_write_sum(dotLocal1, 1e-10, 1e-10)
reg_par_write_sum(dotLocal2, 1e-10, 1e-10)
