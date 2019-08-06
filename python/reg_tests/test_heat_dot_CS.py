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
import subprocess
import sys
import numpy as np
sys.path.append(os.path.abspath('../../'))
# from python.pyADflow import ADFLOW
from python.pyADflow_C import ADFLOW_C as ADFLOW

# ###################################################################
# DO NOT USE THIS IMPORT STRATEGY! THIS IS ONLY USED FOR REGRESSION
# SCRIPTS ONLY. Use 'from adflow import ADFLOW' for regular scripts.

# ###################################################################

# ****************************************************************************
print('Test Heat: MDO tutorial -- Rans -- Scalar JST')

# ****************************************************************************
# hot plate floating in space (only one sym plane)
# gridFile = '../inputFiles/floating_plate_hot.cgns'
gridFile = '../inputFiles/nacelle_coarsestestest_vol.cgns'

# plate with an array of temperatures
# gridFile = '../inputFiles/plate_array_temp.cgns'

# gridFile = '../inputFiles/debug.cgns'
# gridFile = '../inputFiles/debug_cold.cgns'
# gridFile = '../inputFiles/debug_fric.cgns'

# mdo tutorial wing
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
#

res = CFDSolver.getResidual(ap)
# CFDSolver.callMasterRoutine(ap)
# CFDSolver(ap)

# # Create a common set of seeds
wDot = CFDSolver.getStatePerturbation(314)
xVDot = CFDSolver.getSpatialPerturbation(314)
dwBar = CFDSolver.getStatePerturbation(314)
fBar = CFDSolver.getSurfacePerturbation(314)
hfBar = CFDSolver.getSurfacePerturbation(314)[:, 1]

# # wBar = CFDSolver.computeJacobianVectorProductBwd(resBar=dwBar, wDeriv=True)

# CFDSolver.computeDotProductTest(fBar=True, xVDot=True)
# CFDSolver.computeDotProductTest(hfBar=True, xVDot=True)
# CFDSolver.computeDotProductTest(hfBar=True, wDot=True)
# CFDSolver.computeDotProductTest(hfBar=True, wDot=True)
# quit()

# dwDot1, fDot1, hfDot1 = CFDSolver.computeJacobianVectorProductFwd(
#     wDot=wDot, fDeriv=True, hfDeriv=True, residualDeriv=True)
dwDot2, fDot2, hfDot2 = CFDSolver.computeJacobianVectorProductFwd(
    wDot=wDot, fDeriv=True, hfDeriv=True, residualDeriv=True, mode='CS')
# print(np.real(fDot2))
# print(np.real(hfDot2))
