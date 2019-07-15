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
sys.path.append(os.path.abspath('../../'))
from python.pyADflow import ADFLOW

# ###################################################################
# DO NOT USE THIS IMPORT STRATEGY! THIS IS ONLY USED FOR REGRESSION
# SCRIPTS ONLY. Use 'from adflow import ADFLOW' for regular scripts.

# ###################################################################

# ****************************************************************************
print('Test Heat: MDO tutorial -- Rans -- Scalar JST')

# ****************************************************************************
gridFile = '../inputFiles/floating_plate_hot.cgns'

options = {
    'printTiming': False,

    # Common Parameters
    'gridFile': gridFile,
    'outputDirectory': './',

    # 'oversetupdatemode': 'full',
    'volumevariables': ['temp', ],
    'surfacevariables': ['yplus', 'vx', 'vy', 'vz', 'temp', 'heatflux'],
    'monitorVariables': ['heatflux', 'resturb', 'yplus'],

    # Physics Parameters
    # 'equationType': 'euler',
    # 'equationType': 'laminar NS',
    'equationType': 'rans',
    # 'vis2':0.0,
    'liftIndex': 2,
    'CFL': 1.0,

    'useANKSolver': True,
    # 'ANKswitchtol': 1e0,
    # 'ankcfllimit': 1e4,
    'anksecondordswitchtol': 1e-3,
    'ankcoupledswitchtol': 1e-6,

    # NK parameters
    'useNKSolver': True,
    'nkswitchtol': 1.0e-7,
    # 'rkreset': True,
    # 'nrkreset': 20,
    'MGCycle': 'sg',
    # 'MGStart': -1,
    # Convergence Parameters
    'L2Convergence': 1e-12,
    'nCycles': 4000,
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


# Check the residual
res = CFDSolver.getResidual(ap)
totalR0 = CFDSolver.getFreeStreamResidual(ap)
res /= totalR0
print('Norm of residual')
reg_par_write_norm(res, 1e-10, 1e-10)

funcs = {}
CFDSolver.evalFunctions(ap, funcs)
print('Eval Functions:')
reg_root_write_dict(funcs, 1e-10, 1e-10)

wDot = CFDSolver.getStatePerturbation(314)

# print('# ---------------------------------------------------#')
# print('#             Forward mode testing                   #')
# print('# ---------------------------------------------------#')

# resDot, funcsDot, fDot = CFDSolver.computeJacobianVectorProductFwd(
#     wDot=wDot, residualDeriv=True, funcDeriv=True, hfDeriv=True)

# print('||dR/dw * wDot||')
# reg_par_write_norm(resDot, 1e-10, 1e-10)

# print('dFuncs/dw * wDot')
# reg_root_write_dict(funcsDot, 1e-10, 1e-10)

# print('||dF/dw * wDot||')
# reg_par_write_norm(fDot, 1e-10, 1e-10)


# print('-> Derivatives with respect to nodes')
# xVDot = CFDSolver.getSpatialPerturbation(314)

# resDot, funcsDot, fDot = CFDSolver.computeJacobianVectorProductFwd(
#     xVDot=xVDot, residualDeriv=True, funcDeriv=True, fDeriv=True)

# print('||dR/dXv * xVDot||')
# reg_par_write_norm(resDot, 1e-10, 1e-10)

# # These can be finiky sometimes so a bigger tolerance.
# print('dFuncs/dXv * xVDot')
# reg_root_write_dict(funcsDot, 1e-9, 1e-9)

# print('||dF/dXv * xVDot||')
# reg_par_write_norm(fDot, 1e-10, 1e-10)

print('# ---------------------------------------------------#')
print('#             Reverse mode testing                   #')
print('# ---------------------------------------------------#')

# print('-> Res bar Seed')
# dwBar = CFDSolver.getStatePerturbation(314)

# wBar, xVBar, xDvBar = CFDSolver.computeJacobianVectorProductBwd(
#     resBar=dwBar, wDeriv=True, xVDeriv=True, xDvDerivAero=True)

# print('||dwBar^T * dR/dw||')
# reg_par_write_norm(wBar, 1e-10, 1e-10)

# print('||dwBar^T * dR/dXv||')
# norm = CFDSolver.getUniqueSpatialPerturbationNorm(xVBar)
# reg_root_write(norm, 1e-10, 1e-10)

# print('||dwBar^T * dR/xDv||')
# reg_root_write_dict(xDvBar, 1e-10, 1e-10)

# print('-> F Bar Seed')
# fBar = CFDSolver.getSurfacePerturbation(314)

# wBar, xVBar, xDvBar = CFDSolver.computeJacobianVectorProductBwd(
# fBar = fBar, wDeriv = True, xVDeriv = True, xDvDerivAero = True)

# print('||FBar^T * dF/dw||')
# reg_par_write_norm(wBar, 1e-10, 1e-10)

# print('||FBar^T * dF/dXv||')
# norm = CFDSolver.getUniqueSpatialPerturbationNorm(xVBar)
# reg_root_write(norm, 1e-10, 1e-10)

# print('||FBar^T * dF/xDv||')
# reg_root_write_dict(xDvBar, 1e-10, 1e-10)

print('-> HF Bar Seed')
hfBar = CFDSolver.getSurfacePerturbation(314)[:, 1]

wBar, xVBar, xDvBar = CFDSolver.computeJacobianVectorProductBwd(
    hfBar=hfBar, wDeriv=True, xVDeriv=True, xDvDerivAero=True)

print('||HFBar^T * dHF/dw||')
reg_par_write_norm(wBar, 1e-10, 1e-10)

print('-> HF Bar Seed')
hfBar = CFDSolver.getSurfacePerturbation(314)[:, 1]

wBar, xVBar, xDvBar = CFDSolver.computeJacobianVectorProductBwd(
    hfBar=hfBar, wDeriv=True, xVDeriv=True, xDvDerivAero=True)

print('||HFBar^T * dHF/dw||')
reg_par_write_norm(wBar, 1e-10, 1e-10)

# print('||HFBar^T * dHF/dXv||')
# norm = CFDSolver.getUniqueSpatialPerturbationNorm(xVBar)
# reg_root_write(norm, 1e-10, 1e-10)

# print('||HFBar^T * dHF/xDv||')
# reg_root_write_dict(xDvBar, 1e-10, 1e-10)
