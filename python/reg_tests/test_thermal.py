from __future__ import print_function, division, absolute_import

import numpy as np
import unittest

from adflow import ADFLOW
from baseclasses import AeroProblem



ap = AeroProblem(name='fc_thermal_', V=100, T=273, Pr=1.0,
                 P=101e3, areaRef=1.0, chordRef=1.0,  evalFuncs=['cl', 'cd', 'heatflux'])



aeroOptions = {
    # 'printTiming': False,

    # Common Parameters
    'gridFile': gridFile,
    'outputDirectory': './',
    'discretization': 'upwind',

    'volumevariables': ['temp'],
    'surfacevariables': ['temp', 'heatflux', 'ch', 'vx', 'vy', 'vz', 'cf'],
    'monitorVariables':	['heatflux', 'resturb', 'cl', 'cd', 'yplus'],


    # Physics Parameters
    # 'equationType': 'laminar NS',
    'equationType': 'rans',
    # 'vis2':0.0,
    'liftIndex': 2,
    'CFL': 1.0,
    # 'smoother': 'DADI',
    # 'smoother': 'runge',

    'useANKSolver': True,
    # 'ANKswitchtol': 1e0,
    'ankcfllimit': 5e6,
    'anksecondordswitchtol': 5e-3,
    # 'ankcoupledswitchtol': 5e-8,
    # NK parameters
    'useNKSolver': False,
    'nkswitchtol': 1e-11,
    # 'rkreset': True,
    'nrkreset': 40,
    'MGCycle': 'sg',
    # 'MGStart': -1,
    # Convergence Parameters
    'L2Convergence': 1e-9,
    'nCycles': 6000,  # 000,
    'nCyclesCoarse': 50,
    'ankcfllimit': 5e3,
    'nsubiterturb': 5,
    'ankphysicallstolturb': 0.99,
    'anknsubiterturb': 5,
    # 'ankuseturbdadi': False,
    'ankturbkspdebug': True,

    # 'storerindlayer': True,
    # Turbulence model
    'eddyvisinfratio': .210438,
    'useft2SA': False,
    'turbulenceproduction': 'vorticity',
    'useblockettes': False,


    'forcesastractions': False,
}



class TestIsothermalWall(unittest.TestCase):

    def setUp(self):
        self.wall_temps = np.linspace(100, 900, 9)
        self.CFDSolver = ADFLOW(options=aeroOptions)

    def test_getWallTemp(self):
        temps = self.CFDSolver.getWallTemperature()
        np.testing.assert_array_almost_equal(self.wall_temps, temps, decimal=8)

    def test_setWallTemp(self):
        temps = np.linspace(1000, 9000, 9)
        self.CFDSolver.setWallTemperature(temps)

        #to insure won't be overwrited during __call__()
        self.CFDSolver.setAeroProblem(ap)

        np.testing.assert_array_almost_equal(self.CFDSolver.getWallTemperature(), temps, decimal=8)

    # def test_getHeatFluxes(self):
    #     fluxes = self.CFDSolver.getheatFluxes()
    #     np.testing.assert_array_almost_equal(, fluxes, decimal=8)





temp_wall = 323  # kelvin (about 50 c)calc
# CFDSolver(ap, wallTemperature=temp_wall)
CFDSolver(ap)
funcs = {}
CFDSolver.evalFunctions(ap, funcs)
print(funcs)