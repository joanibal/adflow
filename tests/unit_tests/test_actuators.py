from __future__ import print_function
import unittest
import numpy as np
from baseclasses import AeroProblem
import sys
import os
baseDir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(baseDir,'../../'))
from python.pyADflow import ADFLOW
from pprint import pprint as pp
import copy


class ActuatorTests(unittest.TestCase):

    def setUp(self):
        gridFile = os.path.join(baseDir, '../input_files/pipe.cgns')
        self.options = {'gridfile': gridFile,
                        'mgcycle':'sg',
                        'ncycles':1000,
                        'useanksolver':True,
                        'usenksolver':True,
                        'volumevariables': ['temp', 'mach', 'resrho' ],
                        'surfacevariables':['temp', 'vx', 'vy', 'vz', 'p', 'ptloss', 'mach', 'rho'],
                        'equationType':'euler',
                        'l2convergence': 1e-13}
                        
        CFDSolver = ADFLOW(options=self.options)

        # intSurfFile = os.path.join(baseDir, '../input_files/int_plane.xyz')
        # CFDSolver.addIntegrationSurface(intSurfFile, 'plane')
        # CFDSolver.finalizeUserIntegrationSurfaces()

        CFDSolver.addFunction('mdot', 'inlet', name="mdot_in")
        CFDSolver.addFunction('mdot', 'outlet', name="mdot_out")
        # CFDSolver.addFunction('mdot', 'plane', name="mdot_plane")

        CFDSolver.addFunction('mavgptot', 'outlet', name="mavgptot_out")
        CFDSolver.addFunction('mavgptot', 'inlet', name="mavgptot_in")
        # CFDSolver.addFunction('mavgptot', 'plane', name="mavgptot_plane")

        CFDSolver.addFunction('aavgptot', 'outlet', name="aavgptot_out")
        CFDSolver.addFunction('aavgptot', 'inlet', name="aavgptot_in")
        # CFDSolver.addFunction('aavgptot', 'plane', name="aavgptot_plane")

        CFDSolver.addFunction('mavgttot', 'outlet', name="mavgttot_out")
        CFDSolver.addFunction('mavgttot', 'inlet', name="mavgttot_in")
        # CFDSolver.addFunction('mavgttot', 'plane', name="mavgttot_plane")

        CFDSolver.addFunction('mavgps', 'outlet', name="mavgps_out")
        CFDSolver.addFunction('mavgps', 'inlet', name="mavgps_in")
        # CFDSolver.addFunction('mavgps', 'plane', name="mavgps_plane")

        CFDSolver.addFunction('aavgps', 'outlet', name="aavgps_out")
        CFDSolver.addFunction('aavgps', 'inlet', name="aavgps_in")
        # CFDSolver.addFunction('aavgps', 'plane', name="aavgps_plane")

        CFDSolver.addFunction('area', 'inlet', name="area_in")
        CFDSolver.addFunction('area', 'outlet', name="area_out")

        CFDSolver.addFunction('mavgvx', 'inlet', name="vx_in")
        CFDSolver.addFunction('mavgvx', 'outlet', name="vx_out")

        CFDSolver.addFunction('mavgps', 'inlet', name="ps_in")
        CFDSolver.addFunction('mavgps', 'outlet', name="ps_out")

        # CFDSolver.addFunction('area', 'plane', name="area_plane")

        self.CFDSolver = CFDSolver

        
        self.ap = AeroProblem(name='actuator_in_pipe', alpha=00, mach=0.6, altitude=0,
                        areaRef=1.0, chordRef=1.0,
                    evalFuncs=['mdot_in', 'mdot_out', 
                               'mavgptot_in', 'mavgptot_out', 
                               'mavgttot_in', 'mavgttot_out', 
                               'mavgps_in', 'mavgps_out',
                               'area_in', 'area_out', 
                               'aavgps_in', 'aavgps_out',
                               'aavgptot_in', 'aavgptot_out',
                               'ps_in', 'ps_out',
                               'vx_in', 'vx_out'] )


        self.force = 600
        actuatorFile = os.path.join(baseDir, '../input_files/actuator_disk.xyz')
        self.CFDSolver.addActuatorRegion(actuatorFile, np.array([0,0,0]),np.array([1,0,0]), 'actuator', thrust=self.force )

    def test_mom_flux(self):
        "Tests if the correct amount of momentum is added to the flow by the actuator"

        self.CFDSolver(self.ap)
        funcs = {}
        self.CFDSolver.evalFunctions(self.ap, funcs)
        pp(funcs)
            

        force_cfd =  -funcs[self.ap.name + '_mdot_out']*funcs[self.ap.name + '_vx_out'] + funcs[self.ap.name + '_ps_out']*funcs[self.ap.name +'_area_out'] \
                 - ( funcs[self.ap.name + '_mdot_in']*funcs[self.ap.name + '_vx_in'] + funcs[self.ap.name + '_ps_in']*funcs[self.ap.name +'_area_in'])

        # print('analytic', self.force, 'cforced', force_cfd)
        # print('% error', (force_cfd - force)/self.force)

        np.testing.assert_allclose(force_cfd, self.force, rtol=1e-3)

    def test_set_thrust(self):
        "Tests if the correct amount of momentum is added to the flow by the actuator"

        force = 100

        ap = copy.deepcopy(self.ap)
        ap.setActuatorVar('Thrust',  force, 'actuator')
        ap.addDV('Thrust', familyGroup='actuator', name='actuator_thrust')
        self.CFDSolver(ap)
        funcs = {}
        self.CFDSolver.evalFunctions(ap, funcs)
        pp(funcs)
            

        force_cfd =  -funcs[self.ap.name + '_mdot_out']*funcs[self.ap.name + '_vx_out'] + funcs[self.ap.name + '_ps_out']*funcs[self.ap.name +'_area_out'] \
                 - ( funcs[self.ap.name + '_mdot_in']*funcs[self.ap.name + '_vx_in'] + funcs[self.ap.name + '_ps_in']*funcs[self.ap.name +'_area_in'])

        print('analytic', force, 'cforced', force_cfd)
        print('% error', (force_cfd - force)/force)

        np.testing.assert_allclose(force_cfd, force, rtol=1e-3)




if __name__ == '__main__':
    unittest.main()
