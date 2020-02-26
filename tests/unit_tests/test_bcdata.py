from __future__ import print_function
import unittest
import numpy as np
import copy
from baseclasses import AeroProblem
import sys
import os
baseDir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(baseDir,'../../'))
from python.pyADflow import ADFLOW

class BCTests(unittest.TestCase):

    def setUp(self):
        gridFile = os.path.join(baseDir, '../input_files/conic_conv_nozzle_mb_L4.cgns')
        self.options = {'gridfile': gridFile}

        
        # Aerodynamic problem description (taken from test 17)
        ap = AeroProblem(name='nozzle', alpha=90, mach=0.5, altitude=0,
                        areaRef=1.0, chordRef=1.0, R=287.87,
                        evalFuncs=['mdot_up', 'mdot_down', 'mdot_plane',
                                    'mavgptot_up', 'mavgptot_down', 'mavgptot_plane',
                                    'mavgttot_up', 'mavgttot_down', 'mavgttot_plane',
                                    'mavgps_up', 'mavgps_down', 'mavgps_plane',
                                    ])

        # these values come from inspecting the cgns mesh itself
        self.initial_dvs = {
            'downstream': {
                'Pressure': 79326.7,
            },
            'upstream': {
                'PressureStagnation': 100000.0,
                'TemperatureStagnation': 500
            }
        }

        self.set_dvs = {
            'downstream': {
                'Pressure': 89326.7,
            },
            'upstream': {
                'PressureStagnation': 100000.0,
                'TemperatureStagnation': 500
            }
        }

        ap.setBCVar('Pressure',  self.set_dvs['downstream']['Pressure'], 'downstream')
        # ap.addDV('Pressure', familyGroup='downstream')

        ap.setBCVar('PressureStagnation',   self.set_dvs['upstream']['PressureStagnation'], 'upstream')
        # ap.addDV('PressureStagnation', familyGroup='upstream')

        ap.setBCVar('TemperatureStagnation',  self.set_dvs['upstream']['TemperatureStagnation'], 'upstream')
        # ap.addDV('TemperatureStagnation', familyGroup='upstream')
        self.ap = ap 
    
    def test_getBCData(self):
        "Tests if mesh was read properly"
        CFDSolver = ADFLOW(options=self.options)
        CFDSolver.addFamilyGroup('upstream',['inlet'])
        CFDSolver.addFamilyGroup('downstream',['outlet'])

        BCData = CFDSolver.getBCData(groupNames=self.initial_dvs.keys())

        for group in self.initial_dvs:
            for BCVar in self.initial_dvs[group]:
                for patch in BCData[group][BCVar]:
                    np.testing.assert_array_equal(BCData[group][BCVar][patch], \
                                                  np.array([self.initial_dvs[group][BCVar]]))
     

    def test_setBCData(self):
        "Test if bcdata is correctly set to data arrays"

        CFDSolver = ADFLOW(options=self.options)
        CFDSolver.addFamilyGroup('upstream',['inlet'])
        CFDSolver.addFamilyGroup('downstream',['outlet'])

        CFDSolver.setAeroProblem(self.ap)
        BCData = CFDSolver.getBCData(groupNames=self.set_dvs.keys())

        for group in self.set_dvs:
            for BCVar in self.set_dvs[group]:
                for patch in BCData[group][BCVar]:
                    np.testing.assert_array_equal(BCData[group][BCVar][patch], \
                                                  np.array([self.set_dvs[group][BCVar]]))


    def test_getBCData_array(self):
        pass

    def test_setBCData_array(self):
        pass


class BCDVTests(unittest.TestCase):

    def setUp(self):
        gridFile = os.path.join(baseDir, '../input_files/conic_conv_nozzle_mb_L4.cgns')
        self.options = {'gridfile': gridFile}

        ap = AeroProblem(name='nozzle_flow', alpha=90, mach=0.5, altitude=0,
                    areaRef=1.0, chordRef=1.0, R=287.87,
                    evalFuncs=['mdot_up', 'mdot_down', 'mdot_plane',
                                'mavgptot_up', 'mavgptot_down', 'mavgptot_plane',
                                'mavgttot_up', 'mavgttot_down', 'mavgttot_plane',
                                'mavgps_up', 'mavgps_down', 'mavgps_plane',
                                ])

        self.ap = ap 


    def test_addBCVar_in_AP(self):
        CFDSolver = ADFLOW(options=self.options)

        ap = copy.copy(self.ap)

        ap.setBCVar('Pressure',  70000.0, 'outlet')
        ap.addDV('Pressure', familyGroup='outlet', name='outlet_pressure')
        
        DVs = {'outlet_pressure':123000.0}

        ap.setDesignVars(DVs)
        BCData = ap.getBCData()
        assert BCData['outlet']['Pressure'] == DVs['outlet_pressure']


        CFDSolver.setAeroProblem(ap)
        BCData = CFDSolver.getBCData(groupNames=['outlet'])

        for patch in BCData['outlet']['Pressure']:
            np.testing.assert_array_equal(BCData['outlet']['Pressure'][patch], \
                                            np.array([DVs['outlet_pressure']]))

        
    def test_addBCVar_in_AP_dict(self):
        CFDSolver = ADFLOW(options=self.options)

        ap = copy.copy(self.ap)
        group = 'outlet'
        BCVar = 'Pressure'

        BCData = CFDSolver.getBCData(groupNames=[group])
        ap.setBCVar(BCVar,  BCData[group][BCVar], group)
        ap.addDV(BCVar, familyGroup=group, name='outlet_pressure')
        
        DVs = {}
        for ii, DV in enumerate(ap.DVs):
            DVs[DV] = ap.DVs[DV].value + 200 + ii

        ap.setDesignVars(DVs)

        BCData = ap.getBCData()

        for ii, DV in enumerate(ap.DVs):
            key = ap.DVs[DV].dict_key
            assert BCData[group][BCVar][key] == DVs[DV]


        CFDSolver.setAeroProblem(ap)
        BCData = CFDSolver.getBCData(groupNames=[group])

        for ii, DV in enumerate(ap.DVs):
            key = ap.DVs[DV].dict_key
            assert BCData[group][BCVar][key] == DVs[DV]

        
class BCDerivsTests(unittest.TestCase):
    def setUp(self):
        gridFile = os.path.join(baseDir, '../input_files/conic_conv_nozzle_mb_L4.cgns')
        self.options = {'gridfile': gridFile}

        ap = AeroProblem(name='nozzle_flow', alpha=90, mach=0.5, altitude=0,
                    areaRef=1.0, chordRef=1.0, R=287.87,
                    evalFuncs=['mdot_up', 'mdot_down', 'mdot_plane',
                                'mavgptot_up', 'mavgptot_down', 'mavgptot_plane',
                                'mavgttot_up', 'mavgttot_down', 'mavgttot_plane',
                                'mavgps_up', 'mavgps_down', 'mavgps_plane',
                                ])

        self.ap = ap 


    def test_fwd(self):
        CFDSolver = ADFLOW(options=self.options)

        ap = copy.copy(self.ap)

        ap.setBCVar('Pressure',  70000.0, 'outlet')
        ap.addDV('Pressure', familyGroup='outlet', name='outlet_pressure')
        key = 'outlet_pressure'
        xDvDot = {key:1.0}
        CFDSolver.setAeroProblem(ap)
        resDot, funcsDot, fDot = CFDSolver.computeJacobianVectorProductFwd(
            xDvDot=xDvDot, residualDeriv=True, funcDeriv=True, fDeriv=True)

        # handler.root_print('||dR/d%s||'%key)
        # handler.par_add_norm(resDot, 1e-10, 1e-10)

        # handler.root_print('dFuncs/d%s'%key)
        # handler.root_add_dict(funcsDot, 1e-10, 1e-10)

        # handler.root_print('||dF/d%s||'%key)
        # handler.par_add_norm(fDot, 1e-10, 1e-10)



    def test_fwd_dict(self):
        pass

    def test_bwd(self):
        pass


    def test_bwd_dict(self):
        pass



if __name__ == '__main__':
    unittest.main()
