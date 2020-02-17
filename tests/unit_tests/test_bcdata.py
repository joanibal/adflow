from __future__ import print_function
import unittest
import numpy
from baseclasses import AeroProblem
import sys
import os
baseDir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(baseDir,'../../'))
from python.pyADflow import ADFLOW
class BasicTests(unittest.TestCase):

    def setUp(self):
        gridFile = os.path.join(baseDir, '../input_files/conic_conv_nozzle_mb_L4.cgns')
        options = {'gridfile': gridFile}
        self.CFDSolver = ADFLOW(options=options)
        self.CFDSolver.addFamilyGroup('upstream',['inlet'])
        self.CFDSolver.addFamilyGroup('downstream',['outlet'])
        
        # Aerodynamic problem description (taken from test 17)
        ap = AeroProblem(name='nozzle', alpha=90, mach=0.5, altitude=0,
                        areaRef=1.0, chordRef=1.0, R=287.87,
                        evalFuncs=['mdot_up', 'mdot_down', 'mdot_plane',
                                    'mavgptot_up', 'mavgptot_down', 'mavgptot_plane',
                                    'mavgttot_up', 'mavgttot_down', 'mavgttot_plane',
                                    'mavgps_up', 'mavgps_down', 'mavgps_plane',
                                    ])

        # these values come from inspecting the cgns mesh its self
        self.intial_dvs = {
            'downstream': {
                'Pressure': 79326.7,
            },
            'upstream': {
                'PressureStagnation': 100000.0,
                'TemperatureStagnation': 500
            }
        }

        ap.setBCVar('Pressure',  self.intial_dvs['downstream']['Pressure'] + 100 , 'downstream')
        ap.addDV('Pressure', family='downstream')

        ap.setBCVar('PressureStagnation',   self.intial_dvs['upstream']['PressureStagnation']  + 200, 'upstream')
        ap.addDV('PressureStagnation', family='upstream')

        ap.setBCVar('TemperatureStagnation',  self.intial_dvs['upstream']['TemperatureStagnation'] + 50, 'upstream')
        ap.addDV('TemperatureStagnation', family='upstream')
        self.AP = ap 

    def test_getBCData(self):
        "Tests if mesh was read properly"
        print('im, here')

        for group in self.intial_dvs.keys():
            BCDataArray, BCVarNames, BCPatchData = self.CFDSolver.getBCData(groupName=group)

            for idx in range(len(BCDataArray)):
                for var in self.intial_dvs[group].keys():
                    idx_p = BCVarNames[idx].index(var)
                    assert BCDataArray[idx][idx_p] == self.intial_dvs[group][var]
                

if __name__ == '__main__':
    unittest.main()
