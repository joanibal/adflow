from __future__ import print_function
import unittest
import numpy as np
import copy
from baseclasses import AeroProblem
import sys
from pprint import pprint as pp

import os
baseDir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(baseDir,'../../'))
from python.pyADflow import ADFLOW

class BCGetSetTests(unittest.TestCase):

    def setUp(self):
        gridFile = os.path.join(baseDir, '../input_files/conic_conv_nozzle_mb_L4.cgns')
        self.options = {'gridfile': gridFile}

        
        # Aerodynamic problem description (taken from test 17)
        ap = AeroProblem(name='nozzle', alpha=90, mach=0.5, altitude=0,
                        areaRef=1.0, chordRef=1.0, R=287.87)

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
        gridFile = os.path.join(baseDir, '../input_files/conic_conv_nozzle_mb_L4_array.cgns')
        options = copy.deepcopy(self.options)
        options['gridfile'] =  gridFile
        CFDSolver = ADFLOW(options=options)
        CFDSolver.addFamilyGroup('upstream',['inlet'])
        CFDSolver.addFamilyGroup('downstream',['outlet'])

        BCData = CFDSolver.getBCData(groupNames=self.initial_dvs.keys())

        for group in self.initial_dvs:
            for BCVar in self.initial_dvs[group]:
                for patch in BCData[group][BCVar]:
                    np.testing.assert_array_equal(BCData[group][BCVar][patch], \
                                                np.ones(BCData[group][BCVar][patch].size)*
                                                self.initial_dvs[group][BCVar])
     


    def test_setBCData_array(self):
        gridFile = os.path.join(baseDir, '../input_files/conic_conv_nozzle_mb_L4_array.cgns')
        options = copy.deepcopy(self.options)
        options['gridfile'] =  gridFile
        CFDSolver = ADFLOW(options=options)
        CFDSolver.addFamilyGroup('upstream',['inlet'])
        CFDSolver.addFamilyGroup('downstream',['outlet'])

        BCData = CFDSolver.getBCData(groupNames=self.initial_dvs.keys())

        ap = AeroProblem(name='nozzle', alpha=90, mach=0.5, altitude=0,
                        areaRef=1.0, chordRef=1.0, R=287.87)
        for group in self.set_dvs:
            for BCVar in self.set_dvs[group]:
                for patch in BCData[group][BCVar]:
                    BCData[group][BCVar][patch] = np.ones(BCData[group][BCVar][patch].shape)*self.set_dvs[group][BCVar]

            ap.setBCVar(BCVar,  BCData[group][BCVar], group)

        CFDSolver.setAeroProblem(ap)

        for group in self.set_dvs:
            for BCVar in self.set_dvs[group]:
                for patch in BCData[group][BCVar]:
                    np.testing.assert_array_equal(BCData[group][BCVar][patch], \
                                                np.ones(BCData[group][BCVar][patch].size)*
                                                self.set_dvs[group][BCVar])
     

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

        ap = copy.deepcopy(self.ap)

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

        ap = copy.deepcopy(self.ap)
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
        import petsc4py

        petsc4py.init(arch='real-debug-gfortran-3.7.7')
        from petsc4py import PETSc

        gridFile = os.path.join(baseDir, '../input_files/conic_conv_nozzle_mb_L4.cgns')
        intSurfFile = os.path.join(baseDir, '../input_files/integration_plane_viscous.fmt')

        self.options = {'gridfile': gridFile}

        self.options['mgstartlevel'] = 1

        self.options['ncycles'] = 1
        self.options['useanksolver'] = False
        self.options['l2convergence'] = 1e-12
        ap = AeroProblem(name='nozzle_flow', alpha=90, mach=0.5, altitude=0,
                    areaRef=1.0, chordRef=1.0, R=287.87,
                    evalFuncs=['mdot_up', 'mdot_down', 'mdot_plane',
                                'mavgptot_up', 'mavgptot_down', 'mavgptot_plane',
                                'mavgttot_up', 'mavgttot_down', 'mavgttot_plane',
                                'mavgps_up', 'mavgps_down', 'mavgps_plane',
                                ])

        self.ap = ap 

        CFDSolver = ADFLOW(options=self.options)
        CFDSolver.addIntegrationSurface(intSurfFile, 'viscous_plane')
        CFDSolver.finalizeUserIntegrationSurfaces()

        CFDSolver.addFamilyGroup('upstream',['inlet'])
        CFDSolver.addFamilyGroup('downstream',['outlet'])
        CFDSolver.addFamilyGroup('all_flow',['inlet', 'outlet'])
        CFDSolver.addFamilyGroup('output_fam',['all_flow', 'allWalls'])

        CFDSolver.addFunction('mdot', 'upstream', name="mdot_up")
        CFDSolver.addFunction('mdot', 'downstream', name="mdot_down")
        CFDSolver.addFunction('mdot', 'viscous_plane', name="mdot_plane")

        CFDSolver.addFunction('mavgptot', 'downstream', name="mavgptot_down")
        CFDSolver.addFunction('mavgptot', 'upstream', name="mavgptot_up")
        CFDSolver.addFunction('mavgptot', 'viscous_plane', name="mavgptot_plane")

        CFDSolver.addFunction('aavgptot', 'downstream', name="aavgptot_down")
        CFDSolver.addFunction('aavgptot', 'upstream', name="aavgptot_up")
        CFDSolver.addFunction('aavgptot', 'viscous_plane', name="aavgptot_plane")

        CFDSolver.addFunction('mavgttot', 'downstream', name="mavgttot_down")
        CFDSolver.addFunction('mavgttot', 'upstream', name="mavgttot_up")
        CFDSolver.addFunction('mavgttot', 'viscous_plane', name="mavgttot_plane")

        CFDSolver.addFunction('mavgps', 'downstream', name="mavgps_down")
        CFDSolver.addFunction('mavgps', 'upstream', name="mavgps_up")
        CFDSolver.addFunction('mavgps', 'viscous_plane', name="mavgps_plane")

        CFDSolver.addFunction('aavgps', 'downstream', name="aavgps_down")
        CFDSolver.addFunction('aavgps', 'upstream', name="aavgps_up")
        CFDSolver.addFunction('aavgps', 'viscous_plane', name="aavgps_plane")

        self.CFDSolver= CFDSolver

        self.arch = {'real':None,'complex':None}

        self.arch['real'] = os.environ.get("PETSC_ARCH")
        pdir = os.environ.get("PETSC_DIR")
        # Directories in the root of the petsc directory
        dirs = [o for o in os.listdir(pdir) if os.path.isdir(pdir+'/'+o)]
        # See if any one has 'complex' in it...basically we want to see if
        # an architecture has 'complex' in it which means it is (most
        # likely) a complex build
        carch = []
        for d in dirs:
            if 'complex' in d.lower():
                carch.append(d)
        if len(carch) > 0:
            # take the first one if there are multiple
            self.arch['complex'] = carch[0]
        else:
            self.arch['complex'] = None

        print(self.arch)

    def test_fwd(self):
        ap = copy.deepcopy(self.ap)

        ap.setBCVar('Pressure',  79326.7, 'outlet')
        ap.addDV('Pressure', familyGroup='outlet', name='outlet_pressure')
        key = 'outlet_pressure'
        xDvDot = {key:1.0}
        self.CFDSolver(ap)


        resDot, funcsDot, fDot, hfDot = self.CFDSolver.computeJacobianVectorProductFwd(
        xDvDot=xDvDot, residualDeriv=True, funcDeriv=True, fDeriv=True, hfDeriv=True, mode='AD')

        resDot_FD, funcsDot_FD, fDot_FD, hfDot_FD = self.CFDSolver.computeJacobianVectorProductFwd(
        xDvDot=xDvDot, residualDeriv=True, funcDeriv=True, fDeriv=True, hfDeriv=True, mode='FD', h=1e-0)
        
        np.testing.assert_allclose(resDot_FD,resDot, rtol=1e-5)
        for func in funcsDot:
            np.testing.assert_allclose(funcsDot_FD[func],funcsDot[func], rtol=5e-4)
        np.testing.assert_allclose(fDot_FD,fDot, rtol=1e-5)
        np.testing.assert_allclose(hfDot_FD,hfDot, rtol=1e-5)


    # def test_fwd_CS(self):
    #     # arch=self.arch['complex']
    #     import petsc4py
    #     petsc4py.init(arch='complex-debug')
    #     from petsc4py import PETSc
    #     from python.pyADflow_C import ADFLOW_C
        
    #     self.options['useanksolver'] = False
    #     CFDSolver = ADFLOW_C(options=self.options)
    #     intSurfFile = os.path.join(baseDir, '../input_files/integration_plane_viscous.fmt')

    #     CFDSolver.addIntegrationSurface(intSurfFile, 'viscous_plane')
    #     CFDSolver.finalizeUserIntegrationSurfaces()

    #     CFDSolver.addFamilyGroup('upstream',['inlet'])
    #     CFDSolver.addFamilyGroup('downstream',['outlet'])
    #     CFDSolver.addFamilyGroup('all_flow',['inlet', 'outlet'])
    #     CFDSolver.addFamilyGroup('output_fam',['all_flow', 'allWalls'])

    #     CFDSolver.addFunction('mdot', 'upstream', name="mdot_up")
    #     CFDSolver.addFunction('mdot', 'downstream', name="mdot_down")
    #     CFDSolver.addFunction('mdot', 'viscous_plane', name="mdot_plane")

    #     CFDSolver.addFunction('mavgptot', 'downstream', name="mavgptot_down")
    #     CFDSolver.addFunction('mavgptot', 'upstream', name="mavgptot_up")
    #     CFDSolver.addFunction('mavgptot', 'viscous_plane', name="mavgptot_plane")

    #     CFDSolver.addFunction('aavgptot', 'downstream', name="aavgptot_down")
    #     CFDSolver.addFunction('aavgptot', 'upstream', name="aavgptot_up")
    #     CFDSolver.addFunction('aavgptot', 'viscous_plane', name="aavgptot_plane")

    #     CFDSolver.addFunction('mavgttot', 'downstream', name="mavgttot_down")
    #     CFDSolver.addFunction('mavgttot', 'upstream', name="mavgttot_up")
    #     CFDSolver.addFunction('mavgttot', 'viscous_plane', name="mavgttot_plane")

    #     CFDSolver.addFunction('mavgps', 'downstream', name="mavgps_down")
    #     CFDSolver.addFunction('mavgps', 'upstream', name="mavgps_up")
    #     CFDSolver.addFunction('mavgps', 'viscous_plane', name="mavgps_plane")

    #     CFDSolver.addFunction('aavgps', 'downstream', name="aavgps_down")
    #     CFDSolver.addFunction('aavgps', 'upstream', name="aavgps_up")
    #     ap = copy.deepcopy(self.ap)

    #     ap.setBCVar('Pressure',  79326.7, 'outlet')
    #     ap.addDV('Pressure', familyGroup='outlet', name='outlet_pressure')
    #     key = 'outlet_pressure'
    #     xDvDot = {key:1.0}
    #     CFDSolver(ap, writeSolution=False)

    #     resDot, funcsDot, fDot, hfDot = CFDSolver.computeJacobianVectorProductFwd(
    #     xDvDot=xDvDot, residualDeriv=True, funcDeriv=True, fDeriv=True, hfDeriv=True, mode='CS')
    #     pp(funcsDot)
    #     del(CFDSolver)
    #     del(ADFLOW_C)

    #     # print(resDot)
    #     # print(funcsDot)
    #     # print(fDot)
    #     # print(hfDot)



    def test_fwd_dict(self):
        group = 'outlet'
        BCVar = 'Pressure'


        ap = copy.deepcopy(self.ap)
        self.CFDSolver.setAeroProblem(ap)

        BCData = self.CFDSolver.getBCData(groupNames=[group])
        ap.setBCVar(BCVar,  BCData[group][BCVar], group)
        ap.addDV(BCVar, familyGroup=group, name='outlet_pressure')
        
        xDvDot = {}


        self.CFDSolver(ap)
        for DV in ap.DVs:
            xDvDot = {DV:np.array([1.0])}

            resDot, funcsDot, fDot, hfDot = self.CFDSolver.computeJacobianVectorProductFwd(
            xDvDot=xDvDot, residualDeriv=True, funcDeriv=True, fDeriv=True, hfDeriv=True, mode='AD')

            resDot_FD, funcsDot_FD, fDot_FD, hfDot_FD = self.CFDSolver.computeJacobianVectorProductFwd(
            xDvDot=xDvDot, residualDeriv=True, funcDeriv=True, fDeriv=True, hfDeriv=True, mode='FD', h=1)
            
            np.testing.assert_allclose(resDot_FD,resDot, rtol=1e-5)
            for func in funcsDot:
                np.testing.assert_allclose(funcsDot_FD[func],funcsDot[func], rtol=5e-4)
            np.testing.assert_allclose(fDot_FD,fDot, atol=1e-8)
            np.testing.assert_allclose(hfDot_FD,hfDot, atol=1e-8)

    def test_fwd_dict_vs_float(self):

        ap = copy.deepcopy(self.ap)

        ap.setBCVar('Pressure',  79326.7, 'outlet')
        ap.addDV('Pressure', familyGroup='outlet', name='outlet_pressure')
 
        self.CFDSolver(ap)
        resDot, funcsDot, fDot, hfDot = self.CFDSolver.computeJacobianVectorProductFwd(
        xDvDot={'outlet_pressure':1.0}, residualDeriv=True, funcDeriv=True, fDeriv=True, hfDeriv=True, mode='AD')
        
        ap = copy.deepcopy(self.ap)
        group = 'outlet'
        BCVar = 'Pressure'

        BCData = self.CFDSolver.getBCData(groupNames=[group])
        ap.setBCVar(BCVar,  BCData[group][BCVar], group)
        ap.addDV(BCVar, familyGroup=group, name='outlet_pressure')
        
        xDvDot = {}

        for DV in ap.DVs:
            xDvDot[DV] = np.array([1.0])


        self.CFDSolver(ap)
        resDot_dict, funcsDot_dict, fDot_dict, hfDot_dict = self.CFDSolver.computeJacobianVectorProductFwd(
        xDvDot=xDvDot, residualDeriv=True, funcDeriv=True, fDeriv=True, hfDeriv=True, mode='AD')     

        np.testing.assert_array_equal(resDot_dict, resDot)
        for func in funcsDot:
            np.testing.assert_allclose(funcsDot_dict[func],funcsDot[func], rtol=5e-4)
        np.testing.assert_array_equal(fDot_dict, fDot)
        np.testing.assert_array_equal(hfDot_dict, hfDot)


    def test_bwd(self):
        ap = copy.deepcopy(self.ap)

        ap.setBCVar('Pressure',  79326.7, 'outlet')
        ap.addDV('Pressure', familyGroup='outlet', name='outlet_pressure')
        key = 'outlet_pressure'
        xDvDot = {key:1.0}
        self.CFDSolver(ap)


        resDot, _, _, __ = self.CFDSolver.computeJacobianVectorProductFwd(
        xDvDot=xDvDot, residualDeriv=True, funcDeriv=True, fDeriv=True, hfDeriv=True, mode='AD')

        dwBar = self.CFDSolver.getStatePerturbation(314)

        _, _, xDvBar = self.CFDSolver.computeJacobianVectorProductBwd(
            resBar=dwBar, wDeriv=True, xVDeriv=True, xDvDerivAero=True)


        np.testing.assert_array_almost_equal(xDvBar[key], np.dot(resDot, dwBar), decimal=14)
        


    def test_bwd_dict(self):
        group = 'outlet'
        BCVar = 'Pressure'


        ap = copy.deepcopy(self.ap)
        self.CFDSolver.setAeroProblem(ap)

        BCData = self.CFDSolver.getBCData(groupNames=[group])
        ap.setBCVar(BCVar,  BCData[group][BCVar], group)
        ap.addDV(BCVar, familyGroup=group, name='outlet_pressure')
        
        xDvDot = {}


        self.CFDSolver(ap)
        for ii, DV in enumerate(ap.DVs):
            xDvDot = {DV:np.array([1.0])}

            resDot, funcsDot, fDot, hfDot = self.CFDSolver.computeJacobianVectorProductFwd(
            xDvDot=xDvDot, residualDeriv=True, funcDeriv=True, fDeriv=True, hfDeriv=True, mode='AD')


            dwBar = self.CFDSolver.getStatePerturbation(ii)

            _, _, xDvBar = self.CFDSolver.computeJacobianVectorProductBwd(
                resBar=dwBar, wDeriv=True, xVDeriv=True, xDvDerivAero=True)
            np.testing.assert_array_almost_equal(xDvBar[DV], np.dot(resDot, dwBar), decimal=14)


class BCArrayDerivsTests(unittest.TestCase):
    def setUp(self):
        import petsc4py

        petsc4py.init(arch='real-debug-gfortran-3.7.7')
        from petsc4py import PETSc

        gridFile = os.path.join(baseDir, '../input_files/conic_conv_nozzle_mb_L4_array.cgns')
        intSurfFile = os.path.join(baseDir, '../input_files/integration_plane_viscous.fmt')

        self.options = {'gridfile': gridFile}

        self.options['mgstartlevel'] = 1

        self.options['ncycles'] = 1
        self.options['useanksolver'] = False
        self.options['l2convergence'] = 1e-12
        ap = AeroProblem(name='nozzle_flow', alpha=90, mach=0.5, altitude=0,
                    areaRef=1.0, chordRef=1.0, R=287.87,
                    evalFuncs=['mdot_up', 'mdot_down', 'mdot_plane',
                                'mavgptot_up', 'mavgptot_down', 'mavgptot_plane',
                                'mavgttot_up', 'mavgttot_down', 'mavgttot_plane',
                                'mavgps_up', 'mavgps_down', 'mavgps_plane',
                                ])

        self.ap = ap 

        CFDSolver = ADFLOW(options=self.options)
        CFDSolver.addIntegrationSurface(intSurfFile, 'viscous_plane')
        CFDSolver.finalizeUserIntegrationSurfaces()

        CFDSolver.addFamilyGroup('upstream',['inlet'])
        CFDSolver.addFamilyGroup('downstream',['outlet'])
        CFDSolver.addFamilyGroup('all_flow',['inlet', 'outlet'])
        CFDSolver.addFamilyGroup('output_fam',['all_flow', 'allWalls'])

        CFDSolver.addFunction('mdot', 'upstream', name="mdot_up")
        CFDSolver.addFunction('mdot', 'downstream', name="mdot_down")
        CFDSolver.addFunction('mdot', 'viscous_plane', name="mdot_plane")

        CFDSolver.addFunction('mavgptot', 'downstream', name="mavgptot_down")
        CFDSolver.addFunction('mavgptot', 'upstream', name="mavgptot_up")
        CFDSolver.addFunction('mavgptot', 'viscous_plane', name="mavgptot_plane")

        CFDSolver.addFunction('aavgptot', 'downstream', name="aavgptot_down")
        CFDSolver.addFunction('aavgptot', 'upstream', name="aavgptot_up")
        CFDSolver.addFunction('aavgptot', 'viscous_plane', name="aavgptot_plane")

        CFDSolver.addFunction('mavgttot', 'downstream', name="mavgttot_down")
        CFDSolver.addFunction('mavgttot', 'upstream', name="mavgttot_up")
        CFDSolver.addFunction('mavgttot', 'viscous_plane', name="mavgttot_plane")

        CFDSolver.addFunction('mavgps', 'downstream', name="mavgps_down")
        CFDSolver.addFunction('mavgps', 'upstream', name="mavgps_up")
        CFDSolver.addFunction('mavgps', 'viscous_plane', name="mavgps_plane")

        CFDSolver.addFunction('aavgps', 'downstream', name="aavgps_down")
        CFDSolver.addFunction('aavgps', 'upstream', name="aavgps_up")
        CFDSolver.addFunction('aavgps', 'viscous_plane', name="aavgps_plane")

        self.CFDSolver= CFDSolver

    def test_fwd(self):
        group = 'outlet'
        BCVar = 'Pressure'


        ap = copy.deepcopy(self.ap)
        self.CFDSolver.setAeroProblem(ap)

        BCData = self.CFDSolver.getBCData(groupNames=[group])
        ap.setBCVar(BCVar,  BCData[group][BCVar], group)
        ap.addDV(BCVar, familyGroup=group, name='outlet_pressure')
        

        self.CFDSolver(ap)
        for DV in ap.DVs:
            for ii in range(ap.DVs[DV].value.size):
                seed = np.zeros(ap.DVs[DV].value.size)
                seed[ii] = 1
                xDvDot = {DV:seed}
                resDot, funcsDot, fDot, hfDot = self.CFDSolver.computeJacobianVectorProductFwd(
                xDvDot=xDvDot, residualDeriv=True, funcDeriv=True, fDeriv=True, hfDeriv=True, mode='AD')

                resDot_FD, funcsDot_FD, fDot_FD, hfDot_FD = self.CFDSolver.computeJacobianVectorProductFwd(
                xDvDot=xDvDot, residualDeriv=True, funcDeriv=True, fDeriv=True, hfDeriv=True, mode='FD', h=1)
                np.testing.assert_allclose(resDot_FD,resDot, rtol=1e-5)
                for func in funcsDot:
                    np.testing.assert_allclose(funcsDot_FD[func],funcsDot[func], rtol=5e-4)
                np.testing.assert_allclose(fDot_FD,fDot, atol=1e-8)
                np.testing.assert_allclose(hfDot_FD,hfDot, atol=1e-8)

    def test_fwd_vs_float(self):

        ap = copy.deepcopy(self.ap)

        ap.setBCVar('Pressure',  79326.7, 'outlet')
        ap.addDV('Pressure', familyGroup='outlet', name='outlet_pressure')
 
        self.CFDSolver(ap)
        resDot, funcsDot, fDot, hfDot = self.CFDSolver.computeJacobianVectorProductFwd(
        xDvDot={'outlet_pressure':1.0}, residualDeriv=True, funcDeriv=True, fDeriv=True, hfDeriv=True, mode='AD')
        
        ap = copy.deepcopy(self.ap)
        group = 'outlet'
        BCVar = 'Pressure'

        BCData = self.CFDSolver.getBCData(groupNames=[group])
        ap.setBCVar(BCVar,  BCData[group][BCVar], group)
        ap.addDV(BCVar, familyGroup=group, name='outlet_pressure')
        
        xDvDot = {}

        for DV in ap.DVs:
            xDvDot[DV] =  np.ones(ap.DVs[DV].value.size)


        self.CFDSolver(ap)
        resDot_dict, funcsDot_dict, fDot_dict, hfDot_dict = self.CFDSolver.computeJacobianVectorProductFwd(
        xDvDot=xDvDot, residualDeriv=True, funcDeriv=True, fDeriv=True, hfDeriv=True, mode='AD')     

        np.testing.assert_array_equal(resDot_dict, resDot)
        for func in funcsDot:
            np.testing.assert_allclose(funcsDot_dict[func],funcsDot[func], rtol=5e-4)
        np.testing.assert_array_equal(fDot_dict, fDot)
        np.testing.assert_array_equal(hfDot_dict, hfDot)


    def test_bwd(self):

        group = 'outlet'
        BCVar = 'Pressure'


        ap = copy.deepcopy(self.ap)
        self.CFDSolver.setAeroProblem(ap)

        BCData = self.CFDSolver.getBCData(groupNames=[group])
        ap.setBCVar(BCVar,  BCData[group][BCVar], group)
        ap.addDV(BCVar, familyGroup=group, name='outlet_pressure')
        

        self.CFDSolver(ap)
        for DV in ap.DVs:
            for ii in range(ap.DVs[DV].value.size):
                seed = np.zeros(ap.DVs[DV].value.size)
                seed[ii] = 1
                xDvDot = {DV:seed}
                resDot, _, _, _ = self.CFDSolver.computeJacobianVectorProductFwd(
                xDvDot=xDvDot, residualDeriv=True, funcDeriv=True, fDeriv=True, hfDeriv=True, mode='AD')

                dwBar = self.CFDSolver.getStatePerturbation(ii)

                _, _, xDvBar = self.CFDSolver.computeJacobianVectorProductBwd(
                    resBar=dwBar, wDeriv=True, xVDeriv=True, xDvDerivAero=True)

                np.testing.assert_array_almost_equal(np.dot(xDvBar[DV],xDvDot[DV]),\
                                                     np.dot(resDot, dwBar), decimal=14)






if __name__ == '__main__':
    unittest.main()
