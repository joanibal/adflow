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
from adflow.pyADflow import ADFLOW
from mpi4py import MPI


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
                'TemperatureStagnation': 500,
                'VelocityUnitVectorX': 0.,
                'VelocityUnitVectorY': 1.,
                'VelocityUnitVectorZ': 0.,
            }
        }

        self.set_dvs = {
            'downstream': {
                'Pressure': 89326.7,
            },
            'upstream': {
                'PressureStagnation': 100001.0,
                'TemperatureStagnation': 450,
                'VelocityUnitVectorX': 0.1,
                'VelocityUnitVectorY': 0.9,
                'VelocityUnitVectorZ': 0.1,
            }
        }

    
    def test_getBCData(self):
        "Tests if mesh was read properly"
        CFDSolver = ADFLOW(options=self.options)
        CFDSolver.addFamilyGroup('upstream',['inlet'])
        CFDSolver.addFamilyGroup('downstream',['outlet'])

        bc_data = CFDSolver.getBCData(groupNames=self.initial_dvs.keys())

        ## check that the mesh bc match 
        for patch in bc_data:
            for arr in patch:
                np.testing.assert_array_equal(arr.data, np.array(self.initial_dvs[patch.fam][arr.name]))


    def test_setBCData(self):
        "Test if bcdata is correctly set to data arrays"
        
        


        CFDSolver = ADFLOW(options=self.options)
        CFDSolver.addFamilyGroup('upstream',['inlet'])
        CFDSolver.addFamilyGroup('downstream',['outlet'])
        bc_data = CFDSolver.getBCData(groupNames=self.set_dvs.keys())
        
        for patch in bc_data:
            if patch.fam in self.set_dvs: 
                for arr in patch:
                    if arr.name in self.set_dvs[patch.fam]:
                        arr.data = self.set_dvs[patch.fam][arr.name]
        
        CFDSolver.setBCData(bc_data)
        # CFDSolver.setAeroProblem(self.ap)

        new_bc_data = CFDSolver.getBCData(groupNames=self.set_dvs.keys())

        ## check that the mesh bc match 
        for patch in new_bc_data:
            for arr in patch:
                np.testing.assert_array_equal(arr.data,
                 np.array(self.set_dvs[patch.fam][arr.name])
                 , err_msg=arr.name)




    def test_getBCData_array(self):
        gridFile = os.path.join(baseDir, '../input_files/conic_conv_nozzle_mb_L4_array.cgns')
        options = copy.deepcopy(self.options)
        options['gridfile'] =  gridFile
        CFDSolver = ADFLOW(options=options)
        CFDSolver.addFamilyGroup('upstream',['inlet'])
        CFDSolver.addFamilyGroup('downstream',['outlet'])


        bc_data = CFDSolver.getBCData(groupNames=self.initial_dvs.keys())

        ## check that the mesh bc match 
        for patch in bc_data:
            for arr in patch:
                np.testing.assert_array_equal(arr.data, np.ones(arr.size)*
                                                self.initial_dvs[patch.fam][arr.name])

     


    def test_setBCData_array(self):
        gridFile = os.path.join(baseDir, '../input_files/conic_conv_nozzle_mb_L4_array.cgns')
        options = copy.deepcopy(self.options)
        options['gridfile'] =  gridFile
        CFDSolver = ADFLOW(options=options)
        CFDSolver.addFamilyGroup('upstream',['inlet'])
        CFDSolver.addFamilyGroup('downstream',['outlet'])

        bc_data = CFDSolver.getBCData(groupNames=self.set_dvs.keys())
        for patch in bc_data:
            if patch.fam in self.set_dvs: 
                for arr in patch:
                    if arr.name in self.set_dvs[patch.fam]:
                        arr.data = self.set_dvs[patch.fam][arr.name]*np.ones(arr.size)
        
        CFDSolver.setBCData(bc_data)
        bc_data = CFDSolver.getBCData(groupNames=self.set_dvs.keys())

        for patch in bc_data:
            for arr in patch:
                np.testing.assert_array_equal(arr.data, np.ones(arr.size)*
                                                self.set_dvs[patch.fam][arr.name])

     
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


    def test_ap_bcvar_scalar(self):
        CFDSolver = ADFLOW(options=self.options)
        CFDSolver.addFamilyGroup('upstream',['inlet'])
        CFDSolver.addFamilyGroup('downstream',['outlet'])

        ap = self.ap
        group = 'outlet'
        BCVar = 'Pressure'
        # ap.bc_data = bc_data

        ap.setBCVar(BCVar,  70000.0, group)
        ap.addDV(BCVar,  name='outlet_pressure',  familyGroup=group,)
        
        DVs = {'outlet_pressure':123000.0}

        ap.setDesignVars(DVs)

        assert ap.getBCData()[group][BCVar] == DVs['outlet_pressure']

        CFDSolver.setAeroProblem(ap)
        bc_data = CFDSolver.getBCData(groupNames=[group])
        
        for d in bc_data.getBCArraysFlatData(BCVar, familyGroup=group):
            assert d == DVs['outlet_pressure']



    def test_ap_bcvar_array(self):
        CFDSolver = ADFLOW(options=self.options)
        CFDSolver.addFamilyGroup('upstream',['inlet'])
        CFDSolver.addFamilyGroup('downstream',['outlet'])

        ap = self.ap
        group = 'outlet'
        BCVar = 'Pressure'
        # ap.bc_data = bc_data

        ap.setBCVar(BCVar,  70000.0, group)
        ap.addDV(BCVar,  name='outlet_pressure',  familyGroup=group,)
        
        # DVs = {'outlet_pressure':123000.0}
        DVs = {'outlet_pressure':np.arange(1,6)*10**5}

        ap.setDesignVars(DVs)
        print(ap.getBCData()[group][BCVar] == DVs['outlet_pressure'])
        assert (ap.getBCData()[group][BCVar] == DVs['outlet_pressure']).all()

        CFDSolver.setAeroProblem(ap)
        bc_data = CFDSolver.getBCData(groupNames=[group])
        for ii, d in enumerate(bc_data.getBCArraysFlatData(BCVar, familyGroup=group)):
            assert d == DVs['outlet_pressure'][ii]



    def test_ap_bcvar_nodal_array(self):
        gridFile = os.path.join(baseDir, '../input_files/conic_conv_nozzle_mb_L4_array.cgns')
        options = copy.deepcopy(self.options)
        options['gridfile'] =  gridFile
        CFDSolver = ADFLOW(options=options)
        CFDSolver.addFamilyGroup('upstream',['inlet'])
        CFDSolver.addFamilyGroup('downstream',['outlet'])

  
        ap = self.ap
        group = 'outlet'
        BCVar = 'Pressure'
        # ap.bc_data = bc_data

        ap.setBCVar(BCVar,  70000.0, group)
        ap.addDV(BCVar,  name='outlet_pressure',  familyGroup=group,)
        
        # DVs = {'outlet_pressure':123000.0}
        DVs = {'outlet_pressure':np.arange(1,5*9+1)*10**5}

        ap.setDesignVars(DVs)
        assert (ap.getBCData()[group][BCVar] == DVs['outlet_pressure']).all()

        CFDSolver.setAeroProblem(ap)
        bc_data = CFDSolver.getBCData(groupNames=[group])
        for ii, d in enumerate(bc_data.getBCArraysFlatData(BCVar, familyGroup=group)):
            assert d == DVs['outlet_pressure'][ii]

class BCDerivsTests(unittest.TestCase):
    def setUp(self):
        # import petsc4py

        # petsc4py.init(arch='real-debug-gfortran-3.7.7')
        # from petsc4py import PETSc

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

        # self.arch = {'real':None,'complex':None}

        # self.arch['real'] = os.environ.get("PETSC_ARCH")
        # pdir = os.environ.get("PETSC_DIR")
        # # Directories in the root of the petsc directory
        # dirs = [o for o in os.listdir(pdir) if os.path.isdir(pdir+'/'+o)]
        # # See if any one has 'complex' in it...basically we want to see if
        # # an architecture has 'complex' in it which means it is (most
        # # likely) a complex build
        # carch = []
        # for d in dirs:
        #     if 'complex' in d.lower():
        #         carch.append(d)
        # if len(carch) > 0:
        #     # take the first one if there are multiple
        #     self.arch['complex'] = carch[0]
        # else:
        #     self.arch['complex'] = None

        # print(self.arch)

    def test_fwd(self):
        ap = copy.deepcopy(self.ap)
        group = 'outlet'
        BCVar = 'Pressure'

        ap.setBCVar('Pressure',  79326.7, 'outlet')
        ap.addDV(BCVar,  name='outlet_pressure',  familyGroup=group,)
        key = 'outlet_pressure'
        xDvDot = {key:1.0}
        self.CFDSolver(ap)


        resDot, funcsDot, fDot, hfDot = self.CFDSolver.computeJacobianVectorProductFwd(
        xDvDot=xDvDot, residualDeriv=True, funcDeriv=True, fDeriv=True, hfDeriv=True, mode='AD')

        resDot_FD, funcsDot_FD, fDot_FD, hfDot_FD = self.CFDSolver.computeJacobianVectorProductFwd(
        xDvDot=xDvDot, residualDeriv=True, funcDeriv=True, fDeriv=True, hfDeriv=True, mode='FD', h=1e-0)
        np.testing.assert_allclose(resDot_FD,resDot, atol=1e-8)
        # atol used here because the FD can return an error order 1e-15, but the rel error is nan since it should be 0
        for func in funcsDot:
            np.testing.assert_allclose(funcsDot_FD[func],funcsDot[func],  atol=1e-5)
        
        np.testing.assert_allclose(fDot_FD,fDot, atol=1e-10)
        np.testing.assert_allclose(hfDot_FD,hfDot, atol=1e-10)


#     # def test_fwd_CS(self):
#     #     # arch=self.arch['complex']
#     #     import petsc4py
#     #     petsc4py.init(arch='complex-debug')
#     #     from petsc4py import PETSc
#     #     from python.pyADflow_C import ADFLOW_C
        
#     #     self.options['useanksolver'] = False
#     #     CFDSolver = ADFLOW_C(options=self.options)
#     #     intSurfFile = os.path.join(baseDir, '../input_files/integration_plane_viscous.fmt')

#     #     CFDSolver.addIntegrationSurface(intSurfFile, 'viscous_plane')
#     #     CFDSolver.finalizeUserIntegrationSurfaces()

#     #     CFDSolver.addFamilyGroup('upstream',['inlet'])
#     #     CFDSolver.addFamilyGroup('downstream',['outlet'])
#     #     CFDSolver.addFamilyGroup('all_flow',['inlet', 'outlet'])
#     #     CFDSolver.addFamilyGroup('output_fam',['all_flow', 'allWalls'])

#     #     CFDSolver.addFunction('mdot', 'upstream', name="mdot_up")
#     #     CFDSolver.addFunction('mdot', 'downstream', name="mdot_down")
#     #     CFDSolver.addFunction('mdot', 'viscous_plane', name="mdot_plane")

#     #     CFDSolver.addFunction('mavgptot', 'downstream', name="mavgptot_down")
#     #     CFDSolver.addFunction('mavgptot', 'upstream', name="mavgptot_up")
#     #     CFDSolver.addFunction('mavgptot', 'viscous_plane', name="mavgptot_plane")

#     #     CFDSolver.addFunction('aavgptot', 'downstream', name="aavgptot_down")
#     #     CFDSolver.addFunction('aavgptot', 'upstream', name="aavgptot_up")
#     #     CFDSolver.addFunction('aavgptot', 'viscous_plane', name="aavgptot_plane")

#     #     CFDSolver.addFunction('mavgttot', 'downstream', name="mavgttot_down")
#     #     CFDSolver.addFunction('mavgttot', 'upstream', name="mavgttot_up")
#     #     CFDSolver.addFunction('mavgttot', 'viscous_plane', name="mavgttot_plane")

#     #     CFDSolver.addFunction('mavgps', 'downstream', name="mavgps_down")
#     #     CFDSolver.addFunction('mavgps', 'upstream', name="mavgps_up")
#     #     CFDSolver.addFunction('mavgps', 'viscous_plane', name="mavgps_plane")

#     #     CFDSolver.addFunction('aavgps', 'downstream', name="aavgps_down")
#     #     CFDSolver.addFunction('aavgps', 'upstream', name="aavgps_up")
#     #     ap = copy.deepcopy(self.ap)

#     #     ap.setBCVar('Pressure',  79326.7, 'outlet')
#     #     ap.addDV('Pressure', familyGroup='outlet', name='outlet_pressure')
#     #     key = 'outlet_pressure'
#     #     xDvDot = {key:1.0}
#     #     CFDSolver(ap, writeSolution=False)

#     #     resDot, funcsDot, fDot, hfDot = CFDSolver.computeJacobianVectorProductFwd(
#     #     xDvDot=xDvDot, residualDeriv=True, funcDeriv=True, fDeriv=True, hfDeriv=True, mode='CS')
#     #     pp(funcsDot)
#     #     del(CFDSolver)
#     #     del(ADFLOW_C)

#     #     # print(resDot)
#     #     # print(funcsDot)
#     #     # print(fDot)
#     #     # print(hfDot)



    def test_fwd_array(self):



  
        ap = self.ap
        self.CFDSolver.setAeroProblem(ap)
        group = 'outlet'
        BCVar = 'Pressure'
        # ap.bc_data = bc_data



        ap.setBCVar(BCVar,  79326.7, group)
        ap.addDV(BCVar,  name='outlet_pressure',  familyGroup=group,)
        
        # DVs = {'outlet_pressure':123000.0}
        # DVs = {'outlet_pressure':np.arange(1,5*9+1)*10**5}

        xDvDot = {}


        xDvDot = {'outlet_pressure':np.zeros(5)}
        self.CFDSolver(ap)
        for ii in range(5):


            xDvDot['outlet_pressure'][ii] = 1

            resDot, funcsDot, fDot, hfDot = self.CFDSolver.computeJacobianVectorProductFwd(
            xDvDot=xDvDot, residualDeriv=True, funcDeriv=True, fDeriv=True, hfDeriv=True, mode='AD')

            resDot_FD, funcsDot_FD, fDot_FD, hfDot_FD = self.CFDSolver.computeJacobianVectorProductFwd(
            xDvDot=xDvDot, residualDeriv=True, funcDeriv=True, fDeriv=True, hfDeriv=True, mode='FD', h=1)
            
            np.testing.assert_allclose(resDot_FD,resDot, atol=1e-8)
            for func in funcsDot:
                np.testing.assert_allclose(funcsDot_FD[func],funcsDot[func],  atol=1e-5)
            
            np.testing.assert_allclose(fDot_FD,fDot, atol=1e-10)
            np.testing.assert_allclose(hfDot_FD,hfDot, atol=1e-10)

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

        resDot_dict, funcsDot_dict, fDot_dict, hfDot_dict = self.CFDSolver.computeJacobianVectorProductFwd(
        xDvDot={'outlet_pressure': np.ones(5)}, residualDeriv=True, funcDeriv=True, fDeriv=True, hfDeriv=True, mode='AD')     

        np.testing.assert_allclose(resDot_dict,resDot, atol=1e-20)

        for func in funcsDot:
            np.testing.assert_allclose(funcsDot_dict[func],funcsDot[func],  atol=1e-20)

        np.testing.assert_allclose(fDot_dict,fDot, atol=1e-20)
        np.testing.assert_allclose(hfDot_dict,hfDot, atol=1e-20)


            
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


        np.testing.assert_array_almost_equal(np.sum(xDvBar[key]), np.dot(resDot, dwBar), decimal=14)
        


    def test_bwd_array(self):

  
        ap = self.ap
        self.CFDSolver.setAeroProblem(ap)
        group = 'outlet'
        BCVar = 'Pressure'
        # ap.bc_data = bc_data
        DV = 'outlet_pressure'


        ap.setBCVar(BCVar,  np.ones(5)*79326.7, group)
        ap.addDV(BCVar,  name=DV,  familyGroup=group,)
        
        # DVs = {'outlet_pressure':123000.0}
        # DVs = {'outlet_pressure':np.arange(1,5*9+1)*10**5}

        xDvDot = {}


        xDvDot = {DV:np.zeros(5)}
        self.CFDSolver(ap)
        for ii in range(5):


            xDvDot[DV][ii] = 1

            resDot, funcsDot, fDot, hfDot = self.CFDSolver.computeJacobianVectorProductFwd(
            xDvDot=xDvDot, residualDeriv=True, funcDeriv=True, fDeriv=True, hfDeriv=True, mode='AD')


            dwBar = self.CFDSolver.getStatePerturbation(ii)

            _, _, xDvBar = self.CFDSolver.computeJacobianVectorProductBwd(
                resBar=dwBar, wDeriv=True, xVDeriv=True, xDvDerivAero=True)
            
            print(np.dot(xDvBar[DV], xDvDot[DV]), np.dot(resDot, dwBar))
            np.testing.assert_array_almost_equal(np.dot(xDvBar[DV], xDvDot[DV]), np.dot(resDot, dwBar), decimal=14)

class BCNodalDerivsTests_serial(unittest.TestCase):
    N_PROCS = 1
    def setUp(self):

        gridFile = os.path.join(baseDir, '../input_files/conic_conv_nozzle_mb_L4_array.cgns')
        intSurfFile = os.path.join(baseDir, '../input_files/integration_plane_viscous.fmt')

        self.options = {'gridfile': gridFile}

        self.options['mgstartlevel'] = 1
        self.options['mgstartlevel'] = 1
        self.options['mgcycle'] = 'sg'
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
        DV = 'outlet_pressure'
        n =  5*9
        ap = self.ap
        self.CFDSolver.setAeroProblem(self.ap)

        ap.setBCVar(BCVar,  79326.7, group)
        ap.addDV(BCVar,  name=DV,  familyGroup=group)
        

        self.CFDSolver(ap)
        for ii in range(n):
            seed = np.zeros(n)
            seed[ii] = 1
            xDvDot = {DV:seed}

            resDot, funcsDot, fDot, hfDot = self.CFDSolver.computeJacobianVectorProductFwd(
            xDvDot=xDvDot, residualDeriv=True, funcDeriv=True, fDeriv=True, hfDeriv=True, mode='AD')

            resDot_FD, funcsDot_FD, fDot_FD, hfDot_FD = self.CFDSolver.computeJacobianVectorProductFwd(
            xDvDot=xDvDot, residualDeriv=True, funcDeriv=True, fDeriv=True, hfDeriv=True, mode='FD', h=1)
            
            np.testing.assert_allclose(resDot_FD,resDot, atol=1e-8)
            # atol used here because the FD can return an error order 1e-15, but the rel error is nan since it should be 0
            for func in funcsDot:
                np.testing.assert_allclose(funcsDot_FD[func],funcsDot[func],  atol=1e-5)
            
            np.testing.assert_allclose(fDot_FD,fDot, atol=1e-10)
            np.testing.assert_allclose(hfDot_FD,hfDot, atol=1e-10)


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

        resDot_dict, funcsDot_dict, fDot_dict, hfDot_dict = self.CFDSolver.computeJacobianVectorProductFwd(
        xDvDot={'outlet_pressure': np.ones(5*9)}, residualDeriv=True, funcDeriv=True, fDeriv=True, hfDeriv=True, mode='AD')     

        np.testing.assert_allclose(resDot_dict,resDot, atol=1e-20)

        for func in funcsDot:
            np.testing.assert_allclose(funcsDot_dict[func],funcsDot[func],  atol=1e-20)

        np.testing.assert_allclose(fDot_dict,fDot, atol=1e-20)
        np.testing.assert_allclose(hfDot_dict,hfDot, atol=1e-20)


    def test_bwd(self):

        group = 'outlet'
        BCVar = 'Pressure'
        DV = 'outlet_pressure'
        n =  5*9
        ap = self.ap
        self.CFDSolver.setAeroProblem(self.ap)
        
        bc_data = self.CFDSolver.getBCData()
        bc_array = bc_data.getBCArraysFlatData(BCVar, familyGroup=group)
        
        ap.setBCVar(BCVar,  np.ones(bc_array.size)*79326.7, group)
        ap.addDV(BCVar,  name=DV,  familyGroup=group)
        

        self.CFDSolver(ap)
        for ii in range(n):
            seed = np.zeros(n)
            seed[ii] = 1
            xDvDot = {DV:seed}

            resDot, funcsDot, fDot, hfDot = self.CFDSolver.computeJacobianVectorProductFwd(
            xDvDot=xDvDot, residualDeriv=True, funcDeriv=True, fDeriv=True, hfDeriv=True, mode='AD')

            dwBar = self.CFDSolver.getStatePerturbation(ii)

            _, _, xDvBar = self.CFDSolver.computeJacobianVectorProductBwd(
                resBar=dwBar, wDeriv=True, xVDeriv=True, xDvDerivAero=True)

            np.testing.assert_array_almost_equal(np.dot(xDvBar[DV],xDvDot[DV]),\
                                                    np.dot(resDot, dwBar), decimal=14)


class BCNodalDerivsTests_parallel(unittest.TestCase):
    N_PROCS = 2

    def setUp(self):
        # import petsc4py
        # petsc4py.init(arch='real-debug-gfortran-3.7.7')
        # from petsc4py import PETSc

        gridFile = os.path.join(baseDir, '../input_files/conic_conv_nozzle_mb_L4_array.cgns')
        intSurfFile = os.path.join(baseDir, '../input_files/integration_plane_viscous.fmt')

        self.options = {'gridfile': gridFile}

        self.options['mgstartlevel'] = 1
        self.options['mgstartlevel'] = 1
        self.options['mgcycle'] = 'sg'
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

        self.comm = MPI.COMM_WORLD

    def test_fwd(self):
        group = 'outlet'
        BCVar = 'Pressure'
        DV = 'outlet_pressure'
        n =  5*9
        ap = self.ap
        self.CFDSolver.setAeroProblem(self.ap)
        
        bc_data = self.CFDSolver.getBCData()
        bc_array = bc_data.getBCArraysFlatData(BCVar, familyGroup=group)
        n =  bc_array.size
        print('array size', bc_array.size)
        ap.setBCVar(BCVar,  np.ones(bc_array.size)*79326.7, group)
        ap.addDV(BCVar,  name=DV,  familyGroup=group)
        

        self.CFDSolver(ap)

        for irank in range(self.comm.size):
            nnodes = bc_array.size
            nnodes = self.comm.bcast(nnodes, root=irank)
            for ii in range(nnodes):
                seed = np.zeros(n)
                if self.comm.rank == irank:
                    seed[ii] = 1
                
                xDvDot = {DV:seed}
                    

                resDot, funcsDot, fDot, hfDot = self.CFDSolver.computeJacobianVectorProductFwd(
                xDvDot=xDvDot, residualDeriv=True, funcDeriv=True, fDeriv=True, hfDeriv=True, mode='AD')

                resDot_FD, funcsDot_FD, fDot_FD, hfDot_FD = self.CFDSolver.computeJacobianVectorProductFwd(
                xDvDot=xDvDot, residualDeriv=True, funcDeriv=True, fDeriv=True, hfDeriv=True, mode='FD', h=1)
                
                np.testing.assert_allclose(resDot_FD,resDot, atol=1e-8)
                # atol used here because the FD can return an error order 1e-15, but the rel error is nan since it should be 0
                for func in funcsDot:
                    np.testing.assert_allclose(funcsDot_FD[func],funcsDot[func],  atol=1e-5)
                
                np.testing.assert_allclose(fDot_FD,fDot, atol=1e-10)
                np.testing.assert_allclose(hfDot_FD,hfDot, atol=1e-10)

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

        bc_data = self.CFDSolver.getBCData()
        bc_array = bc_data.getBCArraysFlatData(BCVar, familyGroup=group)
        n =  bc_array.size

        resDot_dict, funcsDot_dict, fDot_dict, hfDot_dict = self.CFDSolver.computeJacobianVectorProductFwd(
        xDvDot={'outlet_pressure': np.ones(n)}, residualDeriv=True, funcDeriv=True, fDeriv=True, hfDeriv=True, mode='AD')     

        np.testing.assert_allclose(resDot_dict,resDot, atol=1e-20)

        for func in funcsDot:
            np.testing.assert_allclose(funcsDot_dict[func],funcsDot[func],  atol=1e-20)

        np.testing.assert_allclose(fDot_dict,fDot, atol=1e-20)
        np.testing.assert_allclose(hfDot_dict,hfDot, atol=1e-20)

    def test_bwd(self):

        group = 'outlet'
        BCVar = 'Pressure'
        DV = 'outlet_pressure'
        ap = self.ap
        self.CFDSolver.setAeroProblem(self.ap)
        
        bc_data = self.CFDSolver.getBCData()
        bc_array = bc_data.getBCArraysFlatData(BCVar, familyGroup=group)
        n =  bc_array.size
        print('array size', bc_array.size)
        ap.setBCVar(BCVar,  np.ones(bc_array.size)*79326.7, group)
        ap.addDV(BCVar,  name=DV,  familyGroup=group)
        

        self.CFDSolver(ap)

        for irank in range(self.comm.size):
            nnodes = bc_array.size
            nnodes = self.comm.bcast(nnodes, root=irank)
            for ii in range(nnodes):
                seed = np.zeros(n)
                if self.comm.rank == irank:
                    seed[ii] = 1
                
                xDvDot = {DV:seed}
                    


                resDot, funcsDot, fDot, hfDot = self.CFDSolver.computeJacobianVectorProductFwd(
                xDvDot=xDvDot, residualDeriv=True, funcDeriv=True, fDeriv=True, hfDeriv=True, mode='AD')

                dwBar = self.CFDSolver.getStatePerturbation(ii)

                _, _, xDvBar = self.CFDSolver.computeJacobianVectorProductBwd(
                    resBar=dwBar, wDeriv=True, xVDeriv=True, xDvDerivAero=True)

                np.testing.assert_array_almost_equal(np.dot(xDvBar[DV],xDvDot[DV]),\
                                                        np.dot(resDot, dwBar), decimal=14)



if __name__ == '__main__':
    unittest.main()
