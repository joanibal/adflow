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



class BaseTest():
    "class to hold the common setup method for the test"

    def setUp(self):
        gridFile = os.path.join(baseDir, '../input_files/n0012_hot_000_vol.cgns')
        aero_options = {
            # I/O Parameters
            'gridfile': gridFile,
            "restartfile": gridFile,
            "equationType": "RANS",
            'monitorvariables':['cpu', 'resrho', 'resturb', 'cl', 'cd', "totheattransfer"],
            "writevolumesolution": True,
            'writesurfacesolution':False,


            "MGCycle": "sg",
            "MGStartLevel": -1,
            "useANKSolver": True,
            "useNKSolver": True,

            'L2Convergence':1e-14,
            'nCycles':4000,
            'adjointl2convergence': 1e-14,
            'useblockettes':False
        }


        # Aerodynamic problem description (taken from test 17)
        self.ap = AeroProblem(
            name="n0012_hot",
            V=32,  # m/s
            T=273 + 60,  # kelvin
            P=93e3,  # pa
            areaRef=1.0,  # m^2
            chordRef=1.0,  # m^2
            evalFuncs=["cd", "totheattransfer", "havg", "area", "hot_area"],
            alpha=0.0,
            beta=0.00,
            xRef=0.0,
            yRef=0.0,
            zRef=0.0,
        )

        # these values come from inspecting the cgns mesh itself
        self.CFDSolver = ADFLOW(options=aero_options)
        self.CFDSolver.addFunction('area', 'heated_wall', name="hot_area")

        self.CFDSolver.getResidual(self.ap)
        # self.CFDSolver(self.ap)


# @unittest.skip("")
class HeatXferTests(BaseTest, unittest.TestCase):

    def test_bcdata(self):
        "Check that the temperature data set on the surface is read correctly"
        bc_data = self.CFDSolver.getBCData(['heated_wall'])
        print(bc_data)
        for patch in bc_data:
            for arr in patch:
                np.testing.assert_array_equal(arr.data, np.ones(arr.size)*273)

    def test_havg(self):
        "Check that the func havg (average heattransfer coefficent) is calculated right"
        bc_data = self.CFDSolver.getBCData(['heated_wall'])


        temp = bc_data.getBCArraysData('Temperature')

        self.CFDSolver(self.ap)

        funcs = {}
        self.CFDSolver.evalFunctions(self.ap, funcs)
        h_avg = funcs[self.ap.name + '_havg']
        q = funcs[self.ap.name + '_totheattransfer']
        A = funcs[self.ap.name + '_hot_area']
        dt =  self.ap.T - np.mean(temp)


        np.testing.assert_allclose(q/(A*dt), h_avg,  rtol=1e-7)



    def test_heatxfer(self):
        "Check that the func totheattransfer is calculated right"

        self.CFDSolver(self.ap)

        funcs = {}
        self.CFDSolver.evalFunctions(self.ap, funcs)

        Q_dot  = funcs[self.ap.name + '_totheattransfer']
        q = self.CFDSolver.getHeatXferRates(groupName='heated_wall')


        np.testing.assert_allclose(Q_dot, np.sum(q),  rtol=1e-15)


    def test_heatfluxes(self):
        "Check that the heat fluxes are calculated right"
        bc_data = self.CFDSolver.getBCData(['heated_wall'])
        temp = bc_data.getBCArraysData('Temperature')


        dt =  self.ap.T - np.mean(temp)

        self.CFDSolver(self.ap)

        funcs = {}
        self.CFDSolver.evalFunctions(self.ap, funcs)
        A = funcs[self.ap.name + '_hot_area']

        Q_dot  = funcs[self.ap.name + '_totheattransfer']
        q = self.CFDSolver.getHeatFluxes(groupName='heated_wall')

        # this one is a regression test, becuase I could think of something to compare with
        np.testing.assert_allclose(20378.471567736087, np.sum(q),  rtol=5e-14)



#region derivative tests

# @unittest.skip("")
class HeatXferBCDVDerivsTests(BaseTest, unittest.TestCase):

    def test_fwd_BCDV(self):
        "test the FWD derivatives of the functional"


        group = 'heated_wall'
        BCVar = 'Temperature'


        self.ap.setBCVar(BCVar,  300.0, group)
        self.ap.addDV(BCVar,  name='wall_temp',  familyGroup=group)
        # self.CFDSolver(self.ap)
        self.CFDSolver.getResidual(self.ap)

        xDvDot = {"wall_temp":1.0}

        resDot_FD, funcsDot_FD, fDot_FD, hfDot_FD = self.CFDSolver.computeJacobianVectorProductFwd(
        xDvDot=xDvDot, residualDeriv=True, funcDeriv=True, fDeriv=True, hfDeriv=True, mode='FD', h=1e-6)

        resDot, funcsDot, fDot, hfDot = self.CFDSolver.computeJacobianVectorProductFwd(
        xDvDot=xDvDot, residualDeriv=True, funcDeriv=True, fDeriv=True, hfDeriv=True, mode='AD')


        np.testing.assert_allclose(resDot_FD,resDot, atol=5e-9)


        for func in funcsDot:
            np.testing.assert_allclose(funcsDot_FD[func],funcsDot[func], err_msg=func,  atol=1e-5)

        np.testing.assert_allclose(fDot_FD,fDot, atol=1e-10)

        np.testing.assert_allclose(hfDot_FD,hfDot, atol=1e-10)

    @unittest.skip("")
    def test_fwd_BCDV_CS(self):
        "test the FWD derivatives of the functional"

        from adflow.pyADflow_C import ADFLOW_C
        self.CFDSolver = ADFLOW_C(options=self.aero_options)


        group = 'heated_wall'
        BCVar = 'Temperature'


        self.ap.setBCVar(BCVar,  300.0, group)
        self.ap.addDV(BCVar,  name='wall_temp',  familyGroup=group)
        self.CFDSolver.getResidual(self.ap)

        xDvDot = {"wall_temp":1.0}

        # resDot, funcsDot, fDot, hfDot = self.CFDSolver.computeJacobianVectorProductFwd(
        # xDvDot=xDvDot, residualDeriv=True, funcDeriv=True, fDeriv=True, hfDeriv=True, mode='AD')

        resDot_CS, funcsDot_CS, fDot_CS, hfDot_CS = self.CFDSolver.computeJacobianVectorProductFwd(
        xDvDot=xDvDot, residualDeriv=True, funcDeriv=True, fDeriv=True, hfDeriv=True, mode='CS', h=1e-200)


    def test_bwd_BCDV(self):
        "test the FWD derivatives of the functional"


        group = 'heated_wall'
        BCVar = 'Temperature'
        func = 'wall_temp'

        self.ap.setBCVar(BCVar,  300.0, group)
        self.ap.addDV(BCVar,  name=func,  familyGroup=group)
        # self.CFDSolver(self.ap)
        self.CFDSolver.getResidual(self.ap)


        xDvDot = {func:1.0}

        resDot, funcsDot, fDot, hfDot = self.CFDSolver.computeJacobianVectorProductFwd(
        xDvDot=xDvDot, residualDeriv=True, funcDeriv=True, fDeriv=True, hfDeriv=True, mode='AD')



        dwBar = self.CFDSolver.getStatePerturbation(314)
        fBar = self.CFDSolver.getSurfacePerturbation(314)
        hxferBar = fBar[:,0]


        funcsBar = {}
        for key in self.ap.evalFuncs:
            funcsBar[key] =  1.0




        xVBar = self.CFDSolver.computeJacobianVectorProductBwd(
            resBar=dwBar, fBar=fBar, hfBar=hxferBar, funcsBar=funcsBar ,  xDvDeriv=True)

        # we have to add up all the parts
        _sum = 0
        _sum += np.dot(resDot, dwBar)
        _sum += np.dot(fDot.flatten(), fBar.flatten())
        _sum += np.dot(hfDot.flatten(), hxferBar.flatten())
        _sum += np.sum([x for x in funcsDot.values()])


        np.testing.assert_array_almost_equal(np.dot(xDvDot[func], xVBar[func]), _sum, decimal=14)

# @unittest.skip("")
class HeatXferXVDerivsTests(BaseTest, unittest.TestCase):

    # @unittest.skip("")
    def test_fwd_XV(self):
        "test the FWD derivatives wrt to nodal displacements"

        xVDot = self.CFDSolver.getSpatialPerturbation(321)


        resDot_FD, funcsDot_FD, fDot_FD, hfDot_FD = self.CFDSolver.computeJacobianVectorProductFwd(
        xVDot=xVDot, residualDeriv=True, funcDeriv=True, fDeriv=True, hfDeriv=True, mode='FD', h=1e-10)


        resDot, funcsDot, fDot, hfDot = self.CFDSolver.computeJacobianVectorProductFwd(
        xVDot=xVDot, residualDeriv=True, funcDeriv=True, fDeriv=True, hfDeriv=True, mode='AD')

        # print(resDot_FD[0:10])
        print(resDot[3446])
        print(funcsDot)
        # print(resDot[0:10])


        # these are checked against CS becuase the derivates appear to be poorly conditioned
        import pickle
        with open('resDot_CS.p', 'rb') as f:
            resDot_CS = pickle.load(f)

        np.testing.assert_allclose(resDot_CS, resDot, rtol=5e-8)


        for func in funcsDot:
            np.testing.assert_allclose(funcsDot_FD[func],funcsDot[func], err_msg=func,  rtol=1e-5)

        np.testing.assert_allclose(fDot_FD,fDot, rtol=5e-5)

        np.testing.assert_allclose(hfDot_FD[hfDot_FD != 0],hfDot[hfDot_FD != 0], rtol=5e-4)

    @unittest.skip("")
    def test_fwd_XV_CS(self):
        "test the FWD derivatives of the functional"

        from adflow.pyADflow_C import ADFLOW_C
        self.CFDSolver = ADFLOW_C(options=self.aero_options)
        self.CFDSolver.getResidual(self.ap)
        # self.CFDSolver(self.ap)

        xVDot = self.CFDSolver.getSpatialPerturbation(321)

        resDot_CS, funcsDot_CS, fDot_CS, hfDot_CS = self.CFDSolver.computeJacobianVectorProductFwd(
        xVDot=xVDot, residualDeriv=True, funcDeriv=True, fDeriv=True, hfDeriv=True, mode='CS', h=1e-200)

    # @unittest.skip("")
    def test_bwd_XV(self):
        "test the FWD derivatives of the functional"



        xVDot = self.CFDSolver.getSpatialPerturbation(321)

        resDot, funcsDot, fDot, hfDot = self.CFDSolver.computeJacobianVectorProductFwd(
        xVDot=xVDot, residualDeriv=True, funcDeriv=True, fDeriv=True, hfDeriv=True, mode='AD')



        dwBar = self.CFDSolver.getStatePerturbation(314)
        fBar = self.CFDSolver.getSurfacePerturbation(314)
        hxferBar = fBar[:,0]


        funcsBar = {}
        for key in self.ap.evalFuncs:
            funcsBar[key] =  1.0


        xVBar = self.CFDSolver.computeJacobianVectorProductBwd(
            resBar=dwBar, fBar=fBar, hfBar=hxferBar, funcsBar=funcsBar ,  xVDeriv=True)

        # we have to add up all the parts
        _sum = 0
        _sum += np.dot(resDot, dwBar)
        _sum += np.dot(fDot.flatten(), fBar.flatten())
        _sum += np.dot(hfDot.flatten(), hxferBar.flatten())
        _sum += np.sum([x for x in funcsDot.values()])


        np.testing.assert_allclose(np.dot(xVDot, xVBar), _sum, rtol=1e-14)



class HeatXferWDerivsTests(BaseTest, unittest.TestCase):

    # @unittest.skip("")
    def test_fwd_W(self):
        "test the FWD derivatives wrt to the states"

        # xVDot = self.CFDSolver.getSpatialPerturbation(321)
        wDot = self.CFDSolver.getStatePerturbation(314)


        resDot, funcsDot, fDot, hfDot = self.CFDSolver.computeJacobianVectorProductFwd(
        wDot=wDot, residualDeriv=True, funcDeriv=True, fDeriv=True, hfDeriv=True, mode='AD')

        resDot_FD, funcsDot_FD, fDot_FD, hfDot_FD = self.CFDSolver.computeJacobianVectorProductFwd(
        wDot=wDot, residualDeriv=True, funcDeriv=True, fDeriv=True, hfDeriv=True, mode='FD', h=1e-10)



        # print(resDot_FD[0:10])
        print(resDot[1944])
        # print(funcsDot)
        # print(resDot[0:10])
        import pickle
        with open('resDot_CS_heatxfer.p', 'rb') as f:
            resDot_CS = pickle.load(f)

        # these are checked against CS becuase the derivates appear to be poorly conditioned
        # import pickle
        # with open('resDot_CS.p', 'rb') as f:
        #     resDot_CS = pickle.load(f)

        res = (resDot_CS - resDot )/resDot
        idx =  np.argmax(res)
        print(np.max(res),idx, resDot[idx], resDot_CS[idx])

        np.testing.assert_allclose(resDot_CS, resDot, rtol=5e-12)


        for func in funcsDot:
            np.testing.assert_allclose(funcsDot_FD[func],funcsDot[func], err_msg=func,  rtol=1e-5)

        np.testing.assert_allclose(fDot_FD,fDot, rtol=5e-5)

        np.testing.assert_allclose(hfDot_FD[hfDot_FD != 0],hfDot[hfDot_FD != 0], rtol=5e-4)

    @unittest.skip("")
    def test_fwd_W_CS(self):
        "test the FWD derivatives of the functional"
        from adflow.pyADflow_C import ADFLOW_C
        self.CFDSolver = ADFLOW_C(options=self.aero_options)
        self.CFDSolver.getResidual(self.ap)
        # self.CFDSolver(self.ap)

        wDot = self.CFDSolver.getStatePerturbation(314)

        # resDot, funcsDot, fDot, hfDot = self.CFDSolver.computeJacobianVectorProductFwd(
        # xDvDot=xDvDot, residualDeriv=True, funcDeriv=True, fDeriv=True, hfDeriv=True, mode='AD')

        resDot_CS, funcsDot_CS, fDot_CS, hfDot_CS = self.CFDSolver.computeJacobianVectorProductFwd(
        wDot=wDot, residualDeriv=True, funcDeriv=True, fDeriv=True, hfDeriv=True, mode='CS', h=1e-200)
        print(resDot_CS[1944])

        import pickle
        with open('resDot_CS_heatxfer.p', 'wb') as f:
            pickle.dump(resDot_CS, f)

    # @unittest.skip("")
    def test_bwd_W(self):
        "test the BWD derivatives wrt the states"



        wDot = self.CFDSolver.getStatePerturbation(314)

        resDot, funcsDot, fDot, hfDot = self.CFDSolver.computeJacobianVectorProductFwd(
        wDot=wDot, residualDeriv=True, funcDeriv=True, fDeriv=True, hfDeriv=True, mode='AD')



        dwBar = self.CFDSolver.getStatePerturbation(314)
        fBar = self.CFDSolver.getSurfacePerturbation(314)
        hxferBar = fBar[:,0]


        funcsBar = {}
        for key in self.ap.evalFuncs:
            funcsBar[key] =  1.0


        wBar = self.CFDSolver.computeJacobianVectorProductBwd(
            resBar=dwBar, fBar=fBar, hfBar=hxferBar, funcsBar=funcsBar ,  wDeriv=True)

        # we have to add up all the parts
        _sum = 0
        _sum += np.dot(resDot, dwBar)
        _sum += np.dot(fDot.flatten(), fBar.flatten())
        _sum += np.dot(hfDot.flatten(), hxferBar.flatten())
        _sum += np.sum([x for x in funcsDot.values()])


        np.testing.assert_allclose(np.dot(wDot, wBar), _sum, rtol=5e-14)







#endregion


if __name__ == '__main__':
    unittest.main()
