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


# def delete_module(modname, paranoid=None):
#     from sys import modules
#     try:
#         thismod = modules[modname]
#     except KeyError:
#         raise ValueError(modname)
#     these_symbols = dir(thismod)
#     if paranoid:
#         try:
#             paranoid[:]  # sequence support
#         except:
#             raise ValueError('must supply a finite list for paranoid')
#         else:
#             these_symbols = paranoid[:]
#     del modules[modname]
#     for mod in modules.values():
#         try:
#             delattr(mod, modname)
#         except AttributeError:
#             pass
#         if paranoid:
#             for symbol in these_symbols:
#                 if symbol[:2] == '__':  # ignore special symbols
#                     continue
#                 try:
#                     delattr(mod, symbol)
#                 except AttributeError:
#                     pass


class AreaClacTests(unittest.TestCase):

    def setUp(self):
        # the mesh is an extruded naca 0012 of unit chord and ref area
        gridFile = os.path.join(baseDir, '../input_files/cube_hot.cgns')

        self.aero_options = {
            # I/O Parameters
            'gridfile': gridFile,
            # "restartfile": gridFile,
            "equationType": "RANS",

            "writevolumesolution": False,
            'writesurfacesolution':False,


            "MGCycle": "sg",
            "MGStartLevel": -1,
            "useANKSolver": True,
            "useNKSolver": True,

            'L2Convergence':1e-14,
            'nCycles':4000,
            'adjointl2convergence': 1e-14,

            'useblockettes': False,

        }



        self.ap = AeroProblem(
            name="n0012",
            V=32,  # m/s
            T=273 + 60,  # kelvin
            P=93e3,  # pa
            areaRef=1.0,  # m^2
            chordRef=1.0,  # m^2
            evalFuncs=["cd", "totheattransfer", "havg", "area" ],
            alpha=0.0,
            beta=0.00,
            xRef=0.0,
            yRef=0.0,
            zRef=0.0,
        )

        # these values come from inspecting the cgns mesh itself
        self.CFDSolver = ADFLOW(options=self.aero_options)
        # self.CFDSolver.addFunction('area', 'heated_wall', name="hot_area")

        self.CFDSolver.getResidual(self.ap)
        # self.CFDSolver(self.ap)

    # @unittest.skip("")
    def test_area(self):
        "Tests the area is being calculated correctly"


        def area_of_n0012(a, b):


            def f(x):
                return 0.6*(0.2969*np.sqrt(x) - 0.1260*x - 0.3516*x**2 + 0.2843*x**3 - 0.1015*x**4)

            # break the segment up into 100 pts
            x = np.linspace(a,b)
            y = f(x)


            dx = x[1:] - x[:-1]
            dy = y[1:] - y[:-1]

            s = np.sqrt(dx**2 + dy**2)

            l = np.sum(s)*2 # x2 becuase there is a top and bottom

            return l


        funcs = {}
        self.CFDSolver.evalFunctions(self.ap, funcs)


        # low tolerance here because we are comparing an approx with another approx
        np.testing.assert_allclose(funcs['n0012_area'],area_of_n0012(0,1 ), rtol=5e-3)
        np.testing.assert_allclose(funcs['n0012_hot_area'],area_of_n0012(0.1, 0.6 ), rtol=1e-4)

    # @unittest.skip("")
    def test_fwd(self):

        # self.CFDSolver(self.ap)
        xVDot = self.CFDSolver.getSpatialPerturbation(123)

        resDot, funcsDot, fDot, hfDot = self.CFDSolver.computeJacobianVectorProductFwd(
        xVDot=xVDot, residualDeriv=True, funcDeriv=True, fDeriv=True, hfDeriv=True, mode='AD')

        resDot_FD, funcsDot_FD, fDot_FD, hfDot_FD = self.CFDSolver.computeJacobianVectorProductFwd(
            xVDot=xVDot, residualDeriv=True, funcDeriv=True, fDeriv=True, hfDeriv=True, mode='FD', h=1e-8
        )


        # res not accurate becuase of if-statement branching caused by perturbations
        # np.testing.assert_allclose(resDot_FD,resDot, atol=1e-8)

        for func in funcsDot:
            np.testing.assert_allclose(funcsDot_FD[func],funcsDot[func], err_msg=func,  rtol=1e-5)
            # print('func',func,  (funcsDot[func] - funcsDot_CS[func]), funcsDot[func])
            # print('func FD',func,  (funcsDot_FD[func] - funcsDot_CS[func]), funcsDot_FD[func])


        np.testing.assert_allclose(fDot_FD,fDot, atol=1e-6)
        np.testing.assert_allclose(hfDot_FD,hfDot, atol=1e-6)

    @unittest.skip("cmplx adflow needed")
    def test_fwd_CS(self):

        from adflow.pyADflow_C import ADFLOW_C


        # --- create a compolex solver ---
        self.CFDSolver = ADFLOW_C(options=self.aero_options)
        # self.CFDSolver.addFunction('area', 'heated_wall', name="hot_area")
        self.CFDSolver.getResidual(self.ap)

        # self.CFDSolver(ap)
        xVDot = self.CFDSolver.getSpatialPerturbation(123)

        resDot_CS, funcsDot_CS, fDot_CS, hfDot_CS = self.CFDSolver.computeJacobianVectorProductFwd(
            xVDot=xVDot, residualDeriv=True, funcDeriv=True, fDeriv=True, hfDeriv=True, mode='CS', h=1e-40
        )

        # res = (resDot - resDot_CS)/resDot
        import pickle
        with open('resDot_CS.p', 'wb') as f:
            pickle.dump(resDot_CS, f)
        with open('funcsDot_CS.p', 'wb') as f:
            pickle.dump(funcsDot_CS, f)
        with open('fDot_CS.p', 'wb') as f:
            pickle.dump(fDot_CS, f)
        with open('hfDot_CS.p', 'wb') as f:
            pickle.dump(hfDot_CS, f)

        # print(resDot_CS[109])
        # # max', np.max(res), np.argmax(res), 'res[109]', res[109] )
        # print(funcsDot_CS['area'])
        # import ipdb; ipdb.set_trace()
        # # quit()

        # # import ipdb; ipdb.set_trace()
        # np.testing.assert_allclose(resDot_CS,resDot, atol=1e-8)


        # for func in funcsDot:
        #     np.testing.assert_allclose(funcsDot_CS[func],funcsDot[func], err_msg=func,  rtol=1e-5)

        # # np.testing.assert_allclose(fDot_CS,fDot, atol=1e-10)
        # # np.testing.assert_allclose(hfDot_CS,hfDot, atol=1e-10)
        # CFDSolver(ap, writeSolution=False)

        # resDot, funcsDot, fDot, hfDot = CFDSolver.computeJacobianVectorProductFwd(
        # xDvDot=xDvDot, residualDeriv=True, funcDeriv=True, fDeriv=True, hfDeriv=True, mode='CS')
        # pp(funcsDot)
        # del(CFDSolver)
        # del(ADFLOW_C)

        # print(resDot)
        # print(funcsDot)
        # print(fDot)
        # print(hfDot)


    def test_bwd(self):

        xVDot = self.CFDSolver.getSpatialPerturbation(321)


        resDot, _, _, __ = self.CFDSolver.computeJacobianVectorProductFwd(
        xVDot=xVDot, residualDeriv=True, funcDeriv=True, fDeriv=True, hfDeriv=True, mode='AD')

        dwBar = self.CFDSolver.getStatePerturbation(314)
        xVBar = self.CFDSolver.computeJacobianVectorProductBwd(

if __name__ == '__main__':
    unittest.main()
