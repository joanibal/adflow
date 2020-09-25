# built-ins 
import unittest
import numpy
import os
import sys
import copy
from parameterized import parameterized_class

# MACH classes
from adflow import ADFLOW

# import the testing utilities that live a few directories up 
import reg_test_utils as utils
from reg_default_options import adflowDefOpts
from reg_aeroproblems import ap_tutorial_wing 
from reg_test_classes import test_objects


baseDir = os.path.dirname(os.path.abspath(__file__))


@parameterized_class([
    # scalar JST
    { "name": 'euler_scalar_JST_tut_wing',
        "options": {
            'gridfile':os.path.join(baseDir, '../../inputFiles/mdo_tutorial_euler_scalar_jst.cgns'),
            'l2convergence': 1e-14,
            'mgcycle':'2w',
            'ncyclescoarse':250,
            'usenksolver':True,
            'useblockettes': False,
            'nkswitchtol':1e-2,
            },
        "ref_file": 'solve_euler_scalar_jst_tut_wing.json', 
        'ap':ap_tutorial_wing,
    }, 
    # Matrix JST
    { "name": 'euler_matrix_JST_tut_wing',
        "options": {
        'gridfile':os.path.join(baseDir, '../../inputFiles/mdo_tutorial_euler_matrix.cgns'),

        'mgcycle':'2w',
        'ncyclescoarse':250,
        'usenksolver':True,
        'nkswitchtol':1e-2,

        'vis4':0.1,
        'discretization':'central plus matrix dissipation',
        'coarsediscretization':'central plus matrix dissipation',
        'l2convergence':1e-14,
        'useblockettes': False,

        },
        "ref_file":'euler_matrix_jst_tut_wing.json', 
        'ap':ap_tutorial_wing,
    },  
    # Tutorial wing RANS
    { 'name':'rans_tut_wing',
        "options": {
        'gridfile':os.path.join(baseDir, '../../inputFiles/mdo_tutorial_rans_scalar_jst.cgns'),
        
        'mgcycle':'sg',
        'equationtype':'RANS',
        'smoother':'dadi',
        'cfl':1.5,
        'cflcoarse':1.25,
        'resaveraging':'noresaveraging',
        'nsubiter':3,
        'nsubiterturb':3,
        'ncyclescoarse':100,
        'ncycles':1000,
        'monitorvariables':['cpu', 'resrho','resturb','cl','cd','cmz','yplus','totalr'],
        'usenksolver':True,
        'l2convergence':1e-14,
        'l2convergencecoarse':1e-4,
        'nkswitchtol':1e-3,
        'adjointl2convergence': 1e-14,
        'frozenturbulence':False,

        },
        "ref_file": 'solve_rans_tut_wing.json', 
        'ap':ap_tutorial_wing,
    }, 

])
class TestSolve(test_objects.RegTest):
    '''
    Tests that ADflow can converge the given test problems to the given tolerance.
    
    Compares the norm of the residuals, norm of states, and the functionals.

    based on the old regression tests 9, 11, 12
    '''
    N_PROCS = 4

    name = None
    
    def setUp(self, train=False):
        if self.name == None:
                # return immediately when the setup method is being called on the based class and NOT the 
                # classes created using parametrized 
                # this will happen when training, but will hopefully be fixed down the line
                return

        super().setUp()
        
        options = copy.copy(adflowDefOpts)
        options['outputdirectory'] = os.path.join(baseDir, options['outputdirectory'])
        options.update(self.options)

        # Create the solver
        self.CFDSolver = ADFLOW(options=options, debug=False)
    
    def test_solve(self):

        # do the solve
        self.CFDSolver(self.ap)

        # check its accuracy
        utils.assert_functions_allclose(self.handler, self.CFDSolver, self.ap)
        utils.assert_states_allclose(self.handler, self.CFDSolver)
        utils.assert_residuals_allclose(self.handler, self.CFDSolver, self.ap)


if __name__ == '__main__':
    unittest.main()