from __future__ import print_function
import unittest
import numpy

class BasicTests(unittest.TestCase):

    def setup(self):
        from ... import ADFLOW
        gridFile = 'input_files/mdo_tutorial_euler.cgns'
        options = {'gridfile': gridFile}
        self.CFDSolver = ADFLOW(options=options)

    def test_import(self):
        "Tests if mesh was read properly"
        nstate = self.CFDSolver.getStateSize()
        assert(nstate == 60480)
