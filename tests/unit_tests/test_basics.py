from __future__ import print_function
import unittest
import numpy
baseDir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(baseDir,'../../'))
from python.pyADflow import ADFLOW

class BasicTests(unittest.TestCase):

    def setup(self):
        gridFile = 'input_files/mdo_tutorial_euler.cgns'
        options = {'gridfile': gridFile}
        self.CFDSolver = ADFLOW(options=options)

    def test_import(self):
        "Tests if mesh was read properly"
        nstate = self.CFDSolver.getStateSize()
        assert(nstate == 60480)

if __name__ == '__main__':
    unittest.main()
