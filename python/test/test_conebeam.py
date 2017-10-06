import unittest
import numpy as np
from ccpi.reconstruction.conebeam import alg as ccpi_reconstruction

class TestConebeamReconstruction(unittest.TestCase):
    """
    Test the cone beam reconstruction
    """
    def setUp(self):
        pass
    
    def tearDown(self):
        pass
        
    def test_create_phantom(self):
        data = ccpi_reconstruction.create_phantom()
        self.assertEqual(data[0].shape[0], 250)
        self.assertEqual(data[0].shape[1], 512)
        self.assertEqual(data[0].shape[2], 512)
        self.assertEqual(data[1].shape[0], 250)
        self.assertEqual(data[2], 250.0)
        self.assertEqual(data[3], 987.0)
        self.assertEqual(data[4], 0.390625)
        self.assertEqual(data[5], 0.390625)
        self.assertEqual(data[6], 25.151547886198443)

    def test_cgls_test(self):
        #load phantom data
        #projection = np.load('cone_phantom_data.npy')
        #angles = np.load('cone_phantom_angles.npy')
        
        source_x = 250.0
        detector_x = 987.0
        h_pixel_size = 0.390625
        v_pixel_size = 0.390625
        mask_radius = 16.0156431850479
        