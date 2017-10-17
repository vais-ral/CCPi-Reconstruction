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
        data = ccpi_reconstruction.create_phantom()
        projection = data[0] #projections
        angles = data[1] #angles
        h_offsets = np.zeros(1, dtype='float32')
        v_offsets = np.zeros(1, dtype='float32')        
        source_x = data[2] # source_x
        detector_x = data[3] #detector_x
        h_pixel_size = data[4] #h pixel size
        v_pixel_size = data[5] #v pixel size
        mask_radius = data[6] #mask radius
        full_vox_origin = np.zeros(3, dtype='float32')
        voxel_size = np.ones(3,dtype='float32')
        niterations=10
        nthreads=16
        recon = ccpi_reconstruction.cgls(projection, angles, h_offsets, v_offsets, 1, source_x, detector_x,h_pixel_size,v_pixel_size,mask_radius,False,full_vox_origin,voxel_size,niterations,nthreads,True)
        self.assertNotEqual(recon[256,256,256],0)
        
if __name__ == '__main__':
    unittest.main()
