import unittest
import numpy as np
import tomophantom
import tomophantom.phantom3d
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
        
    @unittest.skip("takes a long to process")
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
        
    def prepare_phantom(self):
        num_angles   = 250
        num_h_pixels = 512
        num_v_pixels = 512        
        
        source_x = -250.0
        source_y = 0
        source_z = 0
        detector_x = 987.0-250.0
        
        xpixel_size  = 0.390625
        ypixel_size  = 0.390625
        
        h_pixels = np.arange(num_h_pixels, dtype='float32')
        v_pixels = np.arange(num_v_pixels, dtype='float32')
        pixel_base = -((num_h_pixels-1)*xpixel_size/2.0)        
        h_pixels = h_pixels * xpixel_size + pixel_base 
        pixel_base = -((num_v_pixels-1)*ypixel_size/2.0)
        v_pixels = v_pixels * ypixel_size + pixel_base           
        
        mask_radius = -source_x * np.sin(np.arctan(h_pixels[num_h_pixels - 1] / (detector_x- source_x)))
        
        grid_offset=np.zeros(3, dtype='float32')
        voxel_size=np.ones(3, dtype='float32')
        image_vol = np.zeros(3, dtype='float32')        
        voxel_size[0] = (2 * mask_radius / num_h_pixels);
        voxel_size[1] = xpixel_size;
        voxel_size[2] = ypixel_size;

        image_vol[0] = voxel_size[0] * num_h_pixels;
        image_vol[1] = voxel_size[1] * num_h_pixels;
        image_vol[2] = voxel_size[2] * num_h_pixels;
        grid_offset[0] = -image_vol[0] / 2;
        grid_offset[1] = -image_vol[1] / 2;
        grid_offset[2] = -image_vol[2] / 2;
        
        #create a phantom dataset
        angles = np.deg2rad(np.linspace(0,360, num_angles, dtype='float32'))
        params = np.array([(1, 1.00, 0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 0.0),], dtype=[('Obj', np.int), ('C0', np.float32), ('x0', np.float32), ('y0',np.float32), ('z0', np.float32),('a',np.float32), ('b', np.float32), ('c', np.float32), ('phi', np.float32)])
        voxel_data = tomophantom.phantom3d.build_volume_phantom_3d_params(num_h_pixels, params)
        
        return voxel_data, angles, h_pixels, v_pixels, source_x, source_y, source_z, dectector_x, xpixel_size, ypixel_size, grid_offset, voxel_size
        
    def test_create_sinogram(self):
#prepare the phantom dataset
        params = self.prepare_phantom()
        source_x = params[4]
        source_y = params[5]
        source_z = params[6]
        dectector_x = params[7]
        h_pixels = params[2]
        v_pixels = params[3]
        angles   = params[1]
        voxel_data = params[0]
        grid_offset = params[8]
        voxel_size = params[9]
#forward project to create a sinogram        
        pixels = ccpi_reconstruction.create_sinogram(source_x, source_y, source_z, detector_x, h_pixels, v_pixels, angles, voxel_data, grid_offset, voxel_size)        
#reconstruct to generate 3d volume 
        
#compare the results
        
if __name__ == '__main__':
    unittest.main()
