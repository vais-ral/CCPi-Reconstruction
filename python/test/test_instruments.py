'''
Unit tests for instrument classes 
'''

import unittest
from ccpi.instrument import Xtek
    
class TestInstruments(unittest.TestCase):    
    def test_read_xtekdata(self):
        xtek = Xtek()
        xtek.read('data/xtek/SophiaBeads_64_averaged.xtekct')
        self.assertEqual(xtek.source_x, -80.6392412185669)
        self.assertEqual(xtek.detector_x, 926.3667587814331)
        self.assertEqual(xtek.num_of_vertical_pixels, 2000)
        self.assertEqual(xtek.num_of_horizontal_pixels, 2000) 
        self.assertEqual(xtek.h_pixels.size, xtek.num_of_horizontal_pixels)
        self.assertEqual(xtek.v_pixels.size, xtek.num_of_vertical_pixels)        
        self.assertNotEqual(xtek.h_pixels[0], 0)
        self.assertNotEqual(xtek.v_pixels[0], 0)
        self.assertEqual(xtek.angles.size, 63)
        self.assertEqual(xtek.pixels.shape, (2000,2000, 63))
        
        
if __name__ == "__main__":
    unittest.main()