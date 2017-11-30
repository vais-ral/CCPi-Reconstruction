'''
This is CT instrument information. it includes cone beam and parallel beam instrument information
'''
import numpy as np
import sys
import os.path
from PIL import Image
from ccpi.reconstruction.FindCenterOfRotation import find_center_vo
from ccpi.reconstruction.parallelbeam import alg

class Instrument(object):
    """This is a base class for the instrument type. This will act as an interface between the C/C++ source
    for the reconstruction and the python wrapper. 
    
    Attributes:
        pixels (numpy.array): This is to store the pixel information of the project data in numpy array
        angles (numpy.array): This is to store the projection angles in radians        
        
    """
    def __init__(self, pixels, angles):
        """Initialises the pixels and angles for the instrument.            
        """
        self.description = 'Base Instrument Class'
        self.pixels = pixels
        self.angles = angles
        self.data_v_offset = 0
        self.acceptedInputKeywords = ['pixels', 'angles', 'volume']
        self.pars = {'pixels':pixels, 'angles':angles ,'data_v_offset':self.data_v_offset}

    def read(filename):
        """This is a base method to read the instrument setup, pixels data and angles data
        required for the reconstruction and populate the object.
            
        Args:
            filename (str): Name of the configuration file/data file
        """
        raise NotImplementedError('Read method is not implemented for this instrument')
                
    def is_cone_beam(self):
        """This is to check if the instrument type is cone beam or parallel beam
        This will be overwritten by the child class
        """
        raise NotImplementedError('is_cone_beam needs to be overwritten for particular instrument')        
        
    def calc_v_alighment(self, num_v_pixels, pixels_per_voxel):
        """This will calculate the number of pixels based on the number of vertical pixels and pixels
        per voxel
        
        Args:
            num_v_pixels (int): The number of pixels in vertical dimension
            pixels_per_voxels (int): The number of pixels in a voxel.
        Returns:
            int value corresponding to new number of pixels in vertical dimension.
        """
        nvox = int(num_v_pixels/pixels_per_voxel)
        if nvox%pixels_per_voxel != 0:
            nvox+=1
        if self.is_cone_beam():
            if nvox%2 !=0:
                nvox+=1
        npix = nvox * pixels_per_voxel
        data_v_offset = (npix-num_v_pixels)/2
        data_v_size = num_v_pixels
        return npix          
        
    def __str__(self):
        return "Instrument: %s" % (self.description)
        
    def doForwardProject(self):
        raise NotImplementedError('doForwardProject must be implemented in the concrete class')
    
    def doBackwardProject(self):
        raise NotImplementedError('doBackwardProject must be implemented in the concrete class')
        
    def setParameter(self, **kwargs):
        '''set named parameter for the instrument
        
        raises Exception if the named parameter is not recognized
        
        '''
        for key , value in kwargs.items():
            if key in self.acceptedInputKeywords:
                self.pars[key] = value
            else:
                raise Exception('Wrong parameter {0} for '.format(key) +
                                'reconstructor')
    # setParameter
    def getParameter(self, key):
        '''Get a named parameter or list of parameter'''
        if type(key) is str:
            if key in self.acceptedInputKeywords:
                return self.pars[key]
            else:
                raise Exception('Unrecongnised parameter: {0} '.format(key) )
        elif type(key) is list:
            outpars = []
            for k in key:
                outpars.append(self.getParameter(k))
            return outpars
        else:
            raise Exception('Unhandled input {0}' .format(str(type(key))))
        
class Diamond(Instrument):
    """This represents diamond instrument (parallel beam). 
    
    TODO: Implement the reading and populating the object to be passed as argument to the reconstruction.
    """
    def __init__(self, pixels=None, angles=None):
        Instrument.__init__(self, pixels, angles)
        self.acceptedInputKeywords.append('pixels_per_voxel')
        self.acceptedInputKeywords.append('center_of_rotation')
        self.pars['pixels_per_voxel'] = 1
    
    def is_cone_beam(self):
        return False
        
    def read(filename):
        pass
    
    def doForwardProject(self, volume, angles, pixel_per_voxel=1, 
                         normalized=False, negative=False):
        '''Performs a forward projection'''
        self.setParameter(angles=angles)
        pixels = alg.pb_forward_project(volume, 
                                              angles , 
                                              pixel_per_voxel)
        self.setParameter(pixels=pixels)
        if normalized:
            m = pixels.min()
            M = pixels.max()
            scale = 1 / (M-m)
            shift = -m 
            print ("m,M,scale,shift" , m,M,scale,shift)
            pixels = pixels * scale + shift
            if negative:
                pixels = 1-pixels
        return pixels
    
    def doBackwardProject(self, center_of_rotation=None, pixel_per_voxel=1):
        '''Does a backward projection'''
        self.pixels , self.angles = self.getParameter(['pixels','angles'])
        if center_of_rotation is None:
            center_of_rotation = find_center_vo(self.pixels)
            print (self.acceptedInputKeywords)
            self.setParameter(center_of_rotation=center_of_rotation)
        
        back = alg.pb_backward_project(self.pixels, 
                                   self.angles, 
                                   center_of_rotation, 
                                   pixel_per_voxel)
        return back
    def find_center_of_rotation(self, pixels=None):
        if pixels is None:
            pixels = self.getParameter('pixels')
            return self.find_center_of_rotation(pixels)
        else:
            return find_center_vo(pixels)
        
        
class Xtek(Instrument):
    """This represents Xtek instrument. This class implements the reading of data from the Xtek output directory
    This class inherits Instrument class. 
    eg: >> x = Xtek()
        >> x.read('data/xtek/test.xtekct')
    """
    def __init__(self, pixels=None, angles=None, pixels_per_voxel=1):
        """
        """
        super(Xtek, self).__init__(pixels, angles)
        self.experiment_name = ''
        self.pixels_per_voxel = pixels_per_voxel #passed as argument
        self.source_x = 0.0
        self.detector_x = 0.0
        self.pixel_h_size = 1.0
        self.pixel_v_size = 1.0
        self.m_radius = 1.0
        self.num_of_vertical_pixels = 0
        self.num_of_horizontal_pixels = 0
        self.has_offsets = False        
        
    def is_cone_beam(self):
        """
        Returns: True since Xtek is cone beam based instrument.
        """
        return True
        
        
    def read_angles(self, filename, initial_angle, number_of_projections):
        """
        """
        input_path = os.path.dirname(filename)
        angles_ctdta_file = os.path.join(input_path, '_ctdata.txt')
        angles_named_file = os.path.join(input_path, self.experiment_name+'.ang')
        self.angles = np.zeros(number_of_projections,dtype='f')
        #look for _ctdata.txt
        if os.path.exists(angles_ctdta_file):
            #read txt file with angles
            with open(angles_file) as f:
                content = f.readlines()
            #skip firt three lines
            #read the middle value of 3 values in each line as angles in degrees
            index = 0
            for line in content[3:]:
                self.angles[index]=float(line.split(' ')[1])
                index+=1
            self.angles = np.deg2rad(self.angles+initial_angle);
        elif os.path.exists(angles_named_file):
            #read the angles file which is text with first line as header
            with open(angles_named_file) as f:
                content = f.readlines()
            #skip first line
            index = 0
            for line in content[1:]:
                self.angles[index] = float(line.split(':')[1])
                index+=1
            self.angles = np.flipud(self.angles+initial_angle) #angles are in the reverse order
        else:
            raise RuntimeError("Can't find angles file")
        pass
     
    def read_images(self, filename, number_of_projections, white_level, scattering):
        """
        """
        input_path = os.path.dirname(filename)
        self.pixels = np.zeros((number_of_projections, self.num_of_horizontal_pixels, self.num_of_vertical_pixels), dtype='float32')
        for i in range(1, number_of_projections+1):
            im = Image.open(os.path.join(input_path,self.experiment_name+"_%04d"%i+".tif"))
            self.pixels[i-1,:,:] = np.fliplr(np.transpose(np.array(im))) ##Not sure this is the correct way to populate the image
            
        max_v = np.amax(self.pixels)
        self.pixels = self.pixels - (white_level*scattering)/100.0
        self.pixels[self.pixels < 1.0] = 0.000001
        self.pixels = - np.log(self.pixels/max_v)
        
    def read(self, filename):
        """
        """
        content = []
        xpixel_size = 0
        ypixel_size = 0
        num_projections = 0
        initial_angle = 0
        scattering = 0
        white_level = 0
        with open(filename) as f:
            content = f.readlines()
        content = [x.strip() for x in content]
        for line in content:
            if line.startswith("SrcToObject"):
                self.source_x = float(line.split('=')[1])
            elif line.startswith("SrcToDetector"):
                self.detector_x = float(line.split('=')[1])
            elif line.startswith("DetectorPixelsY"):
                self.num_of_vertical_pixels = int(line.split('=')[1])
                self.num_of_vertical_pixels = self.calc_v_alighment(self.num_of_vertical_pixels, self.pixels_per_voxel)
            elif line.startswith("DetectorPixelsX"):
                self.num_of_horizontal_pixels = int(line.split('=')[1])
            elif line.startswith("DetectorPixelSizeX"):
                xpixel_size = float(line.split('=')[1])
            elif line.startswith("DetectorPixelSizeY"):
                ypixel_size = float(line.split('=')[1])   
            elif line.startswith("Projections"):
                num_projections = int(line.split('=')[1])
            elif line.startswith("InitialAngle"):
                initial_angle = float(line.split('=')[1])
            elif line.startswith("Name"):
                self.experiment_name = line.split('=')[1]
            elif line.startswith("Scattering"):
                scattering = float(line.split('=')[1])
            elif line.startswith("WhiteLevel"):
                white_level = float(line.split('=')[1])                
            elif line.startswith("MaskRadius"):
                self.mask_radius = float(line.split('=')[1])
            
            
                
        self.h_pixels = np.arange(self.num_of_horizontal_pixels, dtype='f')
        self.v_pixels = np.arange(self.num_of_vertical_pixels, dtype='f')
        pixel_base = -((self.num_of_horizontal_pixels-1)*xpixel_size/2.0)        
        self.h_pixels = self.h_pixels * xpixel_size + pixel_base 
        pixel_base = -((self.num_of_vertical_pixels-1)*ypixel_size/2.0)
        self.v_pixels = self.v_pixels * ypixel_size + pixel_base         
        
        self.xpixel_size = xpixel_size
        self.ypixel_size = ypixel_size
        self.read_angles(filename, initial_angle, num_projections)
        self.read_images(filename, num_projections, white_level, scattering)