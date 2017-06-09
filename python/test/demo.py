# -*- coding: utf-8 -*-
"""
Created on Fri Jun  2 15:48:00 2017

@author: ofn77899
"""

import numpy
import h5py
from ccpi.reconstruction.parallelbeam import alg
import edo

from enum import Enum



class Reconstructor:
    
    class Algorithm(Enum):
        CGLS = alg.cgls
        CGLS_CONV = alg.cgls_conv
        SIRT = alg.sirt
        MLEM = alg.mlem
        CGLS_TICHONOV = alg.cgls_tikhonov
        CGLS_TVREG = alg.cgls_TVreg
        
    def __init__(self, algorithm = None, projection_data = None,
                 angles = None, center_of_rotation = None , 
                 flat_field = None, dark_field = None, 
                 iterations = None, resolution = None, isLogScale = False, threads = None, 
                 normalized_projection = None):
    
        self.pars = dict()
        self.pars['algorithm'] = algorithm
        self.pars['projection_data'] = projection_data
        self.pars['normalized_projection'] = normalized_projection
        self.pars['angles'] = angles
        self.pars['center_of_rotation'] = numpy.double(center_of_rotation)
        self.pars['flat_field'] = flat_field
        self.pars['iterations'] = iterations
        self.pars['dark_field'] = dark_field
        self.pars['resolution'] = resolution
        self.pars['isLogScale'] = isLogScale
        self.pars['threads'] = threads
        
        if projection_data != None and dark_field != None and flat_field != None:
            norm = self.normalize(projection_data, dark_field, flat_field, 0.1)
            self.pars['normalized_projection'] = norm
            
    
    def setPars(self, parameters):
        keys = ['algorithm','projection_data' ,'normalized_projection', \
                'angles' , 'center_of_rotation' , 'flat_field', \
                'iterations','dark_field' , 'resolution', 'isLogScale' , \
                'threads']
        
        for k in keys:
            if k not in parameters.keys():
                self.pars[k] = None
            else:
                self.pars[k] = parameters[k]
                
        
    def sanityCheck(self):
        projection_data = self.pars['projection_data']
        dark_field = self.pars['dark_field']
        flat_field = self.pars['flat_field']
        angles = self.pars['angles']
        
        if projection_data != None and dark_field != None and \
            angles != None and flat_field != None:
            data_shape =  numpy.shape(projection_data)
            angle_shape = numpy.shape(angles)
            
            if angle_shape[0] != data_shape[0]:
                #raise Exception('Projections and angles dimensions do not match: %d vs %d' % \
                #                (angle_shape[0] , data_shape[0]) )
                return (False , 'Projections and angles dimensions do not match: %d vs %d' % \
                                (angle_shape[0] , data_shape[0]) )
            
            if data_shape[1:] != numpy.shape(flat_field):
                #raise Exception('Projection and flat field dimensions do not match')
                return (False , 'Projection and flat field dimensions do not match')
            if data_shape[1:] != numpy.shape(dark_field):
                #raise Exception('Projection and dark field dimensions do not match')
                return (False , 'Projection and dark field dimensions do not match')
            
            return (True , '' )
        elif self.pars['normalized_projection'] != None:
            data_shape =  numpy.shape(self.pars['normalized_projection'])
            angle_shape = numpy.shape(angles)
            
            if angle_shape[0] != data_shape[0]:
                #raise Exception('Projections and angles dimensions do not match: %d vs %d' % \
                #                (angle_shape[0] , data_shape[0]) )
                return (False , 'Projections and angles dimensions do not match: %d vs %d' % \
                                (angle_shape[0] , data_shape[0]) )
            else:
                return (True , '' )
        else:
            return (False , 'Not enough data')
            
    def reconstruct(self, parameters = None):
        if parameters != None:
            self.setPars(parameters)
        
        go , reason = self.sanityCheck()
        if go:
            return self._reconstruct()
        else:
            raise Exception(reason)
            
            
    def _reconstruct(self, parameters=None):
        if parameters!=None:
            self.setPars(parameters)
        parameters = self.pars
        
        if parameters['algorithm'] != None and \
           parameters['normalized_projection'] != None and \
           parameters['angles'] != None and \
           parameters['center_of_rotation'] != None and \
           parameters['iterations'] != None and \
           parameters['resolution'] != None and\
           parameters['threads'] != None and\
           parameters['isLogScale'] != None:
               
               
           if parameters['algorithm'] in (Reconstructor.Algorithm.CGLS,
                        Reconstructor.Algorithm.MLEM, Reconstructor.Algorithm.SIRT):
               #store parameters
               self.pars = parameters
               result = parameters['algorithm'](
                           parameters['normalized_projection'] ,
                           parameters['angles'],
                           parameters['center_of_rotation'],
                           parameters['resolution'],
                           parameters['iterations'],
                           parameters['threads'] ,
                           parameters['isLogScale']
                           )
               return result
           else:
               pass
        else:
           if parameters['projection_data'] != None and \
                     parameters['dark_field'] != None and \
                     parameters['flat_field'] != None:
               norm = self.normalize(parameters['projection_data'],
                                   parameters['dark_field'], 
                                   parameters['flat_field'], 0.1)
               self.pars['normalized_projection'] = norm
               return self._reconstruct(parameters)
              
                
                
    def _normalize(self, projection, dark, flat, def_val=0):
        a = (projection - dark)
        b = (flat-dark)
        with numpy.errstate(divide='ignore', invalid='ignore'):
            c = numpy.true_divide( a, b )
            c[ ~ numpy.isfinite( c )] = def_val  # set to not zero if 0/0 
        return c
    
    def normalize(self, projections, dark, flat, def_val=0):
        norm = [self._normalize(projection, dark, flat, def_val) for projection in projections]
        return numpy.asarray (norm, dtype=numpy.float32)
        


def getEntry(location):
    for item in nx[location].keys():
        print (item)


print ("Loading Data")

##fname = "D:\\Documents\\Dataset\\IMAT\\20170419_crabtomo\\crabtomo\\Sample\\IMAT00005153_crabstomo_Sample_000.tif"
####ind = [i * 1049 for i in range(360)]
#### use only 360 images
##images = 200
##ind = [int(i * 1049 / images) for i in range(images)]
##stack_image = dxchange.reader.read_tiff_stack(fname, ind, digit=None, slc=None)

#fname = "D:\\Documents\\Dataset\\CGLS\\24737_fd.nxs"
fname = "C:\\Users\\ofn77899\\Documents\\CCPi\\CGLS\\24737_fd_2.nxs"
nx = h5py.File(fname, "r")

# the data are stored in a particular location in the hdf5
for item in nx['entry1/tomo_entry/data'].keys():
    print (item)

data = nx.get('entry1/tomo_entry/data/rotation_angle')
angles = numpy.zeros(data.shape)
data.read_direct(angles)
print (angles)
# angles should be in degrees

data = nx.get('entry1/tomo_entry/data/data')
stack = numpy.zeros(data.shape)
data.read_direct(stack)
print (data.shape)

print ("Data Loaded")


# Normalize
data = nx.get('entry1/tomo_entry/instrument/detector/image_key')
itype = numpy.zeros(data.shape)
data.read_direct(itype)
# 2 is dark field
darks = [stack[i] for i in range(len(itype)) if itype[i] == 2 ]
dark = darks[0]
for i in range(1, len(darks)):
    dark += darks[i]
dark = dark / len(darks)
#dark[0][0] = dark[0][1]

# 1 is flat field
flats = [stack[i] for i in range(len(itype)) if itype[i] == 1 ]
flat = flats[0]
for i in range(1, len(flats)):
    flat += flats[i]
flat = flat / len(flats)
#flat[0][0] = dark[0][1]


# 0 is projection data
proj = [stack[i] for i in range(len(itype)) if itype[i] == 0 ]
angle_proj = [angles[i] for i in range(len(itype)) if itype[i] == 0 ]
angle_proj = numpy.asarray (angle_proj)
angle_proj = angle_proj.astype(numpy.float32)

# normalized data are
# norm = (projection - dark)/(flat-dark)

def normalize(projection, dark, flat, def_val=0.1):
    a = (projection - dark)
    b = (flat-dark)
    with numpy.errstate(divide='ignore', invalid='ignore'):
        c = numpy.true_divide( a, b )
        c[ ~ numpy.isfinite( c )] = def_val  # set to not zero if 0/0 
    return c
    

norm = [normalize(projection, dark, flat) for projection in proj]
norm = numpy.asarray (norm)
norm = norm.astype(numpy.float32)

#recon = Reconstructor(algorithm = Algorithm.CGLS, normalized_projection = norm,
#                 angles = angle_proj, center_of_rotation = 86.2 , 
#                 flat_field = flat, dark_field = dark, 
#                 iterations = 15, resolution = 1, isLogScale = False, threads = 3)
recon = Reconstructor(algorithm = Reconstructor.Algorithm.CGLS, projection_data = proj,
                 angles = angle_proj, center_of_rotation = 86.2 , 
                 flat_field = flat, dark_field = dark, 
                 iterations = 15, resolution = 1, isLogScale = False, threads = 3)
img_cgls = recon.reconstruct()

pars = dict()
pars['algorithm'] = Reconstructor.Algorithm.SIRT
pars['projection_data'] = proj
pars['angles'] = angle_proj
pars['center_of_rotation'] = numpy.double(86.2)
pars['flat_field'] = flat
pars['iterations'] = 15
pars['dark_field'] = dark
pars['resolution'] = 1
pars['isLogScale'] = False
pars['threads'] = 3

img_sirt = recon.reconstruct(pars)

recon.pars['algorithm'] = Reconstructor.Algorithm.MLEM
img_mlem = recon.reconstruct()

##numpy.save("cgls_recon.npy", img_data)
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1,3,sharey=True)
ax[0].imshow(img_cgls[80])
ax[0].axis('off')  # clear x- and y-axes
ax[1].imshow(img_sirt[80])
ax[1].axis('off')  # clear x- and y-axes
ax[2].imshow(img_mlem[80])
ax[2].axis('off')  # clear x- and y-axesplt.show()
plt.show()

#viewer = edo.CILViewer()
#viewer.setInputAsNumpy(img_cgls2)
#viewer.displaySliceActor(0)
#viewer.startRenderLoop()

import vtk

def NumpyToVTKImageData(numpyarray):
    if (len(numpy.shape(numpyarray)) == 3):
        doubleImg = vtk.vtkImageData()
        shape = numpy.shape(numpyarray)
        doubleImg.SetDimensions(shape[0], shape[1], shape[2])
        doubleImg.SetOrigin(0,0,0)
        doubleImg.SetSpacing(1,1,1)
        doubleImg.SetExtent(0, shape[0]-1, 0, shape[1]-1, 0, shape[2]-1)
        #self.img3D.SetScalarType(vtk.VTK_UNSIGNED_SHORT, vtk.vtkInformation())
        doubleImg.AllocateScalars(vtk.VTK_DOUBLE,1)
        
        for i in range(shape[0]):
            for j in range(shape[1]):
                for k in range(shape[2]):
                    doubleImg.SetScalarComponentFromDouble(
                        i,j,k,0, numpyarray[i][j][k])
    #self.setInput3DData( numpy_support.numpy_to_vtk(numpyarray) )
        # rescale to appropriate VTK_UNSIGNED_SHORT
        stats = vtk.vtkImageAccumulate()
        stats.SetInputData(doubleImg)
        stats.Update()
        iMin = stats.GetMin()[0]
        iMax = stats.GetMax()[0]
        scale = vtk.VTK_UNSIGNED_SHORT_MAX / (iMax - iMin)

        shiftScaler = vtk.vtkImageShiftScale ()
        shiftScaler.SetInputData(doubleImg)
        shiftScaler.SetScale(scale)
        shiftScaler.SetShift(iMin)
        shiftScaler.SetOutputScalarType(vtk.VTK_UNSIGNED_SHORT)
        shiftScaler.Update()
        return shiftScaler.GetOutput()
        
#writer = vtk.vtkMetaImageWriter()
#writer.SetFileName(alg + "_recon.mha")
#writer.SetInputData(NumpyToVTKImageData(img_cgls2))
#writer.Write()
