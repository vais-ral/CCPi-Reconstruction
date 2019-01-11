# -*- coding: utf-8 -*-
"""
Created on Fri Jun  2 15:48:00 2017

@author: ofn77899
"""

import numpy
#import h5py
from ccpi.reconstruction.parallelbeam import alg

from enum import Enum
from ccpi.viewer.CILViewer2D import Converter, CILViewer2D
#from FindCenterOfRotation import *
import vtk

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
        if (iterations != None):
            self.pars['iterationValues'] = numpy.zeros((iterations)) 
        
        if projection_data != None and dark_field != None and flat_field != None:
            norm = self.normalize(projection_data, dark_field, flat_field, 0.1)
            self.pars['normalized_projection'] = norm
            
    
    def setPars(self, parameters):
        keys = ['algorithm','projection_data' ,'normalized_projection', \
                'angles' , 'center_of_rotation' , 'flat_field', \
                'iterations','dark_field' , 'resolution', 'isLogScale' , \
                'threads' , 'iterationValues', 'regularize']
        
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
               self.pars = parameters
               result = parameters['algorithm'](
                           parameters['normalized_projection'] ,
                           parameters['angles'],
                           parameters['center_of_rotation'],
                           parameters['resolution'],
                           parameters['iterations'],
                           parameters['threads'] ,
                           parameters['regularize'],
                           numpy.zeros((parameters['iterations'])),
                           parameters['isLogScale']
                           )
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
        


class ProjectionPreprocessor():
    def __init__(self):
        self.flatsFileNames = []
        self.darksFileNames = []
        self.projectionsFileNames = []
        self.normalizedProjections = None
    
    def setFlatsFileNames(self, filenames):
        self.flatsFileNames = [i for i in filenames]
        self.normalizedProjections = None
    # setFlatsFileNames

    def setDarksFileNames(self, filenames):
        self.darksFileNames = [i for i in filenames]
        self.normalizedProjections = None
    # setFlatsFileNames

    def setProjectionsFileNames(self, filenames):
        self.projectionsFileNames = [i for i in filenames]
        self.normalizedProjections = None
    # setFlatsFileNames

    def avgProjections(self, filenames):
        '''Returns the average of the input images
        
        Input is passed as filenames pointing to 2D TIFF
        Output is numpy array
        '''
        num = len(filenames)
        stack = Converter.tiffStack2numpy(filenames=filenames , extent=extent,
                                          sampleRate=(1,1,1))
        avg = (stack.T[0] / float(num) )
        
        for i in range(1,num):
        	avg += (stack.T[i] / float(num))
        
        return avg
    
    def avgProjectionsVtk(self, filenames , extent):
        '''Returns the average of the input images
        
        Input is passed as filenames pointing to 2D TIFF
        Output is a vtkImageData
        gives the same result as avgProjections within machine precision.
        '''
        num = len(filenames)
        readers = [vtk.vtkTIFFReader() for i in range(num)]
        scalers = [vtk.vtkImageShiftScale() for i in range(num)]
        vois = [vtk.vtkExtractVOI() for i in range(num)]
        avger = vtk.vtkImageWeightedSum()
        counter = 0
        for reader, voi, scaler, filename in zip(readers, vois , scalers ,filenames):
            reader.SetFileName(filename)
            reader.Update()
            print ("reading {0}".format(filename))
            voi.SetInputData(reader.GetOutput())
            voi.SetVOI(extent)
            voi.Update()
            print ("extracting VOI")
            scaler.SetInputData(voi.GetOutput())
            scaler.SetOutputScalarTypeToDouble()
            scaler.SetScale(1)
            scaler.SetShift(0)
            scaler.Update()
            avger.AddInputData(scaler.GetOutput())
            avger.SetWeight(counter , 1/float(num))
            counter = counter + 1
            print ("adding input connection {0}".format(counter))
            
        
        avger.Update()
        return avger.GetOutput()
    
    def normalizeProjections(self, projectionFileNames=None, 
                             flatFileNames=None,
                             darkFileNames=None,
                             extent=None, sampleRate=None):
        if projectionFileNames is not None:
            self.setProjectionsFileNames(projectionFileNames)
        if flatFileNames is not None:
            self.setFlatsFileNames(flatFileNames)
        if darkFileNames is not None:
            self.setDarksFileNames(darkFileNames)
            
        flat = self.avgProjections(self.flatsFileNames)
        dark = self.avgProjections(self.darksFileNames)
        
        self.normalizedProjections = Converter.tiffStack2numpy(filenames=self.projectionsFileNames, 
                                         darkField=dark, flatField=flat, 
                                         extent=extent, sampleRate=sampleRate)
        
    
    def getNormalizedProjections(self):
        gg = numpy.frompyfunc(ProjectionPreprocessor.greater,2,1)
        ss = numpy.frompyfunc(ProjectionPreprocessor.smaller,2,1)
        eq = numpy.frompyfunc(ProjectionPreprocessor.equal,2,1)
        
        # set value > 1 to 1
        gt = ss(self.normalizedProjections, 1)
        self.normalizedProjections = numpy.where(gt, self.normalizedProjections, 1)
        
        # set values < 0 to 0
        lt = gg(self.normalizedProjections, 0)
        self.normalizedProjections = numpy.where(lt, self.normalizedProjections, 0.1)
        
        # remove 0 values
        et = eq(self.normalizedProjections,0)
        self.normalizedProjections = numpy.where(et , self.normalizedProjections, 0.1)
        
        
        return self.normalizedProjections.astype(numpy.float32)
    
    @staticmethod    
    def greater(x,y):
        return x>y
    @staticmethod
    def smaller(x,y):
        return x<y
    @staticmethod
    def equal(x,y):
        return x == y


def getEntry(location):
    for item in nx[location].keys():
        print (item)


print ("Loading Data")
 


#recon = Reconstructor(algorithm = Algorithm.CGLS, normalized_projection = norm,
#                 angles = angle_proj, center_of_rotation = 86.2 , 
#                 flat_field = flat, dark_field = dark, 
#                 iterations = 15, resolution = 1, isLogScale = False, threads = 3)

#recon = Reconstructor(algorithm = Reconstructor.Algorithm.CGLS, projection_data = proj,
#                 angles = angle_proj, center_of_rotation = 86.2 , 
#                 flat_field = flat, dark_field = dark, 
#                 iterations = 15, resolution = 1, isLogScale = False, threads = 3)
#img_cgls = recon.reconstruct()
#
#pars = dict()
#pars['algorithm'] = Reconstructor.Algorithm.SIRT
#pars['projection_data'] = proj
#pars['angles'] = angle_proj
#pars['center_of_rotation'] = numpy.double(86.2)
#pars['flat_field'] = flat
#pars['iterations'] = 15
#pars['dark_field'] = dark
#pars['resolution'] = 1
#pars['isLogScale'] = False
#pars['threads'] = 3
#
#img_sirt = recon.reconstruct(pars)
#
#recon.pars['algorithm'] = Reconstructor.Algorithm.MLEM
#img_mlem = recon.reconstruct()

############################################################


tot_images = 2203

## create a reduced 3D stack
nreduced = 360
directory = r'D:\Documents\Dataset\Chadwick_Flange_Tomo'

indices = [int(float(num) / float(nreduced) * float(tot_images)) for num in range(nreduced)]
angle_proj = numpy.asarray(indices) * 180 / tot_images

#load data file
reader = vtk.vtkMetaImageReader()
reader.SetFileName("Chadwick_flange_norm_360.mha")
reader.Update()

norm = Converter.vtk2numpy(reader.GetOutput())

pp = ProjectionPreprocessor()
pp.normalizedProjections = norm
norm = pp.getNormalizedProjections()

cor = 175.75
############################################################
#recon.pars['algorithm'] = Reconstructor.Algorithm.CGLS_CONV
#recon.pars['regularize'] = numpy.double(0.1)
#img_cgls_conv = recon.reconstruct()

niterations = 15
threads = 4

data = (norm.T).copy()
del norm
del reader
norm = data

print ("Launching CGLS")
img_cgls = alg.cgls(norm, angle_proj, numpy.double(cor), 1 , niterations, threads, False)
numpy.save("chadwick_chamber.npy", img_cgls)
print ("CGLS finished")

v = CILViewer2D()
v.setInputAsNumpy(img_cgls)
v.startRenderLoop()
##img_mlem = alg.mlem(norm, angle_proj, numpy.double(cor), 1 , niterations, threads, False)
##img_sirt = alg.sirt(norm, angle_proj, numpy.double(cor), 1 , niterations, threads, False)
##
##iteration_values = numpy.zeros((niterations,))
##img_cgls_conv = alg.cgls_conv(norm, angle_proj, numpy.double(cor), 1 , niterations, threads,
##                              iteration_values, False)
##print ("iteration values %s" % str(iteration_values))
##
##iteration_values = numpy.zeros((niterations,))
##img_cgls_tikhonov = alg.cgls_tikhonov(norm, angle_proj, numpy.double(cor), 1 , niterations, threads,
##                                      numpy.double(1e-5), iteration_values , False)
##print ("iteration values %s" % str(iteration_values))
##iteration_values = numpy.zeros((niterations,))
##img_cgls_TVreg = alg.cgls_TVreg(norm, angle_proj, numpy.double(cor), 1 , niterations, threads,
##                                      numpy.double(1e-5), iteration_values , False)
##print ("iteration values %s" % str(iteration_values))
##
##
####numpy.save("cgls_recon.npy", img_data)
##import matplotlib.pyplot as plt
##fig, ax = plt.subplots(1,6,sharey=True)
##ax[0].imshow(img_cgls[80])
##ax[0].axis('off')  # clear x- and y-axes
##ax[1].imshow(img_sirt[80])
##ax[1].axis('off')  # clear x- and y-axes
##ax[2].imshow(img_mlem[80])
##ax[2].axis('off')  # clear x- and y-axesplt.show()
##ax[3].imshow(img_cgls_conv[80])
##ax[3].axis('off')  # clear x- and y-axesplt.show()
##ax[4].imshow(img_cgls_tikhonov[80])
##ax[4].axis('off')  # clear x- and y-axesplt.show()
##ax[5].imshow(img_cgls_TVreg[80])
##ax[5].axis('off')  # clear x- and y-axesplt.show()
##
##
###plt.show()
##fig.savefig('recon3d.png')
##plt.close(fig)
##
