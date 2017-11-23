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

from FindCenterOfRotation import *

import os

class Reconstructor:
    
    class Algorithm(Enum):
        CGLS = alg.cgls
        CGLS_CONV = alg.cgls_conv
        SIRT = alg.sirt
        MLEM = alg.mlem
        CGLS_TIKHONOV = alg.cgls_tikhonov
        CGLS_TVREG = alg.cgls_TVreg
        
    def __init__(self, algorithm = None, projection_data = None,
                 angles = None, center_of_rotation = None , 
                 flat_field = None, dark_field = None, 
                 iterations = None, resolution = None, isLogScale = False, threads = None, 
                 normalized_projection_data = None):

        self.reconstructedVolume = numpy.asarray([])
        self.pars = dict()
        self.pars['algorithm'] = algorithm
        self.pars['projection_data'] = projection_data
        self.pars['normalized_projection_data'] = normalized_projection_data
        self.pars['angles'] = angles
        self.pars['center_of_rotation'] = numpy.double(center_of_rotation)
        self.pars['flat_field'] = flat_field
        self.pars['iterations'] = iterations
        self.pars['dark_field'] = dark_field
        self.pars['resolution'] = resolution
        self.pars['isLogScale'] = isLogScale
        self.pars['threads'] = threads
        self.pars['input_volume'] = None
        if (iterations != None):
            self.pars['iterationValues'] = numpy.zeros((iterations))

        self.preprocessor = ProjectionPreprocessor()
        
        if projection_data != None and dark_field != None and flat_field != None:
            norm = self.normalize(projection_data, dark_field, flat_field, 0.1)
            self.preprocessor.normalizedProjections = norm
            norm = self.preprocessor.clipNormalizedProjections()
            self.pars['normalized_projection_data'] = norm
            
    
    def setPars(self, parameters):
        keys = ['algorithm','projection_data' ,'normalized_projection_data', \
                'angles' , 'center_of_rotation' , 'flat_field', \
                'iterations','dark_field' , 'resolution', 'isLogScale' , \
                'threads' , 'iterationValues', 'regularize',
                'input_volume']
        
        for k in keys:
            if k not in parameters.keys():
                self.pars[k] = None
            else:
                self.pars[k] = parameters[k]
                
    def setParameter(self, **kwargs):
        '''set named parameter for the regularization engine
        
        raises Exception if the named parameter is not recognized
        Typical usage is:
            
        reg = Regularizer(Regularizer.Algorithm.SplitBregman_TV)
        reg.setParameter(input=u0)    
        reg.setParameter(regularization_parameter=10.)
        
        it can be also used as
        reg = Regularizer(Regularizer.Algorithm.SplitBregman_TV)
        reg.setParameter(input=u0 , regularization_parameter=10.)
        '''
        acceptedkeys = ['algorithm','projection_data' ,'normalized_projection_data', \
                'angles' , 'center_of_rotation' , 'flat_field', \
                'iterations','dark_field' , 'resolution', 'isLogScale' , \
                'threads' , 'iterationValues', 'regularize', 'input_volume']
        for key , value in kwargs.items():
            if key in acceptedkeys:
                self.pars[key] = value
            else:
                raise Exception('Wrong parameter {0} for '.format(key) +
                                'Reconstruction algorithm')
    # setParameter
	
    def getParameter(self, parameter_name):
        ret = {}
        if parameter_name in self.pars.keys():
            ret[parameter_name] = self.pars[parameter_name]
        else:
            raise Exception('Parameter {0} not currently '.format(parameter_name) +
                            'saved for regularizer algorithm')
    # setParameter    
    def sanityCheck(self):
        projection_data = self.pars['projection_data']
        dark_field = self.pars['dark_field']
        flat_field = self.pars['flat_field']
        angles = self.pars['angles']
        
        if projection_data is not None and dark_field is not None and \
            angles is not None and flat_field is not None:
            data_shape =  numpy.shape(projection_data)
            angle_shape = numpy.shape(angles)
            
            if angle_shape[0] != data_shape[0]:
                #raise Exception('Projections and angles dimensions do not match: %d vs %d' % \
                #                (angle_shape[0] , data_shape[0]) )
                return (False , 'Projections and angles dimensions do not match:'+
                                ' %d vs %d' % (angle_shape[0] , data_shape[0]) )
            
            if data_shape[1:] != numpy.shape(flat_field):
                #raise Exception('Projection and flat field dimensions do not match')
                return (False , 'Projection and flat field dimensions do not match')
            if data_shape[1:] != numpy.shape(dark_field):
                #raise Exception('Projection and dark field dimensions do not match')
                return (False , 'Projection and dark field dimensions do not match')
            
            return (True , '' )
        elif self.pars['normalized_projection_data'] is not None:
            data_shape =  numpy.shape(self.pars['normalized_projection_data'])
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
            
    def reconstruct(self, parameters = None , **kwargs):
        if kwargs is not None:
            self.setParameter(**kwargs)
            
        go , reason = self.sanityCheck()
        if go:
            return self._reconstruct()
        else:
            raise Exception(reason)
            
            
    def _reconstruct(self, parameters=None):
        if parameters is not None:
            self.setPars(parameters)
        parameters = self.pars
        
        if parameters['algorithm'] is not None and \
           parameters['normalized_projection_data'] is not None and \
           parameters['angles'] is not None and \
           parameters['center_of_rotation'] is not None and \
           parameters['iterations'] is not None and \
           parameters['resolution'] is not None and\
           parameters['threads'] is not None and\
           parameters['isLogScale'] is not None:
               
               
           if parameters['algorithm'] in (Reconstructor.Algorithm.CGLS,
                                          Reconstructor.Algorithm.MLEM,
                                          Reconstructor.Algorithm.SIRT,
                                          Reconstructor.Algorithm.CGLS_CONV):
               #store parameters
               
               self.pars = parameters
               if parameters['input_volume'] is not None:
                   if parameters['algorithm'] == Reconstructor.Algorithm.CGLS:
                       stepalg = alg.cgls_step
                   elif parameters['algorithm'] == Reconstructor.Algorithm.SIRT:
                       stepalg = alg.sirt_step
                   elif parameters['algorithm'] == Reconstructor.Algorithm.MLEM:
                       stepalg = alg.mlem_step
                   elif parameters['algorithm'] == Reconstructor.Algorithm.CGLS_CONV:
                       stepalg = alg.cgls_conv_step
                   
                   result = stepalg(
                           parameters['normalized_projection_data'] ,
                           parameters['angles'],
                           parameters['center_of_rotation'],
                           parameters['resolution'],
                           parameters['iterations'],
                           parameters['threads'] ,
                           parameters['isLogScale'],
                           parameters['input_volume']
                           )
               else:
                   result = parameters['algorithm'](
                           parameters['normalized_projection_data'] ,
                           parameters['angles'],
                           parameters['center_of_rotation'],
                           parameters['resolution'],
                           parameters['iterations'],
                           parameters['threads'] ,
                           parameters['isLogScale']
                           )
               self.reconstructedVolume = result
               return result
           else:
               self.pars = parameters
               if parameters['input_volume'] is not None:
                   if parameters['algorithm'] == Reconstructor.Algorithm.CGLS_TIKHONOV:
                       stepalg = alg.cgls_tikhonov_step
                   elif parameters['algorithm'] == Reconstructor.Algorithm.CGLS_TVREG:
                       stepalg = alg.cgls_TVreg_step

                   result = stepalg(
                           parameters['normalized_projection_data'] ,
                           parameters['angles'],
                           parameters['center_of_rotation'],
                           parameters['resolution'],
                           parameters['iterations'],
                           parameters['threads'] ,
                           parameters['regularize'],
                           numpy.zeros((parameters['iterations'])),
                           parameters['isLogScale'],
                           parameters['input_volume']
                           )
               else:
                   result = parameters['algorithm'](
                               parameters['normalized_projection_data'] ,
                               parameters['angles'],
                               parameters['center_of_rotation'],
                               parameters['resolution'],
                               parameters['iterations'],
                               parameters['threads'] ,
                               parameters['regularize'],
                               numpy.zeros((parameters['iterations'])),
                               parameters['isLogScale']
                               )
               self.reconstructedVolume = result
               return result
        else:
           if parameters['projection_data'] is not None and \
                     parameters['dark_field'] is not None and \
                     parameters['flat_field'] is not None:
               norm = self.normalize(parameters['projection_data'],
                                   parameters['dark_field'], 
                                   parameters['flat_field'], 0.1)
               self.pars['normalized_projection_data'] = norm
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
        
    def nextIteration(self, num_of_iterations=1):
        self.setParameter(iterations=num_of_iterations)
        self.setParameter(input_volume = self.reconstructedVolume)
        return self._reconstruct()
        
        
        

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
        
    
    def clipNormalizedProjections(self, projections=None,
                                 greater_than_one=1, smaller_than_zero=0):

        if projections is None:
            projections = self.normalizedProjections
        
        gg = numpy.frompyfunc(ProjectionPreprocessor.greater,2,1)
        ss = numpy.frompyfunc(ProjectionPreprocessor.smaller,2,1)
        #eq = numpy.frompyfunc(ProjectionPreprocessor.equal,2,1)
        
        # set value > 1 to 1
        gt = ss(projections, 1)
        projections = numpy.where(gt, projections, greater_than_one)
        
        # set values < 0 to 0
        lt = gg(projections, 0)
        projections = numpy.where(lt, projections, smaller_than_zero)
        
        
        return projections
    
    @staticmethod    
    def greater(x,y):
        return x>y
    @staticmethod
    def smaller(x,y):
        return x<y
    @staticmethod
    def equal(x,y):
        return x == y
