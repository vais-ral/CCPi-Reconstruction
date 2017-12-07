# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 09:43:36 2017

@author: ofn77899
"""
import numpy
from ccpi.reconstruction.Reconstructor import IterativeReconstructor
from ccpi.reconstruction.parallelbeam import alg as pbalg
from ccpi.reconstruction.conebeam import alg as cbalg
from ccpi.common import CCPiBaseClass

class TomographyExperiment(CCPiBaseClass):
    
    def __init__(self, **kwargs):
        CCPiBaseClass.__init__(self)
        self.acceptedInputKeywords = ['instrument', 'reconstructor']
        self.available_reconstructors = ['cgls' , 'sirt', 'mlem', 'cgls_conv',
                                         'cgls_tv', 'cgls_tikhonov']
        for key , value in kwargs.items():
            if key in self.acceptedInputKeywords():
                self.setParameter(key=value)
            else:
                print(r'Warning: discarded parameter: "{0}"'.format(key))
    
    def createReconstructor(self, reconstructor_name):
        
        algorithm = self.getAlgorithm(reconstructor_name)
        
        reconstructor = IterativeReconstructor.Factory.create(algorithm)
        self.setParameter(reconstructor=reconstructor)        
    
    
    def getAlgorithm(self, reconstructor_name):
        isConeBeam = self.isConeBeam()
        
        if reconstructor_name == 'cgls':
            if isConeBeam:
                algorithm = cbalg.cgls
            else:
                algorithm = pbalg.cgls
        elif reconstructor_name == 'sirt':
            if isConeBeam:
                algorithm = cbalg.sirt
            else:
                algorithm = pbalg.sirt
        elif reconstructor_name == 'mlem':
            if isConeBeam:
                algorithm = cbalg.mlem
            else:
                algorithm = pbalg.mlem
        elif reconstructor_name == 'cgls_conv':
            if isConeBeam:
                algorithm = cbalg.cgls_conv
            else:
                algorithm = pbalg.cgls_conv
        elif reconstructor_name == 'cgls_tv':
            if isConeBeam:
                algorithm = cbalg.cgls_tv
            else:
                algorithm = pbalg.cgls_tv
        elif reconstructor_name == 'cgls_tikhonov':
            if isConeBeam:
                algorithm = cbalg.cgls_tikhonov
            else:
                algorithm = pbalg.cgls_tikhonov
        return algorithm
    
    def isConeBeam(self):
        instrument = self.getParameter('instrument')
        return instrument.isConeBeam()
    
    def printAvailableReconstructionAlgorithms(self):
        for alg in self.available_reconstructors:
            print (alg)
        
    def configureReconstructor(self, **kwargs):
        '''Configures the reconstructor with the current data present in the experiment
        
        returns True if fully configured, raises exceptions when missing data'''
        
        print ("configureReconstructor {0}".format(kwargs.keys()))

        if 'reconstructor' in kwargs.keys():
            self.setParameter(reconstructor=reconstructor)
            
        reconstructor = self.getParameter('reconstructor')
        
        # pass keyworded arguments
        reconstructor.setParameter(**kwargs)
        
        instrument= self.getParameter('instrument')
        if instrument.isConeBeam():
            pass
        else:
            normalized_projections = instrument.getNormalizedProjections()
            #reconstructor.setParameter(normalized_projection = normalized_projections)
            angles = instrument.getParameter('angles')
            if 'center_of_rotation' in kwargs:
                center_of_rotation = kwargs['center_of_rotation']
            else:
                print ("calling instrument.getCenterOfRotation()")
                center_of_rotation = instrument.getCenterOfRotation()
            
            # pass the data to the reconstructor
            reconstructor.setParameter(normalized_projections = normalized_projections,
                                       angles=angles, 
                                       center_of_rotation=center_of_rotation)
        
        algorithm = reconstructor.getParameter('algorithm')
        
        if algorithm.__name__ == 'cgls' or \
           algorithm.__name__ == 'sirt' or \
           algorithm.__name__ == 'mlem' or \
           algorithm.__name__  == 'cgls_conv':
            if self.isConeBeam():
                pass
            else:
                pass
            
        elif algorithm.__name__  == 'cgls_tv' or \
             algorithm.__name__  == 'cgls_tikhonov':
            if isConeBeam:
                pass
            else:
                try:
                    # this should rise an exception if not configured
                    reg = reconstructor.getParameter('regularization_parameter')
                except Exception():
                    raise Exception('ERROR: Please set the regularization parameter for the reconstructor')
        
        
        return True
        
    def reconstruct(self, algorithm_name=None, 
                    iterations = None, 
                    instrument = None, 
                    regularization_parameter = None,
                    center_of_rotation = None):
        
        
        kwargs = {}
        if algorithm_name is not None:
            algorithm = self.getAlgorithm(algorithm_name)
            kwargs['algorithm']=algorithm
        if iterations is not None:
            #reconstructor.setParameter(iterations=number_of_iterations)
            kwargs['iterations']=iterations
        if regularization_parameter is not None:
            #reconstructor.setParameter(regularization_parameter=regularization_parameter)
            kwargs['regularization_parameter']=regularization_parameter
        if center_of_rotation is not None:
            kwargs['center_of_rotation'] = center_of_rotation
        
        if self.configureReconstructor(**kwargs):
            print ("reconstructor configured correctly")
            reconstructor = self.getParameter('reconstructor')
            check = self.sanityCheck()
            if check[0]:
                return reconstructor.iterate()
            else:
                raise Exception(check[1])
        
            
    def sanityCheck(self):
        instrument = self.getParameter('instrument')
        try:
            projection_data = instrument.getParameter('projections')
            dark_field = instrument.getParameter('dark_field')
            flat_field = instrument.getParameter('flat_field')
        except KeyError:
            norm_proj = instrument.getParameter('normalized_projections')
            projection_data = None
            dark_field = None
            flat_field = None
        angles = instrument.getParameter('angles')
        
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
        elif norm_proj is not None:
            data_shape =  numpy.shape(norm_proj)
            angle_shape = numpy.shape(angles)
            
            print (data_shape)          
            if angle_shape[0] != data_shape[0]:
                #raise Exception('Projections and angles dimensions do not match: %d vs %d' % \
                #                (angle_shape[0] , data_shape[0]) )
                return (False , 'Projections and angles dimensions do not match: %d vs %d' % \
                                (angle_shape[0] , data_shape[0]) )
            else:
                return (True , '' )
        else:
            return (False , 'Not enough data')       
            
    
    
    
    
