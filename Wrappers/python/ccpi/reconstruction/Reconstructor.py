# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 11:45:35 2017

@author: ofn77899
"""

#from common import CCPiBaseClass
from ccpi.reconstruction.parallelbeam import alg as pbalg
from ccpi.reconstruction.conebeam import alg as cbalg
from ccpi.reconstruction.FindCenterOfRotation import find_center_vo
import numpy

class CCPiBaseClass(object):
    def __init__(self, **kwargs):
        self.acceptedInputKeywords = []
        self.pars = {}
    
    def setParameter(self, **kwargs):
        '''set named parameter for the reconstructor engine
        
        raises Exception if the named parameter is not recognized
        
        '''
        for key , value in kwargs.items():
            if key in self.acceptedInputKeywords:
                self.pars[key] = value
            else:
                raise Exception('Wrong parameter "{0}" for {1}'.format(key, self.__class__.__name__ ))
    # setParameter

    def getParameter(self, key):
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
    #getParameter
    

class Reconstructor(CCPiBaseClass):
    '''Abstract class that defines the Reconstruction class for reconstruction

This class defines the methods that must be implemented by concrete classes.
Each concrete sub-class is requested to define a Factory to create an instance
of the specific reconstruction.

    '''
    def __init__(self, **kwargs):
        self.acceptedInputKeywords = ['instrument']
    
    @staticmethod
    def create(type, **kwargs):
        if type in Reconstructor.__subclasses__():
            if kwargs == {}:
                return type.Factory.create()
            else:
                return type.Factory.create(**kwargs)
        else:   
            if type in IterativeReconstructor.__subclasses__():
                print ("found type {0}".format(type.__name__))
                return type.Factory.create(**kwargs)
            raise TypeError(r'Reconstructor "{0}" not available'.format(type))
           
        
    
    def iterate(self, **kwargs):
        return NotImplementedError('iterate must be implemented in the concrete class')
    
    @staticmethod
    def printAvailableAlgorithms():
        msg = 'Available Iterative Reconstruction Algorithms:\n'
        for algs in Reconstructor.__subclasses__():
            if algs != IterativeReconstructor:
                msg = msg + '{0}, instantiate as Reconstructor.create({0}, **kwargs)\n'\
                    .format(algs.__name__)
        for algs in IterativeReconstructor.__subclasses__():
            msg = msg + '{0}, instantiate as Reconstructor.create({0}, **kwargs)\n'\
                .format(algs.__name__)
        return msg
    
    
                
    class Factory():
        @staticmethod
        def create(**kwargs):
            return NotImplemented 
        
class IterativeReconstructor(Reconstructor):
    
    class Factory:
        @staticmethod
        def create(**kwargs):
            if kwargs != {}:
                return IterativeReconstructor(algorithm=pbalg.cgls, **kwargs)
            else:
                return IterativeReconstructor(algorithm=pbalg.cgls)
        
    @staticmethod
    def create(type, **kwargs):
        if type in IterativeReconstructor.__subclasses__():
            if kwargs == {}:
                return type.Factory.create()
            else:
                return type.Factory.create(**kwargs)
        else:
            raise TypeError(r'Reconstructor "{0}" not available'.format(type))
    
    def __init__(self, **kwargs):
        self.acceptedInputKeywords = \
               ['algorithm','projection_data' ,'normalized_projection', \
                'angles' , 'center_of_rotation' , 'flat_field', \
                'iterations','dark_field' , 'resolution', 'isLogScale' , \
                'threads' , 'iterationValues', 'regularization_parameter']
               
        self.pars = {'algorithm' : pbalg.cgls,
                    'resolution' : 1, 
                     'isLogScale' : False,
                     'threads' : 16, 
                     'iterations' : 8}
        
        print ("{0} created".format(__name__))
        if kwargs != {}:
            print (kwargs) 
            self.setParameter(**kwargs)
        
    def iterate(self, **kwargs):
        if kwargs != {}:
            self.setParameter(**kwargs)
        try:
            normalized_projections = self.getParameter('normalized_projection')
        except KeyError:
            raise Exception('Insufficient data. Please pass the normalized projections or flat, dark and projections')
        
        try:
            angles = self.getParameter('angles')
        except KeyError:
            raise Exception('Please pass the angles')
            
        #center of rotation
        try:
            center_of_rotation = self.getParameter('center_of_rotation')
        except KeyError:
            # estimate center of rotation
            center_of_rotation = find_center_vo(normalized_projections)
            self.setParameter(center_of_rotation=center_of_rotation)
        
        resolution , niterations, threads, isPixelDataInLogScale = \
            self.getParameter(['resolution' , 'niterations', 
                               'threads', 'isPixelDataInLogScale'])      
        algorithm = self.getParameter('algorithm')
        
        if algorithm._name__ == "cgls" or \
           algorithm._name__ == "sirt" or \
           algorithm._name__ == "mlem":
            # CGLS
            volume = algorithm(normalized_projections, angles, center_of_rotation , 
                              resolution , niterations, threads, 
                              isPixelDataInLogScale)
        elif algorithm._name__ == "cgls_conv":
            iteration_values1 = numpy.zeros((niterations,))
            volume = algorithm(normalized_projections, angles, 
                                      center_of_rotation , resolution , 
                              niterations , threads,
                              iteration_values1 , isPixelDataInLogScale)
        elif algorithm.__name__ == "cgls_tikhonov" or \
             algorithm.__name__ == "cgls_TVreg":
            iteration_values1 = numpy.zeros((niterations,))
            regularization_parameter = self.getParameter('regularization_parameter')
            volume = algorithm(normalized_projections, angles,
                               center_of_rotation , resolution , 
                               niterations, threads,
                               regularization_parameter, iteration_values1 , 
                               isPixelDataInLogScale)

            
        return volume

class CGLS(IterativeReconstructor):
    class Factory:
        @staticmethod
        def create(**kwargs):
            if kwargs != {}:
                return IterativeReconstructor(algorithm=pbalg.cgls, **kwargs)
            else:
                return IterativeReconstructor(algorithm=pbalg.cgls)
            
class SIRT(IterativeReconstructor):
    class Factory:
        @staticmethod
        def create(**kwargs):
            if kwargs != {}:
                return IterativeReconstructor(algorithm=pbalg.sirt, **kwargs)
            else:
                return IterativeReconstructor(algorithm=pbalg.sirt)

class MLEM(IterativeReconstructor):
    class Factory:
        @staticmethod
        def create(**kwargs):
            if kwargs != {}:
                return IterativeReconstructor(algorithm=pbalg.sirt, **kwargs)
            else:
                return IterativeReconstructor(algorithm=pbalg.sirt)
        

class CGLS_CONV(IterativeReconstructor):
    class Factory:
        @staticmethod
        def create(**kwargs):
            if kwargs != {}:
                return IterativeReconstructor(algorithm=pbalg.cgls_conv, **kwargs)
            else:
                return IterativeReconstructor(algorithm=pbalg.cgls_conv)
            
            
        
            
## CGLS CONV
#iteration_values1 = numpy.zeros((niterations,))
#img_cgls_conv = alg.cgls_conv(norm, angle_proj, center_of_rotation , 
#                              resolution , 
#                              niterations , threads,
#                              iteration_values1 , isPixelDataInLogScale)
#
##Regularization parameter
#regularization = numpy.double(1e-3)
#
## CGLS TIKHONOV
#iteration_values2 = numpy.zeros((niterations,))
#img_cgls_tikhonov = alg.cgls_tikhonov(norm, angle_proj, center_of_rotation , 
#                                      resolution , niterations, threads,
#                                      regularization, iteration_values2 , 
#                                      isPixelDataInLogScale)
#
## CGLS Total Variation Regularization 
#iteration_values3 = numpy.zeros((niterations,))
#img_cgls_TVreg = alg.cgls_TVreg(norm, angle_proj, center_of_rotation , 
#                                resolution ,  niterations, threads,
#                                      regularization, iteration_values3,
#                                      isPixelDataInLogScale)
if __name__ == "__main__":
    r = Reconstructor.create(CGLS, resolution=1, center_of_rotation=9.3)
    #r.setParameter(resolution=1)
    
#    instrument = Diamond()
#    instrument.readNexus()
#    instrument.setParameter(projections, angles, flats, darks)
#    instrument.getNormalizedData()
#    
#    r = Reconstructor.create('cgls', instrument)
#    volume= r.iterate(allpars)
#    
       