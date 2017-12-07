# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 11:45:35 2017

@author: ofn77899
"""

#from common import CCPiBaseClass
from ccpi.reconstruction.parallelbeam import alg as pbalg
from ccpi.reconstruction.conebeam import alg as cbalg
from ccpi.reconstruction.FindCenterOfRotation import find_center_vo
from ccpi.common import CCPiBaseClass
import numpy


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
        for algs in IterativeReconstructor.getAvailableAlgorithms():
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
        def create(algorithm, **kwargs):
            if kwargs != {}:
                return IterativeReconstructor(algorithm=algorithm, **kwargs)
            else:
                return IterativeReconstructor(algorithm=algorithm)
        
    
    def __init__(self, **kwargs):
        self.acceptedInputKeywords = \
               ['algorithm','projection_data' ,'normalized_projections', \
                'angles' , 'center_of_rotation' , 'flat_field', \
                'iterations','dark_field' , 'resolution', 'isLogScale' , \
                'threads' , 'iterationValues', 'regularization_parameter']
               
        self.pars = {'algorithm' : pbalg.cgls,
                    'resolution' : 1, 
                     'isLogScale' : False,
                     'threads' : 16, 
                     'iterations' : 8}
        
        print ("{0} created".format(self.getParameter('algorithm').__name__))
        if kwargs != {}:
            print (kwargs) 
            self.setParameter(**kwargs)
            
    @staticmethod
    def getAvailableAlgorithms():
        return ['cgls' , 'sirt', 'mlem', 'cgls_conv', 'cgls_TVreg' , 'cgls_tikhonov']
        
    def iterate(self, **kwargs):
        if kwargs != {}:
            self.setParameter(**kwargs)
        try:
            normalized_projections = self.getParameter('normalized_projections')
        except KeyError:
            raise Exception('Insufficient data. Please pass the normalized projections or flat, dark and projections')
        
        try:
            angles = self.getParameter('angles')
        except KeyError:
            raise Exception('Please pass the angles')
            
        #center of rotation
        try:
            center_of_rotation = self.getParameter('center_of_rotation')
            print ("center of rotation found as {0}".format(center_of_rotation))
        except KeyError:
            # estimate center of rotation
            print ("estimating center of rotation".format(center_of_rotation))
            center_of_rotation = find_center_vo(normalized_projections)
            self.setParameter(center_of_rotation=center_of_rotation)
        
        resolution , niterations, threads, isPixelDataInLogScale = \
            self.getParameter(['resolution' , 'iterations', 
                               'threads', 'isLogScale'])      
        algorithm = self.getParameter('algorithm')
        print ("niterations",niterations,"threads",threads)
        if algorithm.__name__ == "cgls" or \
           algorithm.__name__ == "sirt" or \
           algorithm.__name__ == "mlem":
            # CGLS
            volume = algorithm(normalized_projections, angles, center_of_rotation , 
                              resolution , niterations, threads, 
                              isPixelDataInLogScale)
        elif algorithm.__name__ == "cgls_conv":
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

