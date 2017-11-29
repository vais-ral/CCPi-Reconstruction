# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 09:43:36 2017

@author: ofn77899
"""

from ccpi.reconstruction.Reconstructor import Reconstructor, IterativeReconstructor
from ccpi.reconstruction.parallelbeam import alg as pbalg
from ccpi.reconstruction.conebeam import alg as cbalg


class TomographyExperiment(CCPiBaseClass):
    
    def __init__(self, **kwargs):
        self.acceptedInputKeywords = ['instrument', 'reconstructor']
        self.available_reconstructors = ['cgls' , 'sirt', 'mlem', 'cgls_conv',
                                         'cgls_tv', 'cgls_tikhonov']
        for key , value in kwargs.items():
            if key in self.acceptedInputKeywords()
                self.setParameter(key=value)
            else:
                print(r'Warning: discarded parameter: "{0}"'.format(key))
    
    def createReconstructor(self, reconstructor_name):
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
        
        reconstructor = Reconstructor.create(algorithm=algorithm)
        self.setParameter(reconstructor=reconstructor)        
        
    def isConeBeam(self):
        instrument = self.getParameter('instrument')
        return instrument.isConeBeam()
    
    def printAvailableReconstructionAlgorithms(self):
        for alg in self.available_reconstructors:
            print (alg)
        
    def configureReconstructor(self, reconstructor=None , **kwargs):
        if reconstructor is None:
            reconstructor = self.getParameter('reconstructor')
        
        # pass evenctual 
        
        instrument= self.getParameter('instrument')
        if instrument.isConeBeam():
            pass
        else:
            normalized_projections = instrument.getNormalizedProjections()
            angles = instrument.getParameter('angles')
            center_of_rotation = instrument.getCenterOfRotation()
        
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
                if kwargs is not {}:
                    if 'regularization_parameter' in kwargs.keys():
                        reconstructor.setParameter(
                            regularization_parameter=kwargs['regularization_parameter']
                            )
                    else:
                        try:
                            # this should rise an exception
                            reg = reconstructor.getParameter('regularization_parameter')
                        except Exception():
                            raise Exception('ERROR: Please set the regularization parameter for the reconstructor')
        
        reconstructor.iterate()
        
        return True
        
    
    
    
    