# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 11:45:35 2017

@author: ofn77899
"""
from abc import ABCMeta, abstractmethod
from enum import Enum


class Reconstructor(metaclass=ABCMeta):
    '''Abstract class that defines the Reconstruction class for reconstruction

This class defines the methods that must be implemented by concrete classes.
Each concrete sub-class is requested to define a Factory to create an instance
of the specific reconstruction.

    '''
    def __init__(self, **kwargs):
        return NotImplemented
    
    @staticmethod
    def create(type, **kwargs):
        if type in Reconstructor.__subclasses__():
            if kwargs == {}:
                return type.Factory.create()
            else:
                return type.Factory.create(**kwargs)
        else:
            raise TypeError(r'Reconstructor "{0}" not available'.format(type))
           
    def setParameter(self, **kwargs):
        '''set named parameter for the reconstructor engine
        
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
                
    class Factory():
        @staticmethod
        def create(**kwargs):
            return NotImplemented 
        
class CGLS(Reconstructor):
    
    class Factory:
        def create(**kwargs):
            if kwargs != {}:
                return CGLS(**kwargs)
            else:
                return CGLS()
    
    def __init__(self, **kwargs):
        self.acceptedInputKeywords = \
               ['algorithm','projection_data' ,'normalized_projection', \
                'angles' , 'center_of_rotation' , 'flat_field', \
                'iterations','dark_field' , 'resolution', 'isLogScale' , \
                'threads' , 'iterationValues', 'regularize']
               
        self.pars = {}
        print ("{0} created".format(__name__))
        if kwargs != {}:
            print (kwargs) 
        
        
        
if __name__ == "__main__":
    r = Reconstructor.create(CGLS)
    r.setParameter(resolution=1)
       