# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 09:43:36 2017

@author: ofn77899
"""

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
    