# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 17:41:44 2017

@author: ofn77899
"""

from Reconstructor import Reconstructor
from FISTAReconstructor import FISTAReconstructor

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