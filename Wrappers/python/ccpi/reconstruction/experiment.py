# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 09:43:36 2017

@author: ofn77899
"""

class TomographyExperiment(CCPiBaseClass):
    
    def __init__(self, **kwargs):
        self.acceptedInputKeywords = ['projections' , 'normalized_projections',
                                      'angles', 'flat_field', 'dark_field',
                                      'instrument', 'reconstructor']
        
        for key , value in kwargs.items():
            if key in self.acceptedInputKeywords()
                self.setParameter(key=value)
            else:
                print(r'Warning: discarded parameter: "{0}"'.format(key))
    
    def getNormalizedProjections(self, projections=None, flat_field=None, dark_field=None, 
                  angles=None, store=False):
        if projections is None:
            try:
                projections = self.getParameter('projections')
            except Exception:
                raise Exception('Please provide the projection data')
        if flat_field is None:
            try:
                flat_field = self.getParameter('flat_field')
            except Exception:
                raise Exception('Please provide the flat_field data')
        if dark_field is None:
            try:
                dark_field = self.getParameter('dark_field')
            except Exception:
                raise Exception('Please provide the dark_field data')
        
        norm = TomographyExperiment.normalize(projections, 
                                                flat_field, 
                                                dark_field)
        if store:
            self.setParameter(normalized_projections=norm)
        return norm
    
    @staticmethod
    def normalize(projections, flats, darks):
        # 2 is dark field
        #darks = [stack[i] for i in range(len(itype)) if itype[i] == 2 ]
        if len(numpy.shape(darks)) == 3:
            dark = darks[0]
            for i in range(1, len(darks)):
                dark += darks[i]
            dark = dark / len(darks)
        else:
            # supposedly is 2
            dark = darks
        
        
        # 1 is flat field
        #flats = [stack[i] for i in range(len(itype)) if itype[i] == 1 ]
        if len(numpy.shape(flats)) == 3:
            flat = flats[0]
            for i in range(1, len(flats)):
                flat += flats[i]
            flat = flat / len(flats)
        else:
            # supposedly the size is 2: already an image
            dark = darks    
        
        norm = [TomographyExperiment.normalize2D(projection, dark, flat) \
                for projection in projections]
        norm = numpy.asarray (norm, dtype=numpy.float32)
        
        return norm
    
    @staticmethod
    def normalize2D(projection, dark, flat, def_val=0.01):
            a = (projection - dark)
            b = (flat-dark)
            with numpy.errstate(divide='ignore', invalid='ignore'):
                c = numpy.true_divide( a, b )
                c[ ~ numpy.isfinite( c )] = def_val  # set to not zero if 0/0 
            return c
    
        
    
    
    
    