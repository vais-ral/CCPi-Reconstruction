# -*- coding: utf-8 -*-
#   This work is part of the Core Imaging Library developed by
#   Visual Analytics and Imaging System Group of the Science Technology
#   Facilities Council, STFC
#  
#   Copyright 2017 Edoardo Pasca
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

import numpy
import h5py
from ccpi.reconstruction.parallelbeam import alg

def load_data(filename):
    '''Load a dataset stored in a NeXuS file (HDF5)'''
    ###############################################################################
    ## Load a dataset
    print ("Loading Data")
    #fname = "C:\\Users\\ofn77899\\Documents\\CCPi\\CGLS\\24737_fd_2.nxs"
    nx = h5py.File(filename, "r")
    
    data = nx.get('entry1/tomo_entry/data/rotation_angle')
    angles = numpy.zeros(data.shape)
    data.read_direct(angles)
    #print (angles)
    # angles should be in degrees
    
    data = nx.get('entry1/tomo_entry/data/data')
    stack = numpy.zeros(data.shape)
    data.read_direct(stack)
    #print (data.shape)
    
    
    ## Data should be in the range 0-1
    ## in this case we will perform a simple normalization between full field (flat)
    ## and dark field:
    ## norm = (projection - dark)/(flat-dark)
    
    ##
    # Normalize
    data = nx.get('entry1/tomo_entry/instrument/detector/image_key')
    itype = numpy.zeros(data.shape)
    data.read_direct(itype)
    # 2 is dark field
    darks = [stack[i] for i in range(len(itype)) if itype[i] == 2 ]
    dark = darks[0]
    for i in range(1, len(darks)):
        dark += darks[i]
    dark = dark / len(darks)
    #dark[0][0] = dark[0][1]
    
    # 1 is flat field
    flats = [stack[i] for i in range(len(itype)) if itype[i] == 1 ]
    flat = flats[0]
    for i in range(1, len(flats)):
        flat += flats[i]
    flat = flat / len(flats)
    #flat[0][0] = dark[0][1]
    
    
    # 0 is projection data
    proj = [stack[i] for i in range(len(itype)) if itype[i] == 0 ]
    angle_proj = [angles[i] for i in range(len(itype)) if itype[i] == 0 ]
    angle_proj = numpy.asarray (angle_proj)
    angle_proj = angle_proj.astype(numpy.float32)
    
    
    def normalize(projection, dark, flat, def_val=0.1):
        a = (projection - dark)
        b = (flat-dark)
        with numpy.errstate(divide='ignore', invalid='ignore'):
            c = numpy.true_divide( a, b )
            c[ ~ numpy.isfinite( c )] = def_val  # set to not zero if 0/0 
        return c
        
    
    norm = [normalize(projection, dark, flat) for projection in proj]
    norm = numpy.asarray (norm)
    norm = norm.astype(numpy.float32)
    print ("Data Loaded")
    
    return norm, angle_proj



###############################################################################
## 1) load a dataset:

# This dataset is freely available at
# https://github.com/DiamondLightSource/Savu/blob/master/test_data/data/24737_fd.nxs 
    
filename = "C:\\Users\\ofn77899\\Documents\\CCPi\\CGLS\\24737_fd_2.nxs"
norm, angle_proj = load_data(filename)

###############################################################################
## 2) 
## 
## Data can now be passed to the reconstruction algorithms:
## CGLS, MLEM, SIRT, CGLS_CONV, CGLS_TIKHONOV, CGLS_TVregularization

# number of iterations
niterations = 35
# number of threads
threads = 3

# CGLS
img_cgls = alg.cgls(norm, angle_proj, numpy.double(86.2), 1 , niterations, threads, False)
## MLEM
#img_mlem = alg.mlem(norm, angle_proj, numpy.double(86.2), 1 , niterations, threads, False)
## SIRT
#img_sirt = alg.sirt(norm, angle_proj, numpy.double(86.2), 1 , niterations, threads, False)
#
## CGLS CONV
#iteration_values = numpy.zeros((niterations,))
#img_cgls_conv = alg.cgls_conv(norm, angle_proj, numpy.double(86.2), 1 , niterations, threads,
#                              iteration_values, False)
#
## CGLS TIKHONOV
#iteration_values = numpy.zeros((niterations,))
#img_cgls_tikhonov = alg.cgls_tikhonov(norm, angle_proj, numpy.double(86.2), 1 , niterations, threads,
#                                      numpy.double(1e-5), iteration_values , False)
#
## CGLS Total Variation Regularization 
#iteration_values = numpy.zeros((niterations,))
#img_cgls_TVreg = alg.cgls_TVreg(norm, angle_proj, numpy.double(86.2), 1 , niterations, threads,
#                                      numpy.double(1e-5), iteration_values , False)
#
#
#
################################################################################
### 3) Visualize a slice of the reconstructed images 
#
#import matplotlib.pyplot as plt
#fs = 10
#fig, ax = plt.subplots(1,6,sharey=True)
#ax[0].imshow(img_cgls[80])
#ax[0].axis('off')  # clear x- and y-axes
#ax[0].set_title("CGLS" , fontsize = fs)
#
#ax[1].imshow(img_sirt[80])
#ax[1].axis('off')  # clear x- and y-axes
#ax[1].set_title("SIRT" , fontsize = fs)
#
#ax[2].imshow(img_mlem[80])
#ax[2].axis('off')  # clear x- and y-axesplt.show()
#ax[2].set_title("MLEM" , fontsize = fs)
#
#ax[3].imshow(img_cgls_conv[80])
#ax[3].axis('off')  # clear x- and y-axesplt.show()
#ax[3].set_title("CGLS CONV" , fontsize = fs)
#
#ax[4].imshow(img_cgls_tikhonov[80])
#ax[4].axis('off')  # clear x- and y-axesplt.show()
#ax[4].set_title("Tikhonov" , fontsize = fs)
#
#ax[5].imshow(img_cgls_TVreg[80])
#ax[5].axis('off')  # clear x- and y-axesplt.show()
#ax[5].set_title("TV Reg" , fontsize = fs)
#plt.show()
#
