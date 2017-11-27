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
#from ccpi.viewer.CILViewer import CILViewer
from ccpi.viewer.CILViewer2D import *
#from ccpi.segmentation.SimpleflexSegmentor import SimpleflexSegmentor

from ccpi.reconstruction.parallelbeam import alg
from ccpi.reconstruction.FindCenterOfRotation import *
#import vtk
import matplotlib.pyplot as plt
from ccpi.instrument import Diamond

def load_data(filename):
    '''Load a dataset stored in a NeXuS file (HDF5)'''
    ###############################################################################
    ## Load a dataset
    print ("Loading Data")
    nx = h5py.File(filename, "r")
    
    data = nx.get('entry1/tomo_entry/data/rotation_angle')
    angles = numpy.zeros(data.shape)
    data.read_direct(angles)
    print (angles)
    # angles should be in degrees
    
    data = nx.get('entry1/tomo_entry/data/data')
    stack = numpy.zeros(data.shape)
    data.read_direct(stack)
    print (data.shape)
    
    print ("Data Loaded")
    
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
    
    return norm, angle_proj


def normalize(stack):
	m = stack.min()
	M = stack.max()
	scale = 1 / (M-m)
	shift = -m 
	print ("m,M,scale,shift" , m,M,scale,shift)
	return stack * scale + shift

###############################################################################
## 1) load a dataset:

# This dataset is freely available at
# https://github.com/DiamondLightSource/Savu/blob/master/test_data/data/24737_fd.nxs 
    
filename = "C:\\Users\\ofn77899\\Documents\\CCPi\\CGLS\\24737_fd_2.nxs"
norm, angle_proj = load_data(filename)
norm = numpy.ascontiguousarray(norm)

###############################################################################
## 2) 
## 
## Data can now be passed to the reconstruction algorithms:
## CGLS, MLEM, SIRT, CGLS_CONV, CGLS_TIKHONOV, CGLS_TVregularization

#center of rotation
center_of_rotation = numpy.double(86.2)
#resolution
resolution = 1
# number of iterations
niterations = 10
# number of threads
threads = 4
#data are in log scale?
isPixelDataInLogScale = False


# CGLS
img_cgls = alg.cgls(norm, angle_proj, center_of_rotation , resolution , 
                    niterations, threads, isPixelDataInLogScale)


					
#img_cgls = numpy.transpose(img_cgls, axes=[0,2,1])					
# check if padded
shape = numpy.shape(img_cgls)
img_cgls = img_cgls[:,:shape[1]-1,:]
#counter = shape[2] - 1
#while True:
#	if numpy.all(img_cgls[:,:,counter] == 0):
#		counter = counter -1 
#	else:
#		break
#
#if counter < shape[2] - 1:
#	img_cgls = img_cgls[:,:,:counter]

#img_cgls = img_cgls[:160]

## MLEM
#img_mlem = alg.mlem(norm, angle_proj,  center_of_rotation , resolution , 
#                    niterations, threads, isPixelDataInLogScale)
## SIRT
#img_sirt = alg.sirt(norm, angle_proj, center_of_rotation , resolution ,  
#                    niterations, threads, isPixelDataInLogScale)
#
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
#
if False:	
	m = img_cgls.min()
	M = img_cgls.max()
	scale = 1 / (M-m)
	shift = -m 
	print ("rescale CGLS m,M,scale,shift" , m,M,scale,shift)
	img_cgls = img_cgls * scale + shift

fp = alg.pb_forward_project
pixel_per_voxel = resolution

diamond = Diamond()

#stack = diamond.doForwardProject(img_cgls, angles=angle_proj, pixel_per_voxel)
#stack = fp(img_cgls, angle_proj , pixel_per_voxel)
#center_of_rotation = find_center_vo(stack)

#print ("stack: ", numpy.shape(stack), stack.min(), stack.max())


if False:
	m = stack.min()
	M = stack.max()
	scale = 1 / (M-m)
	shift = -m 
	print ("m,M,scale,shift" , m,M,scale,shift)
	stack = stack * scale + shift
	
	# negative the projections
	invnorm = 1-stack

invnorm = diamond.doForwardProject(img_cgls, angles=angle_proj, pixel_per_voxel=								pixel_per_voxel, 
								negative=True , normalized=True)

# the cgls doesn't like 0's and more the 1s	
invnorm [numpy.where(invnorm<=0)] = 1e-3
invnorm [numpy.where(invnorm>1)] = 1 - 1e-3

print ("Orig projections " , numpy.shape(norm))
print ("Projected volume " , numpy.shape(invnorm))
#center_of_rotation = numpy.shape(invnorm)[2]/2.
center_of_rotation = diamond.find_center_of_rotation()
niterations = 5

img_cgls2 = alg.cgls(invnorm, angle_proj, center_of_rotation , resolution , 
                    niterations, threads, isPixelDataInLogScale)

# crop
cinvnorm = invnorm[:,37:183:]
tinvnorm = numpy.ascontiguousarray(numpy.transpose(invnorm, [0,2,1]) )

#center_of_rotation = numpy.shape(tinvnorm)[2]/2
niterations = 2
center_of_rotation = find_center_vo(tinvnorm)
img_cgls3 = alg.cgls(tinvnorm, angle_proj, center_of_rotation , resolution , 
                    niterations, threads, isPixelDataInLogScale)

no=60

cols = 3
rows = 2
current = 1
fig = plt.figure()
# projections row
a=fig.add_subplot(rows,cols,current)
a.set_title('exp projections \n{0}'.format(numpy.shape(norm)))
imgplot = plt.imshow(norm[no])

current = current + 1
a=fig.add_subplot(rows,cols,current)
a.set_title('norm Forward-Projections \n{0}'.format(numpy.shape(invnorm)))
imgplot = plt.imshow(invnorm[no])

current = current + 1
a=fig.add_subplot(rows,cols,current)
a.set_title('transposed norm neg FP\n{0}'.format(numpy.shape(tinvnorm)))
imgplot = plt.imshow(tinvnorm[no])

# volume slices row
current = current + 1
a=fig.add_subplot(rows,cols,current)
a.set_title('reconstructed volume \n{0}, slice {1}'.format(numpy.shape(img_cgls),80))
imgplot = plt.imshow(img_cgls[80])

current = current + 1
a=fig.add_subplot(rows,cols,current)
a.set_title('reconstructed vol norm negative \n{0}'.format(numpy.shape(img_cgls2)))
imgplot = plt.imshow(img_cgls2[80])

#sl = stack[no][37:183].T
#sl = normalize(sl)
#sl = tinvnorm[no]


current = current + 1
a=fig.add_subplot(rows,cols,current)
a.set_title('reconstructed vol norm neg trans fp\n {0}'.format(numpy.shape(img_cgls3)))
imgplot = plt.imshow(img_cgls3[80])
plt.show()

def showInteractive(vol):
	v = CILViewer2D()
	v.setInputAsNumpy(vol)
	v.startRenderLoop()