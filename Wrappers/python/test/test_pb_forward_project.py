# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 09:07:26 2017

@author: ofn77899


Test forward projector
"""
from ccpi.reconstruction.parallelbeam import alg
import h5py
import numpy
import matplotlib.pyplot as plt

nx = h5py.File(r'../../../data/phant3D_256.h5', "r")
ph = numpy.asarray(nx.get('/dataset1'))

#phantom = ph[1:254,2:253,3:252]
phantom = ph.copy()

pX,pY,pZ = numpy.shape(phantom)

#filename = r'/home/ofn77899/Reconstruction/CCPi-FISTA_Reconstruction/demos/DendrData.h5'
#nxa = h5py.File(filename, "r")
##getEntry(nx, '/')
## I have exported the entries as children of /
#entries = [entry for entry in nxa['/'].keys()]
#print (entries)
#
#angles_rad = numpy.asarray(nxa.get('/angles_rad'), dtype="float32")


#device = AstraDevice(
#    DeviceModel.DeviceType.PARALLEL3D.value,
#    [ pX , pY , 1., 1., angles_rad],
#    [ pX, pY, pZ ] )
#
#
#proj = device.doForwardProject(phantom)
#stack = [proj[:,i,:] for i in range(len(angles_rad))]
#stack = numpy.asarray(stack)

#angles = numpy.asarray([i / 3.1415 / 10 for i in range(10)], dtype=numpy.float32)
angles = numpy.linspace(0,180,60, dtype=numpy.float32)

center_of_rotation = 128
pixel_per_voxel = 1
fp = alg.pb_forward_project
stack = fp(phantom, angles , center_of_rotation, pixel_per_voxel)



## and backward
bp = alg.pb_backward_project
backvol = bp(stack, angles, center_of_rotation, pixel_per_voxel)


# CGLS
niterations = 10
threads = 4
isPixelDataInLogScale = False
m = stack.min()
M = stack.max()
scale = 1 / (M-m)
shift = -m + 1e-2
norm = stack * scale + shift
# absorption is from 0 to 1 where 0 is max density
norm = 0.99 - norm
norm = numpy.transpose(norm, axes=[0,2,1]).copy()


img_cgls = alg.cgls(norm, angles, center_of_rotation , \
   pixel_per_voxel ,  niterations, threads, isPixelDataInLogScale)

cols = 4
rows = 1
current = 1
fig = plt.figure()
a=fig.add_subplot(rows,cols,current)
a.set_title('volume')
imgplot = plt.imshow(phantom[64])

current = current + 1
a=fig.add_subplot(rows,cols,current)
a.set_title('stack')
imgplot = plt.imshow(stack[20])

current = current + 1
a=fig.add_subplot(rows,cols,current)
a.set_title('norm')
imgplot = plt.imshow(norm[20])

current = current + 1
a=fig.add_subplot(rows,cols,current)
a.set_title('cgls')
imgplot = plt.imshow(img_cgls[64])


plt.show()

from ccpi.viewer.CILViewer2D import *

v = CILViewer2D()
v.setInputAsNumpy(img_cgls)
v.startRenderLoop()

# pf = h5py.File("phantom3D256_projections.h5" , "w")
# pf.create_dataset("/projections", data=stack)
# #pf.create_dataset("/sinogram", data=proj)
# pf.create_dataset("/angles", data=angles)
# #pf.create_dataset("/reconstruction_volume" , data=numpy.asarray([pX, pY, pZ]))
# pf.create_dataset("/camera/size" , data=numpy.asarray([pX , pY ]))
# pf.create_dataset("/camera/spacing" , data=numpy.asarray([1.,1.]))
# pf.flush()
# pf.close()


