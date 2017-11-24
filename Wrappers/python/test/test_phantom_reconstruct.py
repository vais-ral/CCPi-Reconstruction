# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 09:07:26 2017

@author: ofn77899


Test forward projector
"""
from ccpi.reconstruction.parallelbeam import alg
from ccpi.reconstruction.FindCenterOfRotation import *
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
nangles = 60
angles = numpy.linspace(0,180, nangles, dtype=numpy.float32)

phantom = numpy.asarray(phantom, dtype=numpy.float32)

pixel_per_voxel = 1
fp = alg.pb_forward_project
stack = fp(phantom, angles , pixel_per_voxel)
print ("stack " ,stack.min(), stack.max())
#bp = alg.pb_backward_project(stack, angles, center_of_rotation, pixel_per_voxel)


# CGLS
niterations = 30
threads = 4
isPixelDataInLogScale = False
normalize = True
if normalize:
	m = stack.min()
	M = stack.max()
	scale = 1 / (M-m)
	shift = -m 
	norm = stack * scale + shift
	
	#norm = numpy.transpose(norm , axes=[0,2,1])
# absorption is from 0 to 1 where 0 is max density
#norm = 0.99 - norm


if False:
	# add some noise
	perc = 0.01
	norm = norm + (perc* numpy.random.normal(size=numpy.shape(norm)))


# the cgls doesn't like 0's and more the 1s	

print ("norm " ,norm.min(), norm.max())	
invnorm = 1-norm
delta = 1e-3
invnorm [numpy.where(invnorm==0)] = delta
invnorm [numpy.where(invnorm>=1-delta)] = 1-delta

center_of_rotation = numpy.shape(invnorm)[2]/2
center_of_rotation = find_center_vo(invnorm)

print ("invnorm " ,invnorm.min(), invnorm.max())	
#img_cgls = alg.cgls(numpy.transpose(invnorm,[0,1,2]), angles, center_of_rotation , \
#   pixel_per_voxel ,  niterations, threads, isPixelDataInLogScale)

img_cgls = alg.cgls(numpy.transpose(invnorm,[0,2,1]), angles, center_of_rotation , \
   pixel_per_voxel ,  niterations, threads, isPixelDataInLogScale)

#img_cgls2 = alg.cgls(norm, angles, center_of_rotation , \
#   pixel_per_voxel ,  niterations + 4, threads, isPixelDataInLogScale)

cols = 4
rows = 1
current = 0
fig = plt.figure()


current = current + 1
if cols >= current:
	a=fig.add_subplot(rows,cols,current)
	a.set_title('projections')
	imgplot = plt.imshow(stack[20])
	current = current + 1

if cols >= current:
	a=fig.add_subplot(rows,cols,current)
	a.set_title('normalized projections')
	imgplot = plt.imshow(norm[20])
	current = current + 1

if cols >= current:
	a=fig.add_subplot(rows,cols,current)
	a.set_title('inverted norm projections')
	imgplot = plt.imshow(1-norm[20])
	current = current + 1

if cols >= current:
	a=fig.add_subplot(rows,cols,current)
	a.set_title('CGLS iterations {0}'.format(niterations))
	imgplot = plt.imshow(img_cgls[64])

plt.show()

if True:
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


