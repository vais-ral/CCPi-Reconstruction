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
from ccpi.viewer.CILViewer2D import *
import os

os.chdir('C:/Users/ofn77899/Documents/GitHub/CCPi-Reconstruction/Wrappers/Python/demo')

def display(vol):
	v = CILViewer2D()
	v.setInputAsNumpy(vol)
	v.startRenderLoop()
	return v

nx = h5py.File(r'../../../../../CCPi/CIL-Docs/data/phant3D_256.h5', "r")
ph = numpy.asarray(nx.get('/dataset1'))

#phantom = ph[1:254,2:253,3:252]
phantom = ph[10:250,:,20:240].copy()

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
angles = numpy.linspace(0,360, nangles, dtype=numpy.float32)

phantom = numpy.asarray(phantom, dtype=numpy.float32)

pixel_per_voxel = 1
fp = alg.pb_forward_project
stack = fp(phantom, angles , pixel_per_voxel)
print ("stack " ,stack.min(), stack.max())
center_of_rotation = numpy.shape(stack)[2]/2.
print ("cor: {0} {1}".format(numpy.shape(stack)[1]/2., numpy.shape(stack)[2]/2.))
center_of_rotation = find_center_vo(stack)
print ("cor: {0}".format(center_of_rotation))

back = alg.pb_backward_project(stack, angles, center_of_rotation, pixel_per_voxel)



cols = 2
rows = 1
current = 0
fig = plt.figure()


current = current + 1
if cols >= current:
	a=fig.add_subplot(rows,cols,current)
	a.set_title('forward projections')
	imgplot = plt.imshow(stack[20])
	current = current + 1

if cols >= current:
	a=fig.add_subplot(rows,cols,current)
	a.set_title('back projections')
	imgplot = plt.imshow(back[20])
	current = current + 1


plt.show()

if False:
	v = display(back)
	




