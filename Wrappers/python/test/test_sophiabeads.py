# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import os
from ccpi.instrument import Xtek
from ccpi.reconstruction.conebeam import alg as cbalg
import matplotlib.pyplot as plt
import numpy

filename = os.path.join("/" , "media" , "scratch" , "Data" , 
		"SophiaBeads" , "SophiaBeads_64_averaged" , 
		"SophiaBeads_64_averaged.xtekct")

instr = Xtek()
instr.read(filename)
print (instr.angles)

fig = plt.figure()
a = fig.add_subplot(2,1,1)
plt.imshow(instr.pixels[10])

# reduce size
#1564x1564x n of slices
pixels_nY = 600;
pixels_nZ = 600;
size = numpy.shape(instr.pixels)
nY = int ((size[1] - pixels_nY)/2) #% Cut this many pixels in the Y direction.
nZ = int ((size[2] - pixels_nZ)/2) #% Cut this many pixels in the Z direction.
cutdown3D = instr.pixels[:,nY:nY+pixels_nY,nZ:nZ+pixels_nZ]




a = fig.add_subplot(2,1,2)
plt.imshow(-numpy.log(cutdown3D[10]/cutdown3D.max()))
plt.show()

# cutdown the pixel data and rescale to 0:1 (no shift)

#cutdown3D = [instr.pixels[i][nY:nY+pixels_nY,nZ:nZ+pixels_nZ] for i in range(size[0])]



cutdown3D = numpy.asarray(cutdown3D)

#cmax = cutdown3D.max()
#cutdown3D = cutdown3D / cmax

print (numpy.shape(cutdown3D))

h_offsets = numpy.zeros((1,), dtype='float32')
v_offsets = numpy.zeros((1,), dtype='float32')

beam_harden = False
full_vox_origin = numpy.zeros(3, dtype='float32')
voxel_size = numpy.ones(3,dtype='float32')
num_iterations = 5
nthreads = 10
isPixelInLog = True

print ("launching CGLS")

recon = cbalg.cgls(
		-numpy.log(cutdown3D/cutdown3D.max()), 
		instr.angles, 
		h_offsets,
		v_offsets,
		instr.pixels_per_voxel,
		instr.source_x,
		instr.detector_x,
		instr.pixel_h_size,
		instr.pixel_v_size,
		instr.m_radius,
		beam_harden,
		full_vox_origin,
		voxel_size,
		num_iterations,
		nthreads,
		isPixelInLog
		)
