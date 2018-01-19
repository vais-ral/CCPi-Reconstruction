# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 09:53:27 2018

@author: ofn77899
@editor: mnf98541
"""

from ccpi.viewer.CILViewer2D import Converter
import vtk
import numpy
import h5py
# load your image as imagedata


'''Load a dataset stored in a NeXuS file (HDF5)'''
###############################################################################
## Load a dataset
print ("Loading Data")
# replace with your file location
filename = "/home/mnf98541/24737_fd.nxs"
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

imagedata = numpy.transpose(norm, [2,1,0]).copy()

conv = Converter.numpy2vtkImporter(imagedata)
conv.Update()

writer = vtk.vtkMetaImageWriter()
writer.SetInputData(conv.GetOutput())
writer.SetFileName("ImageData.mha")
writer.Write()

print(angle_proj)