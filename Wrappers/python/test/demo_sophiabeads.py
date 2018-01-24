#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import ccpi.instrument as CI
import astra
import matplotlib.pyplot as pp
import numpy as np

# Set up Xtek instrument.
x = CI.Xtek()

# Load SophiaBeads dataset with 256 projections. Update path for your system.
x.read('/media/jakob/050d8d45-fab3-4285-935f-260e6c5f162c1/Data/SophiaBeads/SophiaBeads_256_averaged/SophiaBeads_256_averaged.xtekct')

# Show a projection as image
pp.imshow(x.pixels[0,])
pp.show()

# Function to convert to ASTRA format for 2D central slice fan beam reconstruction
def Xtek2astra2d(x,pad_size=0):
    
    # number of detector elements (bins)
    numbins = x.num_of_horizontal_pixels
    
    #  % number of pixels on each side of object Should really be geom.voxels(1); 
    N = numbins 
    
    # Geometric magnification
    geomag = np.abs(x.detector_x/x.source_x)
    
    # Physical sizes (mm) of detector pixel and object voxel
    pixel_size_physical = np.abs(x.h_pixels[1] - x.h_pixels[0])
    voxel_size_physical = pixel_size_physical / geomag
    
    # Relative size of detector pixel to object voxel
    det_width = geomag  # i.e. pixel_size_physical / voxel_size_physical;
    
    # Number of detector pixels
    det_count = x.num_of_horizontal_pixels
    
    # The angles, negated and converted to radians
    angles = -x.angles * np.pi/180;
    
    # Source to centre of rotation distance in units of object pixels
    source_origin = abs(x.source_x) / voxel_size_physical;
    
    # Centre of rotation to detector distance in units of object pixels
    origin_det = abs(x.detector_x - x.source_x) / voxel_size_physical;
    
    # Extract center slice of sinogram    
    DATA = x.pixels[:,:,1000]
    
    # Pad either on left or right hand side as simple centering correction
    if pad_size > 0:
        DATA = np.hstack( (np.tile(DATA[:,0], (pad_size,1)).transpose(), \
                       DATA ) )
    else:
        DATA = np.hstack( (DATA, \
                       np.tile(DATA[:,-1], (-pad_size,1)).transpose() ) )
    
    # Update number of detector pixels to included padding.
    det_count_pad = det_count+np.abs(pad_size);
    
    # Set up ASTRA volume and projection geometries
    vol_geom = astra.create_vol_geom(N);
    proj_geom = astra.create_proj_geom('fanflat', det_width, det_count_pad, \
        angles, source_origin, origin_det)
    
    proj_id = astra.create_projector('line_fanflat', proj_geom, vol_geom)
    
    return DATA, vol_geom, proj_geom, proj_id

# Run function to extract center slice data and ASTRA geometries    
DATA, vol_geom, proj_geom, proj_id = Xtek2astra2d(x,pad_size=30)

# Do ASTRA FBP reconstruction on CPU
rec_id1, rec1 = astra.create_reconstruction('FBP',proj_id,DATA)

# Do ASTRA FBP reconstruction on GPU
rec_id2, rec2 = astra.create_reconstruction('FBP_CUDA',proj_id,DATA)

# Do ASTRA CGLS reconstruction on GPU
rec_id3, rec3 = astra.create_reconstruction('CGLS_CUDA',proj_id,DATA,iterations=20)

# Display image and zooms
pp.imshow(rec1)
pp.show()

pp.imshow(rec1[:1000,:1000])
pp.show()

pp.imshow(rec2)
pp.show()

pp.imshow(rec2[:1000,:1000])
pp.show()

pp.imshow(rec3)
pp.show()

pp.imshow(rec3[:1000,:1000])
pp.show()

# Function to convert to ASTRA format for 3D cone beam reconstruction
# Most of this is the same as in the 2D version and perhaps could be merged.
def Xtek2astra3d(x,pad_size=0):
    
    # number of detector elements (bins)
    numbins = x.num_of_horizontal_pixels
    
    #  % number of pixels on each side of object Should really be geom.voxels(1); 
    N = numbins 
    
    # Geometric magnification
    geomag = np.abs(x.detector_x/x.source_x)
    
    # Physical sizes (mm) of detector pixel and object voxel
    pixel_size_physical = np.abs(x.h_pixels[1] - x.h_pixels[0])
    voxel_size_physical = pixel_size_physical / geomag
    
    # Relative size of detector pixel to object voxel
    det_width = geomag  # i.e. pixel_size_physical / voxel_size_physical;
    
    # Number of detector pixels
    det_count_horz = x.num_of_horizontal_pixels
    
    # The angles, negated and converted to radians
    angles = -x.angles * np.pi/180;
    
    # Source to centre of rotation distance in units of object pixels
    source_origin = abs(x.source_x) / voxel_size_physical;
    
    # Centre of rotation to detector distance in units of object pixels
    origin_det = abs(x.detector_x - x.source_x) / voxel_size_physical;
    
    # Number of slices in z-direction
    slices = x.num_of_vertical_pixels
    
    # Permute to ASTRA order: row (v), angle, column (u) 
    # from the instrument order which is angles, column, row.
    DATA = np.transpose(x.pixels, (2, 0, 1))
    
    # Pad either on left or right hand side as simple centering correction
    if pad_size > 0:
        padvals = DATA[:,:,0]
        padvals.shape = (padvals.shape[0],padvals.shape[1],1)
        DATA = np.dstack( (np.tile(padvals,(1,1,pad_size)),DATA))
    else:
        padvals = DATA[:,:,-1]
        padvals.shape = (padvals.shape[0],padvals.shape[1],1)
        DATA = np.dstack( (np.tile(DATA,padvals,(1,1,pad_size))))
    
    # Update number of detector pixels to included padding.
    det_count_horz_pad = det_count_horz+np.abs(pad_size);
    
    # Set up ASTRA volume and projection geometries
    vol_geom = astra.create_vol_geom(N, N, slices);
    proj_geom = astra.create_proj_geom('cone', det_width, det_width, \
        slices, det_count_horz_pad, angles, source_origin, origin_det)
    
    return DATA, vol_geom, proj_geom

# Run function to convert 3D data to ASTRA data format and geometries
DATA3, vol_geom3, proj_geom3 = Xtek2astra3d(x,pad_size=30)

data_id = astra.data3d.create('-proj3d', proj_geom3, DATA3)

# Create a data object for the reconstruction
rec_id = astra.data3d.create('-vol', vol_geom3)

# Set up and run ASTRA FDK on GPU and retrieve result.
cfg = astra.astra_dict('FDK_CUDA')
cfg['ReconstructionDataId'] = rec_id
cfg['ProjectionDataId'] = data_id
alg_id = astra.algorithm.create(cfg)
astra.algorithm.run(alg_id)
rec_3d = astra.data3d.get(rec_id)


pp.imshow(rec_3d[1000,:,:])
pp.show()

pp.imshow(rec_3d[:,1000,:])
pp.show()

pp.imshow(rec_3d[:,:,1000])
pp.show()

# Clean up. Note that GPU memory is tied up in the algorithm object,
# and main RAM in the data objects.
astra.algorithm.delete(alg_id)
astra.data3d.delete(rec_id)
astra.data3d.delete(rec_id1)
astra.data3d.delete(rec_id2)
astra.data3d.delete(rec_id3)
astra.data3d.delete(data_id)