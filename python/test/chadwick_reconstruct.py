# -*- coding: utf-8 -*-
"""
Created on Fri Jun  2 15:48:00 2017

@author: ofn77899
"""

import numpy
#import h5py
from ccpi.reconstruction.parallelbeam import alg
from ccpi.reconstruction.IterativeReconstruction import Reconstructor, ProjectionPreprocessor

from enum import Enum
from ccpi.viewer.CILViewer2D import Converter, CILViewer2D
#from FindCenterOfRotation import *
import vtk

import os




def getEntry(location):
    for item in nx[location].keys():
        print (item)


print ("Loading Data")
 


#recon = Reconstructor(algorithm = Algorithm.CGLS, normalized_projection = norm,
#                 angles = angle_proj, center_of_rotation = 86.2 , 
#                 flat_field = flat, dark_field = dark, 
#                 iterations = 15, resolution = 1, isLogScale = False, threads = 3)

#recon = Reconstructor(algorithm = Reconstructor.Algorithm.CGLS, projection_data = proj,
#                 angles = angle_proj, center_of_rotation = 86.2 , 
#                 flat_field = flat, dark_field = dark, 
#                 iterations = 15, resolution = 1, isLogScale = False, threads = 3)
#img_cgls = recon.reconstruct()
#
#pars = dict()
#pars['algorithm'] = Reconstructor.Algorithm.SIRT
#pars['projection_data'] = proj
#pars['angles'] = angle_proj
#pars['center_of_rotation'] = numpy.double(86.2)
#pars['flat_field'] = flat
#pars['iterations'] = 15
#pars['dark_field'] = dark
#pars['resolution'] = 1
#pars['isLogScale'] = False
#pars['threads'] = 3
#
#img_sirt = recon.reconstruct(pars)
#
#recon.pars['algorithm'] = Reconstructor.Algorithm.MLEM
#img_mlem = recon.reconstruct()

############################################################


tot_images = 2203

## create a reduced 3D stack
nreduced = 360*2
directory = r'D:\Documents\Dataset\Chadwick_Flange_Tomo'

indices = [int(float(num) / float(nreduced) * float(tot_images)) for num in range(nreduced)]
angle_proj = numpy.asarray(indices) * 180. / tot_images




#print ("VOI Dimensions {0}".format(str(voi.GetOutput().GetDimensions())))

# center of rotation 
cor = 508.25
############################################################
#recon.pars['algorithm'] = Reconstructor.Algorithm.CGLS_CONV
#recon.pars['regularize'] = numpy.double(0.1)
#img_cgls_conv = recon.reconstruct()

niterations = 10
threads = 2
recon = Reconstructor(algorithm = Reconstructor.Algorithm.CGLS,
                      center_of_rotation = 508.25 , 
                      iterations = niterations,
                      resolution = 1,
                      isLogScale = False,
                      threads = threads)

##data = (norm[:6].T).copy()
###del norm
##del reader
##

use_proj = [360 * i for i in range(1,5)]

for nproj in [360 * 2]:
    #load data file
    
    indices = [int(float(num) / float(nproj) * float(tot_images)) for num in range(nproj)]
    angle_proj = numpy.asarray(indices) * 180. / tot_images

    norm = numpy.load(r'D:\Documents\Dataset\Chadwick_Flange_Tomo\Chadwick_flange_norm_%d.npy' % nproj)
    print ("Selecting ROI")
    data = norm[:,0:1,:].copy()
    del norm
    print ("Normalizing")
    pp = ProjectionPreprocessor()
    pp.normalizedProjections = data
    data = pp.clipNormalizedProjections()
    print ("Launching CGLS 60-61 {0:02} data proj {1}".format(niterations, nproj))
    
    recon.setParameter(normalized_projection_data=data,
                                 angles=angle_proj)
    img_cgls = recon.reconstruct()

    numpy.save("chadwick_chamber60-61_{1:03}_{0:02}.npy".format(
                niterations , nproj),
               img_cgls)

#v = CILViewer2D()
#v.setInputAsNumpy(img_cgls)


##print ("Launching CGLS 970-990")
##extent[2] = 970
##extent[3] = 990
##voi.SetVOI(extent)
##voi.Update()
##norm = Converter.vtk2numpy(voi.GetOutput())
##pp.normalizedProjections = norm
##norm = pp.getNormalizedProjections()
##
####data = (norm[3:6].T).copy()
##img_cgls = alg.cgls(norm.T.copy(), angle_proj, numpy.double(cor), 1 , niterations, threads, False)
##numpy.save("chadwick_chamber970-990.npy", img_cgls)
##print ("Launching CGLS 990-1019")
##extent[2] = 990
##extent[3] = 1019
##voi.SetVOI(extent)
##voi.Update()
##norm = Converter.vtk2numpy(voi.GetOutput())
##pp.normalizedProjections = norm
##norm = pp.getNormalizedProjections()
##
####data = (norm[3:6].T).copy()
##img_cgls = alg.cgls(norm.T.copy(), angle_proj, numpy.double(cor), 1 , niterations, threads, False)
##numpy.save("chadwick_chamber990-1019.npy", img_cgls)
##print ("Launching CGLS 0-6")

####data = (norm[:6].T).copy()
##
##extent[4] = 0
##extent[5] = 5
##voi.SetVOI(extent)
##voi.Update()
##norm = Converter.vtk2numpy(voi.GetOutput())
##pp = ProjectionPreprocessor()
##pp.normalizedProjections = norm
##norm = pp.getNormalizedProjections()
##
##img_cgls = alg.cgls(norm.T.copy(), angle_proj, numpy.double(cor), 1 , niterations, threads, False)
##numpy.save("chadwick_chamber0-6.npy", img_cgls)
##
##extent[4] = 3
##extent[5] = 5
##voi.SetVOI(extent)
##voi.Update()
##norm = Converter.vtk2numpy(voi.GetOutput())
##pp = ProjectionPreprocessor()
##pp.normalizedProjections = norm
##norm = pp.getNormalizedProjections()
##
##img_cgls = alg.cgls(norm.T.copy(), angle_proj, numpy.double(cor), 1 , niterations, threads, False)
##numpy.save("chadwick_chamber3-6.npy", img_cgls)
##img_mlem = alg.mlem(norm, angle_proj, numpy.double(cor), 1 , niterations, threads, False)
##img_sirt = alg.sirt(norm, angle_proj, numpy.double(cor), 1 , niterations, threads, False)
##
##iteration_values = numpy.zeros((niterations,))
##img_cgls_conv = alg.cgls_conv(norm, angle_proj, numpy.double(cor), 1 , niterations, threads,
##                              iteration_values, False)
##print ("iteration values %s" % str(iteration_values))
##
##iteration_values = numpy.zeros((niterations,))
##img_cgls_tikhonov = alg.cgls_tikhonov(norm, angle_proj, numpy.double(cor), 1 , niterations, threads,
##                                      numpy.double(1e-5), iteration_values , False)
##print ("iteration values %s" % str(iteration_values))
##iteration_values = numpy.zeros((niterations,))
##img_cgls_TVreg = alg.cgls_TVreg(norm, angle_proj, numpy.double(cor), 1 , niterations, threads,
##                                      numpy.double(1e-5), iteration_values , False)
##print ("iteration values %s" % str(iteration_values))
##
##
####numpy.save("cgls_recon.npy", img_data)
##import matplotlib.pyplot as plt
##fig, ax = plt.subplots(1,6,sharey=True)
##ax[0].imshow(img_cgls[80])
##ax[0].axis('off')  # clear x- and y-axes
##ax[1].imshow(img_sirt[80])
##ax[1].axis('off')  # clear x- and y-axes
##ax[2].imshow(img_mlem[80])
##ax[2].axis('off')  # clear x- and y-axesplt.show()
##ax[3].imshow(img_cgls_conv[80])
##ax[3].axis('off')  # clear x- and y-axesplt.show()
##ax[4].imshow(img_cgls_tikhonov[80])
##ax[4].axis('off')  # clear x- and y-axesplt.show()
##ax[5].imshow(img_cgls_TVreg[80])
##ax[5].axis('off')  # clear x- and y-axesplt.show()
##
##
###plt.show()
##fig.savefig('recon3d.png')
##plt.close(fig)
##
