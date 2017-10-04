from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
from scipy import ndimage
#import logging



#logger = logging.getLogger(__name__)
#
#
#__author__ = "Doga Gursoy, Luis Barroso-Luque, Nghia Vo"
#__copyright__ = "Copyright (c) 2015, UChicago Argonne, LLC."
#__docformat__ = 'restructuredtext en'
#__all__ = ['find_center_vo',
#           ]
PI = 3.14159265359


import h5py

#!/usr/bin/env python
# -*- coding: utf-8 -*-

# #########################################################################
# Copyright (c) 2015, UChicago Argonne, LLC. All rights reserved.         #
#                                                                         #
# Copyright 2015. UChicago Argonne, LLC. This software was produced       #
# under U.S. Government contract DE-AC02-06CH11357 for Argonne National   #
# Laboratory (ANL), which is operated by UChicago Argonne, LLC for the    #
# U.S. Department of Energy. The U.S. Government has rights to use,       #
# reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR    #
# UChicago Argonne, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR        #
# ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is     #
# modified to produce derivative works, such modified software should     #
# be clearly marked, so as not to confuse it with the version available   #
# from ANL.                                                               #
#                                                                         #
# Additionally, redistribution and use in source and binary forms, with   #
# or without modification, are permitted provided that the following      #
# conditions are met:                                                     #
#                                                                         #
#     * Redistributions of source code must retain the above copyright    #
#       notice, this list of conditions and the following disclaimer.     #
#                                                                         #
#     * Redistributions in binary form must reproduce the above copyright #
#       notice, this list of conditions and the following disclaimer in   #
#       the documentation and/or other materials provided with the        #
#       distribution.                                                     #
#                                                                         #
#     * Neither the name of UChicago Argonne, LLC, Argonne National       #
#       Laboratory, ANL, the U.S. Government, nor the names of its        #
#       contributors may be used to endorse or promote products derived   #
#       from this software without specific prior written permission.     #
#                                                                         #
# THIS SOFTWARE IS PROVIDED BY UChicago Argonne, LLC AND CONTRIBUTORS     #
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT       #
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS       #
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL UChicago     #
# Argonne, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,        #
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,    #
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;        #
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER        #
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT      #
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN       #
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE         #
# POSSIBILITY OF SUCH DAMAGE.                                             #
# #########################################################################


def as_ndarray(arr, dtype=None, copy=False):
    if not isinstance(arr, np.ndarray):
        arr = np.array(arr, dtype=dtype, copy=copy)
    return arr


def as_dtype(arr, dtype, copy=False):
    if not arr.dtype == dtype:
        arr = np.array(arr, dtype=dtype, copy=copy)
    return arr


def as_float32(arr):
    arr = as_ndarray(arr, np.float32)
    return as_dtype(arr, np.float32)





def find_center_vo(tomo, ind=None, smin=-40, smax=40, srad=10, step=0.5,
                   ratio=2., drop=20):
    """
    Find rotation axis location using Nghia Vo's method. :cite:`Vo:14`.

    Parameters
    ----------
    tomo : ndarray
        3D tomographic data.
    ind : int, optional
        Index of the slice to be used for reconstruction.
    smin, smax : int, optional
        Reference to the horizontal center of the sinogram.
    srad : float, optional
        Fine search radius.
    step : float, optional
        Step of fine searching.
    ratio : float, optional
        The ratio between the FOV of the camera and the size of object.
        It's used to generate the mask.
    drop : int, optional
        Drop lines around vertical center of the mask.

    Returns
    -------
    float
        Rotation axis location.
        
    Notes
    -----
    The function may not yield a correct estimate, if:
    
    - the sample size is bigger than the field of view of the camera. 
      In this case the ``ratio`` argument need to be set larger
      than the default of 2.0.
    
    - there is distortion in the imaging hardware. If there's 
      no correction applied, the center of the projection image may 
      yield a better estimate.
    
    - the sample contrast is weak. Paganin's filter need to be applied 
      to overcome this. 
   
    - the sample was changed during the scan. 
    """
    tomo = as_float32(tomo)

    if ind is None:
        ind = tomo.shape[1] // 2
    _tomo = tomo[:, ind, :]

    

    # Reduce noise by smooth filters. Use different filters for coarse and fine search 
    _tomo_cs = ndimage.filters.gaussian_filter(_tomo, (3, 1))
    _tomo_fs = ndimage.filters.median_filter(_tomo, (2, 2))

    # Coarse and fine searches for finding the rotation center.
    if _tomo.shape[0] * _tomo.shape[1] > 4e6:  # If data is large (>2kx2k)
        #_tomo_coarse = downsample(np.expand_dims(_tomo_cs,1), level=2)[:, 0, :]
        #init_cen = _search_coarse(_tomo_coarse, smin, smax, ratio, drop)
        #fine_cen = _search_fine(_tomo_fs, srad, step, init_cen*4, ratio, drop)
        init_cen = _search_coarse(_tomo_cs, smin, smax, ratio, drop)
        fine_cen = _search_fine(_tomo_fs, srad, step, init_cen, ratio, drop)
    else:
        init_cen = _search_coarse(_tomo_cs, smin, smax, ratio, drop)
        fine_cen = _search_fine(_tomo_fs, srad, step, init_cen, ratio, drop)

    #logger.debug('Rotation center search finished: %i', fine_cen)
    return fine_cen



def _search_coarse(sino, smin, smax, ratio, drop):
    """
    Coarse search for finding the rotation center.
    """
    (Nrow, Ncol) = sino.shape
    centerfliplr = (Ncol - 1.0) / 2.0

    # Copy the sinogram and flip left right, the purpose is to
    # make a full [0;2Pi] sinogram
    _copy_sino = np.fliplr(sino[1:])

    # This image is used for compensating the shift of sinogram 2
    temp_img = np.zeros((Nrow - 1, Ncol), dtype='float32')
    temp_img[:] = sino[-1]

    # Start coarse search in which the shift step is 1
    listshift = np.arange(smin, smax + 1)
    listmetric = np.zeros(len(listshift), dtype='float32')
    mask = _create_mask(2 * Nrow - 1, Ncol, 0.5 * ratio * Ncol, drop)
    for i in listshift:
        _sino = np.roll(_copy_sino, i, axis=1)
        if i >= 0:
            _sino[:, 0:i] = temp_img[:, 0:i]
        else:
            _sino[:, i:] = temp_img[:, i:]
        listmetric[i - smin] = np.sum(np.abs(np.fft.fftshift(
            #pyfftw.interfaces.numpy_fft.fft2(
            #    np.vstack((sino, _sino)))
            np.fft.fft2(np.vstack((sino, _sino)))
            )) * mask)
    minpos = np.argmin(listmetric)
    return centerfliplr + listshift[minpos] / 2.0


def _search_fine(sino, srad, step, init_cen, ratio, drop):
    """
    Fine search for finding the rotation center.
    """
    Nrow, Ncol = sino.shape
    centerfliplr = (Ncol + 1.0) / 2.0 - 1.0
    # Use to shift the sinogram 2 to the raw CoR.
    shiftsino = np.int16(2 * (init_cen - centerfliplr))
    _copy_sino = np.roll(np.fliplr(sino[1:]), shiftsino, axis=1)
    if init_cen <= centerfliplr:
        lefttake = np.int16(np.ceil(srad + 1))
        righttake = np.int16(np.floor(2 * init_cen - srad - 1))
    else:
        lefttake = np.int16(np.ceil(
            init_cen - (Ncol - 1 - init_cen) + srad + 1))
        righttake = np.int16(np.floor(Ncol - 1 - srad - 1))
    Ncol1 = righttake - lefttake + 1
    mask = _create_mask(2 * Nrow - 1, Ncol1, 0.5 * ratio * Ncol, drop)
    numshift = np.int16((2 * srad) / step) + 1
    listshift = np.linspace(-srad, srad, num=numshift)
    listmetric = np.zeros(len(listshift), dtype='float32')
    factor1 = np.mean(sino[-1, lefttake:righttake])
    num1 = 0
    for i in listshift:
        _sino = ndimage.interpolation.shift(
            _copy_sino, (0, i), prefilter=False)
        factor2 = np.mean(_sino[0,lefttake:righttake])
        _sino = _sino * factor1 / factor2
        sinojoin = np.vstack((sino, _sino))
        listmetric[num1] = np.sum(np.abs(np.fft.fftshift(
            #pyfftw.interfaces.numpy_fft.fft2(
            #    sinojoin[:, lefttake:righttake + 1])
            np.fft.fft2(sinojoin[:, lefttake:righttake + 1])
            )) * mask)
        num1 = num1 + 1
    minpos = np.argmin(listmetric)
    return init_cen + listshift[minpos] / 2.0


def _create_mask(nrow, ncol, radius, drop):
    du = 1.0 / ncol
    dv = (nrow - 1.0) / (nrow * 2.0 * PI)
    centerrow = np.ceil(nrow / 2) - 1
    centercol = np.ceil(ncol / 2) - 1
    # added by Edoardo Pasca
    centerrow = int(centerrow)
    centercol = int(centercol)
    mask = np.zeros((nrow, ncol), dtype='float32')
    for i in range(nrow):
        num1 = np.round(((i - centerrow) * dv / radius) / du)
        (p1, p2) = np.int16(np.clip(np.sort(
            (-num1 + centercol, num1 + centercol)), 0, ncol - 1))
        mask[i, p1:p2 + 1] = np.ones(p2 - p1 + 1, dtype='float32')
    if drop < centerrow:
        mask[centerrow - drop:centerrow + drop + 1,
             :] = np.zeros((2 * drop + 1, ncol), dtype='float32')
    mask[:,centercol-1:centercol+2] = np.zeros((nrow, 3), dtype='float32')
    return mask



#fname = "C:\\Users\\ofn77899\\Documents\\CCPi\\CGLS\\24737_fd_2.nxs"
#nx = h5py.File(fname, "r")
#
## the data are stored in a particular location in the hdf5
#for item in nx['entry1/tomo_entry/data'].keys():
#    print (item)
#
#data = nx.get('entry1/tomo_entry/data/rotation_angle')
#angles = np.zeros(data.shape)
#data.read_direct(angles)
#print (angles)
## angles should be in degrees
#
#data = nx.get('entry1/tomo_entry/data/data')
#stack = np.zeros(data.shape)
#data.read_direct(stack)
#print (data.shape)
#
#print ("Data Loaded")
#
#
## Normalize
#data = nx.get('entry1/tomo_entry/instrument/detector/image_key')
#itype = np.zeros(data.shape)
#data.read_direct(itype)
## 2 is dark field
#darks = [stack[i] for i in range(len(itype)) if itype[i] == 2 ]
#dark = darks[0]
#for i in range(1, len(darks)):
#    dark += darks[i]
#dark = dark / len(darks)
##dark[0][0] = dark[0][1]
#
## 1 is flat field
#flats = [stack[i] for i in range(len(itype)) if itype[i] == 1 ]
#flat = flats[0]
#for i in range(1, len(flats)):
#    flat += flats[i]
#flat = flat / len(flats)
##flat[0][0] = dark[0][1]
#
#
## 0 is projection data
#proj = [stack[i] for i in range(len(itype)) if itype[i] == 0 ]
#angle_proj = [angles[i] for i in range(len(itype)) if itype[i] == 0 ]
#angle_proj = np.asarray (angle_proj)
#angle_proj = angle_proj.astype(np.float32)
#
## normalized data are
## norm = (projection - dark)/(flat-dark)
#
#def normalize(projection, dark, flat, def_val=0.1):
#    a = (projection - dark)
#    b = (flat-dark)
#    with np.errstate(divide='ignore', invalid='ignore'):
#        c = np.true_divide( a, b )
#        c[ ~ np.isfinite( c )] = def_val  # set to not zero if 0/0 
#    return c
#    
#
#norm = [normalize(projection, dark, flat) for projection in proj]
#norm = np.asarray (norm)
#norm = norm.astype(np.float32)
#

#cor = find_center_vo (norm)
#print ("Center of rotation %f" % cor)