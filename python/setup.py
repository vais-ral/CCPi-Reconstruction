#!/usr/bin/env python

import setuptools
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import os
import numpy
import platform	

extra_compile_args = ['-fopenmp','-O2', '-funsigned-char', '-Wall', '-Werror']
if platform.system() == 'Windows':
   extra_compile_args[0:] = ['/DWIN32','/EHsc']
setup(
    name='ccpi',
	description='This is a CCPi package for Tomography',
	version='0.1',
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("ccpi.reconstruction",
                             sources=["ccpi_reconstruction.cpp", 
										"mpi.cpp", 
										"utils.cpp",
										"instruments.cpp",
										"results.cpp",
										"voxels.cpp",
										"cgls.cpp",
										"mlem.cpp",
										"sirt.cpp",
										"total_v.cpp",
										"parallel.cpp",
										"cone.cpp",
										"diamond.cpp",
										"xtek.cpp",
										"tv_reg.cpp",
										"tv_core.cpp",
										"p2D.cpp",
										"c2D.cpp",
										"tiff.cpp",
										"timer.cpp",
										"tikhonov.cpp",
										"ui_calls.cpp"],
                             include_dirs=[numpy.get_include(), ], extra_compile_args=extra_compile_args), 
				  Extension("ccpi.filters",
                             sources=["ccpi_filters.cpp", 
										"mpi.cpp", 
										"utils.cpp",
										"instruments.cpp",
										"results.cpp",
										"voxels.cpp",
										"cgls.cpp",
										"mlem.cpp",
										"sirt.cpp",
										"total_v.cpp",
										"parallel.cpp",
										"cone.cpp",
										"diamond.cpp",
										"xtek.cpp",
										"tv_reg.cpp",
										"tv_core.cpp",
										"p2D.cpp",
										"c2D.cpp",
										"tiff.cpp",
										"timer.cpp",
										"tikhonov.cpp",
										"ui_calls.cpp"],
                             include_dirs=[numpy.get_include(), ], extra_compile_args=extra_compile_args)							 ],
	packages = {'ccpi'}
)
