#!/usr/bin/env python

import setuptools
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import os
import numpy
import platform	

extra_include_dirs = [numpy.get_include()]
extra_library_dirs = []
extra_compile_args = ['-fopenmp','-O2', '-funsigned-char', '-Wall', '-Werror']
extra_libraries = []
if platform.system() == 'Windows':
   extra_compile_args[0:] = ['/DWIN32','/EHsc']   
   extra_include_dirs += ['C:\\Apps\\Anaconda3\\envs\\cgls\\Library\\include',"..\\src\\","..\\src\\Algorithms","..\\src\\Readers", "."]
   extra_library_dirs += ['C:\Apps\Anaconda3\envs\cgls\Library\lib']
else:
   extra_libraries = ['boost_python','gomp']

setup(
    name='ccpi',
	description='This is a CCPi package for Tomography',
	version='0.1',
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("ccpi.diamond",
                             sources=[  "src/diamond_module.cpp",
                                        "src/diamond_wrapper.cpp",
										"../src/mpi.cpp", 
										"../src/utils.cpp",
										"../src/instruments.cpp",
										"../src/results.cpp",
										"../src/voxels.cpp",
										"../src/Algorithms/cgls.cpp",
										"../src/Algorithms/mlem.cpp",
										"../src/Algorithms/sirt.cpp",
										"../src/total_v.cpp",
										"../src/parallel.cpp",
										"../src/cone.cpp",
										"../src/diamond.cpp",
										"../src/Readers/xtek.cpp",
										"../src/Algorithms/tv_reg.cpp",
										"../src/tv_core.cpp",
										"../src/p2D.cpp",
										"../src/c2D.cpp",
										"../src/Readers/tiff.cpp",
										"../src/timer.cpp",
										"../src/tikhonov.cpp",
										"../src/ui_calls.cpp"],
                             include_dirs=extra_include_dirs, library_dirs=extra_library_dirs, extra_compile_args=extra_compile_args, libraries=extra_libraries ), ],
	packages = {'ccpi'}
)
