import numpy as np
import h5py
from ccpi.reconstruction.parallelbeam import alg as ccpi_reconstruction
angles=np.load('angles.npy')
sino=np.load('sino.npy')
sino = sino.astype(np.float32)
angles = angles.astype(np.float32)
voxels = ccpi_reconstruction.cgls(sino, angles, np.double(0.0), 1, 5, 1, False) 
#output = h5py.File('output.h5','w')
#odata = output.create_dataset("default", data=voxels)
#output.close()
#output = h5py.File('input.h5','w')
#odata = output.create_dataset("default", data=sino)
#output.close()