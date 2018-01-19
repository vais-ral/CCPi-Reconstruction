## Creating test data for Reconstruction ParaView plugin

The process for setting up a test case takes a couple of steps, due to the fact that 
the example datafile is unable to be read properly by the NexusReader plugin (it cannot 
read the angle array property). Two python scripts are needed, which are in the respository 
in the Wrappers/Paraview/test directory, as well as the original example file found at 
https://github.com/DiamondLightSource/Savu/blob/master/test_data/data/24737_fd.nxs.

First, you need to use preprocess_data.py. This takes the nxs file, normalises it using 
the darks and flats, and finally saves the normalised projected data in an mha file. It also 
prints out an array of numbers, these are the angles associated with the projected data. Copy 
this array into the next file.

Second, you need angles_creator.py. We don't actually run this, instead we use this as 
a Programmable Filter in ParaView. Make sure that angle_arr is the same as the print out 
you got from preprocess_data.py. Have this file open so you can copy its contents into 
ParaView.

The steps for loading the test data are below:

1. Open ImageData.mha
2. Apply a Programmable Filter to ImageData.mha
3. Copy and paste code from angles\_creator.py, this creates our angle array as a vtkImageData object
4. Apply Reconstruction(CCPi) filter, with ImageData.mha as pixels and the Programmable Filter 
as the angles
5. Set the iterations to something like 10, and the rotation centre to 86.2
6. Press apply