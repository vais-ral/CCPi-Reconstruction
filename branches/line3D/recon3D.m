% Example script to reconstruct in 3D using CGLS

% N. Wadeson

% 06/06/2012


% Solving the algebraic problem: Ax = b.
%
% Parameter Definitions:
% Pixels_nY:	The number of detector pixel columns
% Pixels_nZ:    The number of detector pixel rows
%		If only performing a sub-section reconstruction such that
%		not all the rays are required, cut them down here. 
%		The detector array centre remains unchanged. 
%
% iterations:	The number of CGLS iterations to perform.  The solution 
%		vector x is output for each multiple of 10 iterations
%		and/or the maximum iteration.
% 		For example:
%		iterations = 5, x is output at 5i
%		iterations = 15, x is output at 10i and 15i
%		iterations = 26, x is output at 10i, 20i, 26i
%
% voxels:	A 3-vector giving the number of voxels for reconstruction 
%		in x, y and z respectively (with z pointing upwards). 
%		If you are reading in a reconstruction.xtekct file, this
%		value will be overwritten by the number in the file.
%
% output:
% cglsOut: 	The solution vector for iterations as above.
% rho:
% eta:
%
% ***NOTE*** 
% (1) Check which version of load_data you wish to use (below)
% (2) Do you wish to use the beam hardening correction?

%------------------Parameters--------------------------------------------
pixels_nY = 1000;
pixels_nZ = 1000;
iterations = 10;
voxels = [500 500 500];
%------------------------------------------------------------------------

addpath c/
addpath tools/

% Load data and geometrical parameters from file
% Use this if you are using a reconstruction.xtekct file from FDK reconstruction
%[data geom] = load_data(pathname, filename, pathname2, filename2);
% Use this otherwise
[data geom] = load_data(pathname, filename);

% Cut down detector region if only reconstructing a sub-section
%aY = (2000 - pixels_nY)/2;
%aZ = (2000 - pixels_nZ)/2;
%geom.dets.z = geom.dets.z(aZ+1:aZ+pixels_nZ);
%geom.dets.nz = length(geom.dets.z);
%geom.dets.y = geom.dets.y(aY+1:aY+pixels_nY);
%geom.dets.ny = length(geom.dets.y);

% Option to overwrite the original data
% data = foam_minus_casing;

% Cut down data to only include the chosen detectors
%data = data(aY+1:aY+pixels_nY, aZ+1:aZ+pixels_nZ,:);

% Find centre of rotation
geom = centre_geom(data, geom);
b = data(:);

% Optional beam-hardening correction
b = b.^2;

% Call cgls main function
[cglsOut rho eta] = cgls_XTek_single(b, iterations, geom, voxels);
