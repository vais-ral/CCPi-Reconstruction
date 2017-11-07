% Example script to reconstruct in 2D using CGLS

% N. Wadeson

% 06/06/2012

% Solving the algebraic problem: Ax = b.
%
% Parameter Definitions:
% Pixels_nY:	The number of detector pixel columns (max 2000)
% Pixels_nZ:    The number of detector pixel rows (max 2000)
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
%		value will be overwritten by the number in the file with
%		voxel(3) set to 1.
%
% output:
% cglsOut: 	The solution vector for iterations as above.
% rho: 		The residual norm
% eta:		The solution norm
%
% ***NOTE*** 
% (1) Check which version of load_data you wish to use (below)
% (2) Do you wish to use the beam hardening correction?

%------------------Parameters--------------------------------------------
pixels_nY = 2000; 
pixels_nZ = 2000; 
iterations = 10; 
voxels = [2000 2000 1];
%------------------------------------------------------------------------

addpath c/
addpath tools/

% Load data and geometrical parameters from file
% Use this if you are using a reconstruction.xtekct file from FDK reconstruction
%[data geom] = load_data(pathname, filename, pathname2, filename2);
% Use this otherwise
[data geom] = load_data(pathname, filename);

% Cut down detector region if only reconstructing a sub-section
aY = (2000 - pixels_nY)/2;
aZ = (2000 - pixels_nZ)/2;
geom.dets.z = geom.dets.z(aZ+1:aZ+pixels_nZ);
geom.dets.nz = length(geom.dets.z);
geom.dets.y = geom.dets.y(aY+1:aY+pixels_nY);
geom.dets.ny = length(geom.dets.y);

% Option to overwrite the original data
% data = foam_minus_casing_100Projs;

% Cut down data to only include the chosen detectors
data = data(aY+1:aY+pixels_nY, aZ+1:aZ+pixels_nZ,:);

% Convert data to 2D: central slice only
[b geom] = convert2D(data, geom);
% Amend detector z position to ensure rays are 'in plane'
geom.dets.z = 0.0;
b = b(:);

% Optional beam hardening correction
b = b.^2;

% Find centre of rotation
geom = centre_geom(b, geom);

% Call cgls main function
[cglsOut rho eta] = cgls_XTek_single(b, iterations, geom, voxels);

% Set all negative values to zero
cglsOut(find(cglsOut < 0)) = 0;
